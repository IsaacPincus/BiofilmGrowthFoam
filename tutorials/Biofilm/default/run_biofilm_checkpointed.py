#%%
"""
Biofilm growth driver with checkpoint / restart support.

Usage
-----
Set RESUME_MODE below:
    "auto"      -> resume from the latest checkpoint if one exists, else start fresh.
    "fresh"     -> always start from scratch (deletes existing checkpoints).
    "restart"   -> force a resume; raises if no checkpoint is found.
    "bootstrap" -> resume an already-run simulation that has NO checkpoint, by
                   rebuilding the Python-side state from the existing OpenFOAM
                   case folders. B is read from the last completed step; Bdead is
                   read if present, otherwise reset to zero; the RNG starts fresh.
                   A normal checkpoint is then written so future crashes resume
                   via the usual path.

How it works
------------
A small checkpoint (pickle) is written at the END of every completed step into
CHECKPOINT_DIR. It stores only the state that does NOT already live on disk in
the OpenFOAM cases (mainly BDeadField, the 2D fields, the latest case times, and
the RNG states). The OpenFOAM time directories themselves are the heavy state and
are left in place.

On restart we:
  1. prune any time directory newer than the checkpoint (partial output from an
     interrupted step),
  2. rebuild the deterministic mesh<->image mapping from the existing mesh,
  3. restore the Python-side arrays and RNG, and
  4. continue the time loop from the next step.

Checkpoints are written atomically (temp file + os.replace), so a crash during a
write can never leave a corrupted checkpoint behind.
"""

import numpy as np
import foamlib
import pythonScripts.utils as utils
import os
import fluidfoam
from PIL import Image
import shutil
import pickle
import time

# ============================================================
# Restart / checkpoint configuration
# ============================================================
RESUME_MODE = "auto"          # "auto" | "fresh" | "restart" | "bootstrap"
CHECKPOINT_DIR = "./checkpoints"

# Bootstrap only: if True, always reset Bdead to zero instead of trying to read it
# from the existing case folder.
BOOTSTRAP_FORCE_ZERO_BDEAD = False

# ============================================================
# Parallel (MPI) configuration
# ============================================================
# When True, the two OpenFOAM solver calls are run with MPI. The biofilm
# growth/division (numpy) stays serial; the Python loop reconstructs the fields
# back to a single domain after each solve so the image<->cell mapping still
# works. Set False for the original serial behaviour.
#
# IMPORTANT: cpus MUST equal numberOfSubdomains in each case's
# system/decomposeParDict. This script (re)writes that dict to match N_CPUS, so
# just set N_CPUS here. Your custom solvers must be parallel-safe (they are, as
# long as they only use standard fvm/fvc operators and per-cell local logic).
RUN_PARALLEL = True
N_CPUS = 2
DECOMP_METHOD = "scotch"      # "scotch" (no coeffs needed) or "simple", etc.

# ============================================================
# Adaptive yield-stress sub-stepping
# ============================================================
# Instead of one fixed-duration yield solve, relax in small chunks and stop early
# once the biofilm has stopped moving. After each chunk we look at nu2 in the
# biofilm-only region (alpha.water < alpha_biofilm_cut) and decide whether to keep
# going. See run_yield_adaptive() for the stop condition (and how to invert it).
YIELD_ADAPTIVE = True
dt_yield_chunk = 0.1          # seconds of relaxation simulated per chunk
t_yield_max = 0.5             # hard cap on total relaxation per step
nu2_threshold = 3e-3          # convergence threshold on nu2 in the biofilm region
alpha_biofilm_cut = 0.05      # "biofilm only" == alpha.water below this


def _time_dirs(case_path):
    """Return sorted [(time, dirname)] for numeric time directories in a case."""
    out = []
    for name in os.listdir(case_path):
        if os.path.isdir(os.path.join(case_path, name)):
            try:
                out.append((float(name), name))
            except ValueError:
                pass
    return sorted(out)


def prune_time_dirs_after(case_path, max_time, tol=1e-9):
    """Delete any time directory whose time is greater than max_time.

    Used on restart to remove partial output written by an interrupted step."""
    removed = []
    for t, name in _time_dirs(case_path):
        if t > max_time + tol:
            shutil.rmtree(os.path.join(case_path, name))
            removed.append(name)
    if removed:
        print(f"  pruned stray time dirs in {case_path}: {removed}")
    return removed


def detect_last_step_from_case(case_path, step_time, tol=1e-6):
    """Find the largest completed step N for which a folder named N*step_time exists.

    Used by 'bootstrap' mode to figure out how far an un-checkpointed run got.
    The canonical per-step folders are integer multiples of step_time (= dt); the
    small intermediate solver folders (dt*k + 0.25, + 0.5, ...) are ignored."""
    best = 0
    for t, _name in _time_dirs(case_path):
        if t <= 0:
            continue
        k = t / step_time
        if abs(k - round(k)) < tol and round(k) >= 1:
            best = max(best, int(round(k)))
    return best


def latest_checkpoint_path():
    if not os.path.isdir(CHECKPOINT_DIR):
        return None
    cks = [f for f in os.listdir(CHECKPOINT_DIR)
           if f.startswith("checkpoint_") and f.endswith(".pkl")]
    if not cks:
        return None
    cks.sort(key=lambda f: int(f[len("checkpoint_"):-len(".pkl")]))
    return os.path.join(CHECKPOINT_DIR, cks[-1])


def save_checkpoint(step, **state):
    """Atomically write a checkpoint for the just-completed step."""
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    path = os.path.join(CHECKPOINT_DIR, f"checkpoint_{step:05d}.pkl")
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        pickle.dump({"step": step, **state}, fh)
    os.replace(tmp, path)  # atomic swap: a crash mid-write cannot corrupt the file
    print(f"  checkpoint saved: {path}")
    return path


def load_checkpoint(path):
    with open(path, "rb") as fh:
        return pickle.load(fh)


# ============================================================
# Parallel / decomposition helpers
# ============================================================
_DECOMPOSE_DICT_TEMPLATE = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM                                        |
|  \\\\    /   O peration     |                                                 |
|   \\\\  /    A nd           |                                                 |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains {n};
method          {method};

// ************************************************************************* //
"""


def write_decompose_dict(case_path, n, method):
    """(Re)write system/decomposeParDict so numberOfSubdomains matches N_CPUS."""
    path = os.path.join(case_path, "system", "decomposeParDict")
    with open(path, "w") as fh:
        fh.write(_DECOMPOSE_DICT_TEMPLATE.format(n=n, method=method))


def remove_processor_dirs(case_path):
    """Delete all processorN directories (start each run from a clean decomposition)."""
    for name in os.listdir(case_path):
        if name.startswith("processor") and os.path.isdir(os.path.join(case_path, name)):
            shutil.rmtree(os.path.join(case_path, name))


def prepare_parallel(case, case_path, resume_time, n, method):
    """Build a fresh decomposition for `case` from its current serial state.

    Wipes any stale processor dirs, writes the decomposeParDict, decomposes the
    mesh, and guarantees the serial fields at `resume_time` are present on the
    processors (so the first parallel solve starts from the right state)."""
    remove_processor_dirs(case_path)
    write_decompose_dict(case_path, n, method)
    case.run(["decomposePar", "-force"])                       # mesh (+ start fields)
    case.run(["decomposePar", "-fields", "-time", str(resume_time), "-force"])


def push_fields(case, time):
    """Decompose the serial fields at `time` onto the (already meshed) processors."""
    case.run(["decomposePar", "-fields", "-time", str(time), "-force"])


def pull_latest(case):
    """Reconstruct only the newest processor time back to a serial time directory."""
    case.run(["reconstructPar", "-latestTime"])


def run_yield_adaptive(case, t_start, chunk, t_max, nu2_thresh, alpha_cut,
                       parallel, cpus):
    """Relax the yield-stress case from t_start in `chunk`-sized increments, up to
    t_max total, stopping early once the biofilm has settled.

    After every chunk we read nu2 in the biofilm-only region (alpha.water <
    alpha_cut). Returns the total relaxation time actually simulated.

    STOP CONDITION (the one judgement call here):
      We assume small nu2 == still deforming, large nu2 == settled. So we keep
      stepping while nu2 is still below the threshold *somewhere* in the biofilm
      (min over the region < thresh), and stop once it is above the threshold
      *everywhere* (the region has gone static). If your nu2 convention is the
      opposite, flip the comparison on the line marked below (and/or swap .min()
      for .max())."""
    t = 0.0
    k = 0
    while t < t_max - 1e-9:
        k += 1
        t_next = min(chunk * k, t_max)
        with case.control_dict as f:
            f["endTime"] = t_start + t_next
        if parallel:
            case.run("biofilmYieldFoam", parallel=True, cpus=cpus)
            case.run(["reconstructPar", "-latestTime"])  # need serial fields to check
        else:
            case.run("biofilmYieldFoam")
        t = t_next

        try:
            with case[-1]["nu2"] as field:
                nu2 = np.asarray(field.internal_field, dtype=float)
        except Exception as e:
            print(f"    [yield] WARNING: could not read nu2 ({e!r}); "
                  f"relaxing the full {t_max:.2f}s this step")
            continue

        with case[-1]["alpha.water"] as field:
            alpha = np.asarray(field.internal_field, dtype=float)

        region = alpha < alpha_cut
        if not np.any(region):
            print(f"    [yield] t={t:.2f}s: no biofilm-only cells yet")
            continue

        agg = nu2[region].min()
        print(f"    [yield] t={t:.2f}s  min nu2 (biofilm) = {agg:.3e}")
        if agg >= nu2_thresh:          # <-- stop condition: flip if nu2 is inverted
            print(f"    [yield] settled after {t:.2f}s")
            break

    return t


# we're going to create an array and use it as the basis for the calculations
# We have the following, all the same size:
#   mask, 1s for free regions, 0s for grains
#   B, biofilm concentration
#   Bdead, 'dead' biofilm concentration

## constants, setting up initial conditions etc
pathConvection = './convection'
pathYieldStress = './yieldStress'

show_plots = False
rng = np.random.default_rng()

# we read in the mask as an image, so it can be whatever we like
mask, openfoam_image_field, Nx, Ny = utils.image_to_openfoam_mask(image_path='circle_grid.png',
                                                                   method='mean', invert=False)

# true domain size
# Nx = 360
# Ny = 90
L = 700e-6
W = Ny/Nx * L
dx = L / Nx
dy = W / Ny
thickness = 1e-5

# yield stresses, [N/m^2]
tauAlive = 4e-3
tauDead = 4e-4
# # detachment rate constant, [m/s]
# kdet = 1/36000

mu = 6e-4       # maximum growth rate for biofilm, [1/s]
Ks = 0.1           # Monod half-saturation constant, [mol/m^3]
Yield = 0.3       # Yield of biomass
kd = 2e-6          # biomass decay constant, [1/s]
# kd = 0
Diffusivity = 1e-8

# # timestep for each removal/growth step, in seconds
dt = 1800
no_steps = 100

dt_yield = 5e-1
dt_convection = 3

# # new bits of biofilm each timestep
# n_new = 5
n_init = 80

Bcut = 8000

# U_threshold = 1e-5

tau0default = 4e-3


# %%
# ------------------------------------------------------------
# Decide fresh vs restart
# ------------------------------------------------------------
_ckpt_path = latest_checkpoint_path()

bootstrap = False
if RESUME_MODE == "restart":
    if _ckpt_path is None:
        raise FileNotFoundError(f"No checkpoint found in {CHECKPOINT_DIR}")
    restart = True
elif RESUME_MODE == "fresh":
    restart = False
    if os.path.isdir(CHECKPOINT_DIR):
        shutil.rmtree(CHECKPOINT_DIR)   # avoid a stale checkpoint being picked up later
    _ckpt_path = None
elif RESUME_MODE == "bootstrap":
    # resume an already-run simulation that has NO checkpoint, reconstructing the
    # Python-side state from the existing OpenFOAM case folders.
    restart = True
    bootstrap = True
    _ckpt_path = None
else:  # "auto"
    restart = _ckpt_path is not None

print("=" * 60)
if bootstrap:
    print("BOOTSTRAPPING from existing case folders (no checkpoint)")
elif restart:
    print("RESTARTING from checkpoint")
else:
    print("Starting FRESH")
print("=" * 60)

ckpt = None
last_step = None
if bootstrap:
    last_step = detect_last_step_from_case(pathYieldStress, dt)
    if last_step < 1:
        raise RuntimeError(
            f"No completed 'dt*step' time folders found in {pathYieldStress} "
            f"(looking for multiples of dt={dt}). Nothing to bootstrap from."
        )
    print(f"  detected last completed step {last_step} "
          f"(yield time {dt*last_step}, convection time {dt_convection*last_step})")
    # drop any partial output above the last completed step before foamlib scans
    prune_time_dirs_after(pathYieldStress, dt * last_step)
    prune_time_dirs_after(pathConvection, dt_convection * last_step)
    # any pre-existing checkpoints are inconsistent with a manual bootstrap
    if os.path.isdir(CHECKPOINT_DIR):
        print(f"  clearing existing {CHECKPOINT_DIR} (rebuilding state from the case)")
        shutil.rmtree(CHECKPOINT_DIR)
elif restart:
    ckpt = load_checkpoint(_ckpt_path)
    print(f"  loaded {_ckpt_path} (completed through step {ckpt['step']})")
    # remove any partial output written by an interrupted step BEFORE foamlib
    # re-scans the case directories.
    prune_time_dirs_after(pathYieldStress, ckpt["yieldstress_time"])
    prune_time_dirs_after(pathConvection, ckpt["convection_time"])


# %%
# Setup yield stress case
caseYieldStress = foamlib.FoamCase(pathYieldStress)  # Loads the OpenFOAM case

if not restart:
    caseYieldStress.clean()

    # creates blockMeshDict
    output = utils.create_blockmesh_dict(Nx, Ny, L, W,
                                         thickness=thickness,
                                         output_dir=os.path.join(pathYieldStress, "system"))
    caseYieldStress.run("blockMesh")  # do blockMesh

    # make sure that we set all fields back to the uniform value
    with caseYieldStress[0]["B"] as field:
        field.internal_field = 0
    with caseYieldStress[0]["Bdead"] as field:
        field.internal_field = 0
    with caseYieldStress[0]["tau0"] as field:
        field.internal_field = 0
    with caseYieldStress[0]["alpha.water"] as field:
        field.internal_field = 0

    # cut out masked regions using topoSet
    nonzeros = np.nonzero(openfoam_image_field)
    cell_id_grains = nonzeros[0]

    utils.create_topoSetDict(cell_id_grains.tolist(),
                             output_path=os.path.join(pathYieldStress, "system/topoSetDict"))

    caseYieldStress.run("topoSet")
    caseYieldStress.run("subsetMesh grainCells -overwrite -patch walls")

    with caseYieldStress.control_dict as f:
        f["endTime"] = dt_yield
        f["writeInterval"] = dt_yield_chunk

# now also setup the other case, including the same grains as in the yield stress case
caseConvection = foamlib.FoamCase(pathConvection)  # Loads the OpenFOAM case

if not restart:
    caseConvection.clean()
    shutil.copytree(os.path.join(pathYieldStress, "constant/polyMesh"),
                    os.path.join(pathConvection, "constant/polyMesh"), dirs_exist_ok=True)
    with caseConvection.transport_properties as f:
        f["Ks"] = Ks
        f["mu"] = mu
        f["Y"] = Yield
        f["D"] = foamlib.Dimensioned(Diffusivity,
                                     foamlib.DimensionSet(length=2, time=-1))

    with caseConvection.control_dict as f:
        f["endTime"] = dt_convection

# %%
# now we need to re-associate each cell with the pixels in the original image.
# we want a list of cells, each with the i,j location of the pixel in the image, which maps to the mask location
# we also want a 2D image with either the cell #, or -1 for a grain location.
# NOTE: this is deterministic given the (existing) mesh + image, so it is rebuilt
# on every run -- fresh or restart -- rather than checkpointed.
x, y, z = fluidfoam.readmesh(pathYieldStress)
# Ensure x, y are NumPy arrays
x = np.asarray(x)
y = np.asarray(y)
# Compute image indices
xn = (x / dx).astype(int)
yn = (y / dy).astype(int)
# Initialize output arrays
image_cell_values = np.array(mask, dtype=int) - 1
cell_locations_in_image = np.stack((xn, yn), axis=1)
# Assign index i to each corresponding (xn, yn)
# Note: This assumes no duplicate (xn, yn) pairs — last one wins if duplicates exist.
image_cell_values[yn, xn] = np.arange(len(x))
img_array = np.asarray(image_cell_values)

# meshgrid (always needed)
# note that meshgrid in numpy is opposite of matlab, it's (ny,nx) array indexing, like an image
X, Y = np.meshgrid(np.linspace(0, L, Nx), np.linspace(0, W, Ny))

# %%
# ------------------------------------------------------------
# Initial conditions (fresh) or restore state (restart)
# ------------------------------------------------------------
if not restart:
    # Parameters
    r0 = W/25
    # x0 = 3*L/5
    # y0 = W/2
    Bvalue = 6000.0
    Balive2D = np.zeros_like(mask, dtype=float)
    Bdead2D = np.zeros_like(mask, dtype=float)
    # # put a patch of alive biofilm at a circle
    # for ii in range(1):
    #     x0 = np.random.uniform(L/5, 4*L/5)
    #     y0 = np.random.uniform(0, W)
    #     Balive2D[(X - x0)**2 + (Y - y0)**2 < r0**2] = Bvalue
    # only next to mask, in a circle
    for ii in range(n_init):
        # random location
        rand_loc = rng.choice(utils.get_boundaries(mask), axis=0)
        x0 = X[rand_loc[0], rand_loc[1]]
        y0 = Y[rand_loc[0], rand_loc[1]]
        # print(f"x0: {x0:.2g}, y0: {y0:.2g}")
        Balive2D[(X - x0)**2 + (Y - y0)**2 < r0**2] = Bvalue
        # Balive2D[rand_loc[0], rand_loc[1]] = Bvalue
    # # put a patch of alive biofilm at the top and bottom
    # condition = ((Y < W/3) | (Y > 2*W/3)) & (X > L/5)
    # Balive2D[condition] = Bvalue

    # tau0 and alpha.water
    tau02D = np.zeros_like(Balive2D)
    tau02D[Balive2D > 0] = tau0default
    alphaWater = np.ones_like(Balive2D)
    alphaWater[Balive2D > 0] = 0

    # add biofilm to locations we want it, and tau0
    BField = utils.to_1d(Balive2D, image_cell_values)
    BDeadField = utils.to_1d(Bdead2D, image_cell_values)
    tau0Field = utils.to_1d(tau02D, image_cell_values)
    alphaWaterField = utils.to_1d(alphaWater, image_cell_values)

    with caseYieldStress[-1]["B"] as field:
        field.internal_field = BField
    with caseYieldStress[-1]["Bdead"] as field:
        field.internal_field = BDeadField
    with caseYieldStress[-1]["tau0"] as field:
        field.internal_field = tau0Field
    with caseYieldStress[-1]["alpha.water"] as field:
        field.internal_field = alphaWaterField

    start_step = 1

    # Save a "step 0" checkpoint so that an interrupted step 1 can be retried.
    save_checkpoint(
        0,
        Balive2D=Balive2D,
        Bdead2D=Bdead2D,
        BDeadField=BDeadField,
        yieldstress_time=float(caseYieldStress[-1].name),
        convection_time=float(caseConvection[-1].name),
        rng_state=rng.bit_generator.state,
        np_random_state=np.random.get_state(),
    )
elif bootstrap:
    # ------------------------------------------------------------
    # Reconstruct Python-side state from the existing case folders.
    #   - B (alive biofilm) is read from the last completed step folder.
    #   - Bdead is read if present, otherwise reset to zero.
    #   - the RNG has no stored state, so a fresh RNG is used.
    # ------------------------------------------------------------
    with caseYieldStress[-1]["B"] as field:
        BField = field.internal_field
    Balive2D = utils.to_2d(BField, image_cell_values)

    BDeadField = None
    if not BOOTSTRAP_FORCE_ZERO_BDEAD:
        try:
            with caseYieldStress[-1]["Bdead"] as field:
                BDeadField = field.internal_field
            print("  Bdead recovered from the case folder")
        except Exception as e:
            print(f"  could not read Bdead from case ({e!r}); resetting to zero")

    if BDeadField is None:
        BDeadField = np.zeros(len(x))
        print("  Bdead reset to zero")
    else:
        BDeadField = np.asarray(BDeadField, dtype=float)
        if BDeadField.ndim == 0:           # uniform field -> expand to all cells
            BDeadField = np.full(len(x), float(BDeadField))
    Bdead2D = utils.to_2d(BDeadField, image_cell_values)

    # No stored RNG state -> continue with a fresh RNG.
    print("  RNG state not available; continuing with a fresh RNG")

    # Re-establish controlDict end times for a clean resume at step start_step.
    with caseYieldStress.control_dict as f:
        f["endTime"] = dt * last_step + dt_yield
        f["writeInterval"] = dt_yield_chunk
    with caseConvection.control_dict as f:
        f["endTime"] = dt_convection * (last_step + 1)

    start_step = last_step + 1
    print(f"  resuming at step {start_step}")

    # Convert the bootstrapped state into a normal checkpoint so any future
    # interruption resumes via the usual checkpoint path.
    save_checkpoint(
        last_step,
        Balive2D=Balive2D,
        Bdead2D=Bdead2D,
        BDeadField=BDeadField,
        yieldstress_time=float(dt * last_step),
        convection_time=float(dt_convection * last_step),
        rng_state=rng.bit_generator.state,
        np_random_state=np.random.get_state(),
    )
else:
    # restore Python-side state from the checkpoint
    Balive2D = ckpt["Balive2D"]
    Bdead2D = ckpt["Bdead2D"]
    BDeadField = ckpt["BDeadField"]
    rng.bit_generator.state = ckpt["rng_state"]
    np.random.set_state(ckpt["np_random_state"])
    start_step = ckpt["step"] + 1
    print(f"  resuming at step {start_step}")


#%%
######################################################################
# Parallel preparation: build a clean decomposition from the current
# serial state (works the same for fresh / restart / bootstrap, since
# each branch above leaves caseX[-1] at the correct resume time).
######################################################################
if RUN_PARALLEL:
    print(f"Preparing parallel run on {N_CPUS} subdomains ({DECOMP_METHOD}) ...")
    prepare_parallel(caseYieldStress, pathYieldStress,
                     caseYieldStress[-1].name, N_CPUS, DECOMP_METHOD)
    prepare_parallel(caseConvection, pathConvection,
                     caseConvection[-1].name, N_CPUS, DECOMP_METHOD)


#%%
######################################################################
# time loop starts here
######################################################################
if start_step >= no_steps:
    print(f"Nothing to do: start_step={start_step} >= no_steps={no_steps}")

for step in range(start_step, no_steps):
    t_start = time.perf_counter()
    print(f"\n=== step {step} / {no_steps - 1} ===")
    # run the yield stress case (processors already hold the start-time fields:
    # from prepare_parallel on the first iteration, or from the end-of-step
    # push_fields() below on subsequent ones)
    if YIELD_ADAPTIVE:
        t_relax = run_yield_adaptive(
            caseYieldStress, float(caseYieldStress[-1].name),
            dt_yield_chunk, t_yield_max, nu2_threshold, alpha_biofilm_cut,
            RUN_PARALLEL, N_CPUS)
        print(f"  yield relaxed for {t_relax:.2f}s")
    elif RUN_PARALLEL:
        caseYieldStress.run("biofilmYieldFoam", parallel=True, cpus=N_CPUS)
        pull_latest(caseYieldStress)        # reconstruct -> serial for the numpy step
    else:
        caseYieldStress.run("biofilmYieldFoam")

    # Copy B, U values to convection case
    with caseYieldStress[-1]["B"] as field:
        BField = field.internal_field
    with caseYieldStress[-1]["Bdead"] as field:
        BdeadField = field.internal_field
    with caseYieldStress[-1]["U"] as field:
        UField = field.internal_field
        UboundField = field.boundary_field

    with caseConvection[-1]["B"] as field:
        field.internal_field = BField
    with caseConvection[-1]["U"] as field:
        field.internal_field = UField
        field.boundary_field = UboundField

    if RUN_PARALLEL:
        # push the just-written serial B, U onto the processors, then solve
        push_fields(caseConvection, caseConvection[-1].name)
        caseConvection.run("convectionWithBiomassFoam", parallel=True, cpus=N_CPUS)
        pull_latest(caseConvection)
    else:
        caseConvection.run("convectionWithBiomassFoam")

    # now we grow the biofilm. For now we won't have any dead biofilm
    with caseConvection[-1]["C"] as field:  # read concentration field
        CField = field.internal_field

    # for each element of the B field, update based on growth
    solidIndicator2D = ( Balive2D > 10 ) - 0.5
    solidIndicator2DMasked = np.ma.masked_array(solidIndicator2D, ~mask)
    # t = skfmm.travel_time(solidIndicator2DMasked, np.ones_like(solidIndicator2D))

    prefactor = (mu*(CField/(Ks + CField)) - kd)
    BdeadField = BdeadField + kd * BField / prefactor * (np.exp(dt*prefactor) - 1)
    BField = BField * np.exp(prefactor * dt)

    Balive2D = utils.to_2d(BField, image_cell_values)
    Bdead2D = utils.to_2d(BdeadField, image_cell_values)

    # if there is any Balive or Bdead outside of alpha.water, we can get rid of it.
    with caseYieldStress[-1]["alpha.water"] as field:
        alphaWaterField = field.internal_field

    alphaWater2D = utils.to_2d(alphaWaterField, image_cell_values)
    # Balive2D[Balive2D < 1] = 0
    # Bdead2D[Bdead2D < 0.1] = 0
    Balive2D[alphaWater2D > 0.5] = 0
    Bdead2D[alphaWater2D > 0.5] = 0
    # Bdead2D[Bdead2D > 0.5] = 0

    # find all cells with values of Balive greater than some cutoff
    # Go through each cell which is growing.
    # If there's an empty cell adjacent, choose one at random and displace half the biomass there.
    # If there's no empty cell, displace that biomass to a random adjacent cell, and restart
    # here, have some vibe coding, thanks Claude
    Bdivide = (Balive2D > Bcut) & (mask)
    while np.any(Bdivide):
        Bdivide = (Balive2D > Bcut) & (mask)
        print(np.sum(Bdivide))
        target_coords = np.column_stack(np.where(Bdivide))

        # Neighbor offsets (4-connectivity)
        offsets = np.array([[-1,0],  [0,-1],
                            [0,1],   [1,0]])

        results = []
        for i, j in target_coords:
            # cut the current value in half
            Balive2D[i,j] = Balive2D[i,j]/2

            # Calculate all neighbor positions
            neighbor_coords = np.array([i, j]) + offsets

            # Filter valid neighbors (in bounds and not masked)
            valid_mask = ((neighbor_coords[:, 0] >= 0) &
                            (neighbor_coords[:, 0] < Balive2D.shape[0]) &
                            (neighbor_coords[:, 1] >= 0) &
                            (neighbor_coords[:, 1] < Balive2D.shape[1]))
            valid_neighbors = neighbor_coords[valid_mask]

            # Remove masked neighbors
            unmasked = mask[valid_neighbors[:, 0], valid_neighbors[:, 1]]
            final_neighbors = valid_neighbors[unmasked]

            # move half the biomass to an adjacent cell
            random_valid_idx = np.random.choice(len(final_neighbors))
            random_coord = final_neighbors[random_valid_idx]
            Balive2D[random_coord[0], random_coord[1]] += Balive2D[i,j]/2

    # create a new folder in YieldStress at the given timestep
    folderTimeName = caseYieldStress[-1].name
    print(folderTimeName)

    # guard against a leftover destination from an interrupted attempt at this step
    dest = './yieldStress/'+str(dt*step)
    if os.path.exists(dest):
        shutil.rmtree(dest)
    shutil.copytree('./yieldStress/'+folderTimeName, dest)

    # yield stress is given by ratio of alive and dead
    # add in a tiny bit of dead everywhere just so this doesn't lead to divide by zero errors
    Bdead2D[Bdead2D < 1e-5] = 1e-5
    # # linear interpolation
    # tau02D = ((Balive2D * tauAlive + Bdead2D * tauDead)/(Balive2D + Bdead2D))
    # harmonic interpolation
    fracA = Balive2D/(Balive2D + Bdead2D);
    tau02D = 1/(fracA/tauAlive + (1-fracA)/tauDead)
    # set alpha.water back to where the biofilm actually is
    alphaWater = np.ones_like(Balive2D)
    alphaWater[Balive2D > 1] = 0

    BField = utils.to_1d(Balive2D, image_cell_values)
    BDeadField = utils.to_1d(Bdead2D, image_cell_values)
    tau0Field = utils.to_1d(tau02D, image_cell_values)
    alphaWaterField = utils.to_1d(alphaWater, image_cell_values)

    with caseYieldStress[-1]["B"] as field:
        field.internal_field = BField
    with caseYieldStress[-1]["Bdead"] as field:
        field.internal_field = BDeadField
    with caseYieldStress[-1]["tau0"] as field:
        field.internal_field = tau0Field
    with caseYieldStress[-1]["alpha.water"] as field:
        field.internal_field = alphaWaterField

    # copy over C from convection as well for comparison purposes
    with caseConvection[-1]["C"] as field:
        CFieldCopy = field.internal_field
    with caseYieldStress[-1]["C"] as field:
        field.internal_field = CFieldCopy
        field.boundary_field = {
            "walls": {"type": "zeroGradient"},
            "inlet": {"type": "zeroGradient"},
            "outlet": {"type": "zeroGradient"},
            "emptyFaces": {"type": "zeroGradient"},
        }

    with caseYieldStress.control_dict as f:
        f["endTime"] = dt*step + dt_yield
    with caseConvection.control_dict as f:
        f["endTime"] = f["endTime"] + dt_convection

    # push the freshly-written serial yield fields onto the processors so the
    # next step's parallel solve starts from them (mesh is unchanged, so this is
    # a cheap -fields decomposition, not a full re-decomposition)
    if RUN_PARALLEL:
        push_fields(caseYieldStress, caseYieldStress[-1].name)

    # ------------------------------------------------------------
    # checkpoint: written LAST, after every field + controlDict update,
    # so disk state and Python state are consistent for "end of this step".
    # ------------------------------------------------------------
    save_checkpoint(
        step,
        Balive2D=Balive2D,
        Bdead2D=Bdead2D,
        BDeadField=BDeadField,
        yieldstress_time=float(caseYieldStress[-1].name),
        convection_time=float(caseConvection[-1].name),
        rng_state=rng.bit_generator.state,
        np_random_state=np.random.get_state(),
    )

    print(f"  step {step} took {time.perf_counter() - t_start:.1f}s")

    # now let's see if it worked!

print("\nDone.")