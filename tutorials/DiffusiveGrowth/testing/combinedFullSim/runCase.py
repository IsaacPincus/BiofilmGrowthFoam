"""
Biofilm growth + detachment driver -- single shared case directory.

Both solvers run against ./biofilm, alternating on separated timescales:

  biofilmYieldFoam   -- VOF with a spatially varying Herschel-Bulkley yield
                        stress. Deforms and detaches the biofilm under flow.
                        Runs for seconds of real time.

  biomassGrowthFoam  -- transport of C, B, A with degenerate nonlinear diffusion
                        + Monod kinetics. Runs for one dt (~30 min).

Sharing a directory removes the entire copy layer: phi never moves, so there is
no dimensions loss, no ASCII truncation, and no time-name mismatch. The two
solvers hand off through the files they both write.

REQUIRES (solver side):
  * Each solver must carry the other's fields as MUST_READ / AUTO_WRITE, even if
    it never touches them -- otherwise the time dir it writes will be missing
    fields the other solver needs to start. Specifically:
        biomassGrowthFoam  must carry  alpha.water, p_rgh, tau0
        biofilmYieldFoam   must carry  C, A
  * controlDict's `application` entry is metadata only; running either binary in
    this case is fine regardless of what it says.
  * fvSchemes needs div() and laplacian() entries for BOTH solvers' fields, since
    `default none` is set. Each solver ignores what it does not use.

B is the biofilm VOLUME FRACTION, so it and alpha.water are complementary:

    B = 1 - alpha.water

The clock is shared and monotonic: each step advances the case by (t_relax + dt).
Since t_relax ~ 1 s and dt = 1800 s, case time ~ biological time to within 0.1%.

State lives on disk, so there is no checkpointing: to resume, set START_STEP to
the last completed step + 1 and delete any time dirs above it.
"""

import os
import shutil
import time

import numpy as np
import foamlib
import fluidfoam

import pythonScripts.utils as utils


# ============================================================
# Configuration
# ============================================================
casePath = "./caseDir"
SNAPSHOT_DIR = "./snapshots"

IMAGE_PATH = "circle_grid.png"
YIELD_SOLVER = "biofilmYieldFoam"
GROWTH_SOLVER = "biomassGrowthFoam"

# --- geometry ---
L = 700e-6                  # domain length [m]; W follows from the image aspect
thickness = 1e-5            # z-extent, single cell layer

# --- time stepping ---
dt = 1800.0                 # growth step [s]
no_steps = 5
START_STEP = 1              # bump this to resume from an existing case

# --- yield (VOF) solve: Courant-limited, needs a small adjustable step ---
YIELD_ADAPTIVE = False
dt_yield = 1.0              # fixed relax duration when YIELD_ADAPTIVE = False
dt_yield_chunk = 0.1        # seconds of relaxation per chunk when adaptive
t_yield_max = 2.0           # hard cap on relaxation per step
nu2_threshold = 3e-3        # convergence threshold on nu2 in the biofilm region
alpha_biofilm_cut = 0.05    # "biofilm only" == alpha.water below this (B > 0.95)

YIELD_CONTROL = {
    "deltaT": 1e-8,
    "adjustTimeStep": True,
    "maxCo": 0.5,
    "maxAlphaCo": 0.5,
    "writeControl": "adjustableRunTime",
    "maxDeltaT": dt_yield/10.0
}

# --- growth solve: not Courant-limited, takes big steps ---
GROWTH_CONTROL = {
    "deltaT": 10.0,
    "adjustTimeStep": False,
    "writeControl": "runTime",
    "maxDeltaT": 10.0
}

# --- yield stress model -----------------------------------------------------
# tau0 scales with the biofilm volume fraction: dilute biofilm yields more easily.
# Set TAU0_EXPONENT = 0 for a uniform tau0 = tauBiofilm wherever B > B_solid.
tauBiofilm = 4e-3           # yield stress of dense biofilm [N/m^2]
TAU0_EXPONENT = 1.0         # tau0 = tauBiofilm * B**TAU0_EXPONENT
B_solid = 0.05              # below this, treat as fluid (tau0 = 0)

# --- initial seeding ---
n_init = 20                 # number of blobs
r0_frac = 1 / 25            # blob radius as a fraction of W
B_init = 0.9                # seeded volume fraction
SEED = 12345                # fixed, so the seeding is reproducible

# --- parallel ---
RUN_PARALLEL = False
N_CPUS = 4
DECOMP_METHOD = "scotch"


# ============================================================
# Helpers
# ============================================================
_DECOMPOSE_DICT = """\
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}}
numberOfSubdomains {n};
method          {method};
"""


def prepare_parallel(case, case_path, n, method):
    """Wipe stale processor dirs, write decomposeParDict, decompose the mesh."""
    for name in os.listdir(case_path):
        if name.startswith("processor") and os.path.isdir(os.path.join(case_path, name)):
            shutil.rmtree(os.path.join(case_path, name))
    with open(os.path.join(case_path, "system", "decomposeParDict"), "w") as fh:
        fh.write(_DECOMPOSE_DICT.format(n=n, method=method))
    case.run(["decomposePar", "-force"])


def run_solver(case, solver, control):
    """Set the per-solver controlDict entries, then run.

    controlDict is shared, so EVERY entry that differs between the two solvers
    must be set on every call -- not just endTime. A leftover adjustTimeStep=yes
    from the yield solve would silently Courant-limit the growth solve to a
    crawl; a leftover deltaT=10 would blow up the VOF interface.
    """
    with case.control_dict as f:
        for key, value in control.items():
            f[key] = value

    if RUN_PARALLEL:
        case.run(["decomposePar", "-fields", "-time", str(case[-1].name), "-force"])
        case.run(solver, parallel=True, cpus=N_CPUS)
        case.run(["reconstructPar", "-latestTime"])
    else:
        case.run(solver)


def read_field(case, name):
    with case[-1][name] as f:
        return np.asarray(f.internal_field, dtype=float)


def write_field(case, name, values):
    with case[-1][name] as f:
        f.internal_field = values


def tau0_from_B(B):
    """Spatially varying yield stress from the biofilm volume fraction."""
    tau0 = tauBiofilm * np.power(np.clip(B, 0.0, 1.0), TAU0_EXPONENT)
    tau0[B < B_solid] = 0.0
    return tau0

def reset_inherited_deltaT(case, deltaT):
    """Overwrite the deltaT that startFrom=latestTime inherits from uniform/time.

    With adjustTimeStep, OpenFOAM writes the adapted deltaT into <time>/uniform/time,
    and startFrom latestTime reads it back -- overriding controlDict. After a growth
    step (deltaT ~ 10s), that inherited value makes the yield solver compute a first
    Courant number in the thousands and exit. Force it back to the yield start step.
    """
    tf = os.path.join(case.path, case[-1].name, "uniform", "time")
    if os.path.isfile(tf):
        case.file(os.path.join(case[-1].name, "uniform", "time"))["deltaT"] = deltaT

def B_from_alpha(case):
    """Read the biofilm phase fraction back out of the VOF field.

    Warns if alpha.water is out of bounds by more than round-off: MULES is bounded
    to solver tolerance, not exactly, and a large excursion means the interface
    scheme is misbehaving (tighten maxAlphaCo / raise nAlphaSubCycles).
    """
    alpha = read_field(case, "alpha.water")
    over = np.abs(alpha - np.clip(alpha, 0.0, 1.0)).max()
    if over > 1e-6:
        print(f"  WARNING: alpha.water out of bounds by {over:.2e} "
              f"-- check MULES / Courant")
    return np.clip(1.0 - alpha, 0.0, 1.0)


def sync_yield_fields_from_B(case, BField):
    """Push the grown biofilm back into the fields the yield solver reads."""
    write_field(case, "alpha.water", 1.0 - BField)
    write_field(case, "tau0", tau0_from_B(BField))


def reset_zero_dir(case_path):
    """Restore a pristine 0/ from 0.orig.

    Essential before blockMesh + subsetMesh: a leftover nonuniform field sized to
    the previous (subset) mesh makes subsetMesh fail on the field-length check.
    """
    zero = os.path.join(case_path, "0")
    orig = os.path.join(case_path, "0.orig")
    if not os.path.isdir(orig):
        raise FileNotFoundError(f"{orig} missing -- needed to reset {zero}")
    if os.path.isdir(zero):
        shutil.rmtree(zero)
    shutil.copytree(orig, zero)


def run_yield_adaptive(case, t_start):
    """Relax in chunks from t_start, stopping once the biofilm has settled.

    STOP CONDITION: small nu2 == still deforming, large nu2 == settled. So we keep
    stepping while nu2 is below threshold SOMEWHERE in the biofilm, and stop once
    it is above threshold EVERYWHERE. If your nu2 convention is inverted, flip the
    comparison on the line marked below.

    Returns the total relaxation time actually simulated.
    """
    t = 0.0
    k = 0
    while t < t_yield_max - 1e-9:
        k += 1
        t_next = min(dt_yield_chunk * k, t_yield_max)

        control = dict(YIELD_CONTROL)
        control["endTime"] = t_start + t_next
        control["writeInterval"] = dt_yield_chunk
        run_solver(case, YIELD_SOLVER, control)
        t = t_next

        try:
            nu2 = read_field(case, "nu2")
        except Exception as e:
            print(f"    [yield] WARNING: could not read nu2 ({e!r}); "
                  f"relaxing the full {t_yield_max:.2f}s")
            continue

        alpha = read_field(case, "alpha.water")
        region = alpha < alpha_biofilm_cut
        if not np.any(region):
            print(f"    [yield] t={t:.2f}s: no biofilm-only cells yet")
            continue

        agg = nu2[region].min()
        print(f"    [yield] t={t:.2f}s  min nu2 (biofilm) = {agg:.3e}")
        if agg >= nu2_threshold:      # <-- stop condition: flip if nu2 is inverted
            print(f"    [yield] settled after {t:.2f}s")
            break

    return t


# ============================================================
# Mask, mesh, cell <-> image mapping
# ============================================================
rng = np.random.default_rng(SEED)

mask, openfoam_image_field, Nx, Ny = utils.image_to_openfoam_mask(
    image_path=IMAGE_PATH, method="mean", invert=False)

W = Ny / Nx * L
dx = L / Nx
dy = W / Ny

print(f"domain {L*1e6:.0f} x {W*1e6:.0f} um, mesh {Nx} x {Ny}, dx = {dx*1e6:.2f} um")

fresh = START_STEP == 1
case = foamlib.FoamCase(casePath)

if fresh:
    case.clean()
    reset_zero_dir(casePath)

    utils.create_blockmesh_dict(
        Nx, Ny, L, W, thickness=thickness,
        output_dir=os.path.join(casePath, "system"))
    case.run("blockMesh")

    # carve out the pillars
    cell_id_grains = np.nonzero(openfoam_image_field)[0]
    utils.create_topoSetDict(
        cell_id_grains.tolist(),
        output_path=os.path.join(casePath, "system/topoSetDict"))
    case.run("topoSet")
    case.run("subsetMesh grainCells -overwrite -patch walls")

# Cell <-> image mapping. Deterministic given (mesh, image), so it is rebuilt every
# run rather than stored. Only used for visualisation. MUST come after subsetMesh:
# cell indices refer to the reduced mesh.
x, y, z = fluidfoam.readmesh(casePath)
x = np.asarray(x)
y = np.asarray(y)

xn = (x / dx).astype(int)
yn = (y / dy).astype(int)

image_cell_values = np.array(mask, dtype=int) - 1     # -1 marks a grain pixel
image_cell_values[yn, xn] = np.arange(len(x))

X, Y = np.meshgrid(np.linspace(0, L, Nx), np.linspace(0, W, Ny))


# ============================================================
# Seed the biofilm
# ============================================================
if fresh:
    r0 = W * r0_frac
    B2D = np.zeros_like(mask, dtype=float)

    boundary_pixels = utils.get_boundaries(mask)      # fluid pixels touching a grain
    for _ in range(n_init):
        loc = rng.choice(boundary_pixels, axis=0)
        x0 = X[loc[0], loc[1]]
        y0 = Y[loc[0], loc[1]]
        B2D[(X - x0) ** 2 + (Y - y0) ** 2 < r0 ** 2] = B_init

    B2D[~mask] = 0.0
    BField = utils.to_1d(B2D, image_cell_values)

    write_field(case, "B", BField)
    sync_yield_fields_from_B(case, BField)

    print(f"seeded {n_init} blobs, {np.count_nonzero(B2D)} cells, B = {B_init}")

if RUN_PARALLEL:
    prepare_parallel(case, casePath, N_CPUS, DECOMP_METHOD)

os.makedirs(SNAPSHOT_DIR, exist_ok=True)


# ============================================================
# Time loop
# ============================================================
t_bio = 0.0     # accumulated growth time (the biologically meaningful clock)

for step in range(START_STEP, no_steps):
    t_wall = time.perf_counter()
    t_case = float(case[-1].name)
    print(f"\n=== step {step} / {no_steps - 1}  "
          f"(case t = {t_case:.1f}s, bio t = {t_bio/3600:.1f}h) ===")

    # ---- 1. mechanical relaxation: deform / detach under flow ----
    # set the timestep back to 1e-8 first
    reset_inherited_deltaT(case, 1e-8)
    if YIELD_ADAPTIVE:
        t_relax = run_yield_adaptive(case, t_case)
    else:
        control = dict(YIELD_CONTROL)
        control["endTime"] = t_case + dt_yield
        control["writeInterval"] = dt_yield/5.0
        control["startTime"] = t_case
        print(control)
        print(f"  case[-1].name = {case[-1].name}, setting endTime = {t_case + dt_yield}")
        run_solver(case, YIELD_SOLVER, control)
        t_relax = dt_yield
    print(f"  yield relaxed for {t_relax:.2f}s")

    # input("Press Enter to continue...")

    # the biofilm has been advected, so B must be re-read from the VOF field
    BField = B_from_alpha(case)
    write_field(case, "B", BField)

    # ---- 2. grow for one dt ----
    # phi is already in this time dir, written by the yield solver's pressure
    # correction -- nothing to copy, and it is divergence-free by construction.
    t_case = float(case[-1].name)
    control = dict(GROWTH_CONTROL)
    control["endTime"] = t_case + dt
    control["startTime"] = t_case
    control["writeInterval"] = dt          # write only at the end of the step
    reset_inherited_deltaT(case, 10.0)
    run_solver(case, GROWTH_SOLVER, control)

    t_bio += dt

    # ---- 3. push the grown biofilm back into the yield solver's fields ----
    BField = read_field(case, "B")
    sync_yield_fields_from_B(case, BField)

    # ---- snapshot for visualisation ----
    fields2D = {"B": utils.to_2d(BField, image_cell_values)}
    for name in ("C", "A"):
        fields2D[name] = utils.to_2d(read_field(case, name), image_cell_values)

    # np.savez_compressed(
    #     os.path.join(SNAPSHOT_DIR, f"step_{step:05d}.npz"),
    #     mask=mask, X=X, Y=Y, t=t_bio, **fields2D)

    Bmax_now = fields2D["B"][mask].max()
    print(f"  B_max = {Bmax_now:.4f}   B_total = {BField.sum():.4e}   "
          f"C_min = {fields2D['C'][mask].min():.3e}   "
          f"A_max = {fields2D['A'][mask].max():.3e}")

    if Bmax_now > 0.99:
        print("  WARNING: B is at the cap -- check capFrac / Picard convergence")

    print(f"  step took {time.perf_counter() - t_wall:.1f}s")

print("\nDone.")