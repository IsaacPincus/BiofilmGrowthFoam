#%%

import numpy as np
import foamlib
import pythonScripts.utils as utils
import os
import fluidfoam
from PIL import Image
import shutil

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
mask, openfoam_image_field, Nx, Ny = utils.image_to_openfoam_mask(image_path='channels.png',
                                                                   method='mean', invert=False)

# true domain size
# Nx = 360
# Ny = 90
L = 2e-3
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
Diffusivity = 1e-8

# # timestep for each removal/growth step, in seconds
dt = 1200
no_steps = 50

dt_yield = 5e-1
dt_convection = 3

# # new bits of biofilm each timestep
# n_new = 5

Bcut = 8000

# U_threshold = 1e-5

tau0default = 4e-3


# %%
# Setup yield stress case
caseYieldStress = foamlib.FoamCase(pathYieldStress) # Loads the OpenFOAM case
caseYieldStress.clean()

# creates blockMeshDict
output = utils.create_blockmesh_dict(Nx, Ny, L, W, 
                                     thickness=thickness, 
                                     output_dir=os.path.join(pathYieldStress, "system"))
caseYieldStress.run("blockMesh") # do blockMesh

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
    f["writeInterval"] = dt_yield/2

# now also setup the other case, including the same grains as in the yield stress case
caseConvection = foamlib.FoamCase(pathConvection) # Loads the OpenFOAM case
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
# # Optional: output as image to check
# if img_array.max() > 255 or img_array.min() < 0:
#     img_array = 255 * (img_array - img_array.min()) / (img_array.max() - img_array.min())
#     img_array = img_array.astype(np.uint8)
# else:
#     img_array = img_array.astype(np.uint8)
# # Convert to PIL Image and save
# img = Image.fromarray(img_array, mode='L')  # 'L' = 8-bit grayscale
# img.save("image_cell_values.png")

# now, we want to write some other fields. We want each pixel to have a value of Balive and Bdead,
# for the amount of biofilm there.
# note that meshgrid in numpy is opposite of matlab, it's (ny,nx) array indexing, like an image
X, Y = np.meshgrid(np.linspace(0,L,Nx), np.linspace(0,W,Ny))
# Parameters
r0 = W/20
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
for ii in range(20):
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


#%%
######################################################################
# time loop starts here
######################################################################
for step in range(1,no_steps):
    # run the yield stress case
    caseYieldStress.run("biofilmYieldFoam")

    # Copy B, U values to convection case
    with caseYieldStress[-1]["B"] as field:
        BField = field.internal_field
    with caseYieldStress[-1]["U"] as field:
        UField = field.internal_field

    with caseConvection[-1]["B"] as field:
        field.internal_field = BField
    with caseConvection[-1]["U"] as field:
        field.internal_field = UField

    caseConvection.run("convectionWithBiomassFoam")

    # now we grow the biofilm. For now we won't have any dead biofilm
    with caseConvection[-1]["C"] as field: # read concentration field
        CField = field.internal_field

    # for each element of the B field, update based on growth
    solidIndicator2D = ( Balive2D > 10 ) - 0.5
    solidIndicator2DMasked = np.ma.masked_array(solidIndicator2D, ~mask)
    # t = skfmm.travel_time(solidIndicator2DMasked, np.ones_like(solidIndicator2D))

    prefactor = (mu*(CField/(Ks + CField)) - kd)
    BdeadField = BDeadField + kd * BField / prefactor * (np.exp(dt*prefactor) - 1)
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
    Bdead2D[Bdead2D > 0.5] = 0

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

    shutil.copytree('./yieldStress/'+folderTimeName, 
                    './yieldStress/'+str(dt*step))

    # yield stress is given by ratio of alive and dead
    # add in a tiny bit of dead everywhere just so this doesn't lead to divide by zero errors
    Bdead2D[Bdead2D < 1e-5] = 1e-5
    tau02D = ((Balive2D * tauAlive + Bdead2D * tauDead)/(Balive2D + Bdead2D))
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

    # now let's see if it worked!