import numpy as np
from PIL import Image
import os
import skimage

def image_to_openfoam_mask(image_path, Nx=None, Ny=None, method='threshold', 
                          threshold=128, invert=False):
    """
    Convert image to OpenFOAM boolean field.
    
    Returns:
        tuple: (mask_2d, openfoam_field, Nx, Ny)
    """
    # Read image
    image = Image.open(image_path).convert('L')
    
    # Resize if specified
    if Nx and Ny:
        image = image.resize((Nx, Ny))
    else:
        Nx, Ny = image.size
    
    image_array = np.array(image)
    
    # Create mask based on method
    if method == 'threshold':
        thresh_val = threshold
    elif method == 'mean':
        thresh_val = np.mean(image_array)
    
    mask_2d = image_array > thresh_val
    
    if invert:
        mask_2d = ~mask_2d
    
    # Convert to OpenFOAM field (as integers: 1=True, 0=False)
    openfoam_field = np.zeros(Nx * Ny, dtype=int)
    
    for j in range(Ny):
        for i in range(Nx):
            cell_id = i + j * Nx
            openfoam_field[cell_id] = int(mask_2d[j, i])
    
    print(f"Image size: {Nx} x {Ny}")
    print(f"Threshold: {thresh_val:.1f}")
    print(f"True cells: {np.sum(openfoam_field)} / {len(openfoam_field)}")
    
    return mask_2d, openfoam_field, Nx, Ny



def create_blockmesh_dict(Nx, Ny, L, W, thickness=0.1, output_dir="system"):
    """
    Create blockMeshDict for a 2D rectangular mesh.
    
    Args:
        Nx (int): Number of cells in x-direction (length)
        Ny (int): Number of cells in y-direction (width)  
        L (float): Length in x-direction
        W (float): Width in y-direction
        thickness (float): Thickness in z-direction (for 2D case)
        output_dir (str): Directory to write blockMeshDict
        
    Returns:
        str: Content of blockMeshDict file
    """
    
    # Calculate cell sizes
    dx = L / Nx
    dy = W / Ny
    
    # Define vertices (8 vertices for a hexahedral block)
    # Bottom face (z = 0)
    x0, y0, z0 = 0.0, 0.0, 0.0
    x1, y1, z1 = L, 0.0, 0.0
    x2, y2, z2 = L, W, 0.0
    x3, y3, z3 = 0.0, W, 0.0
    
    # Top face (z = thickness)
    x4, y4, z4 = 0.0, 0.0, thickness
    x5, y5, z5 = L, 0.0, thickness
    x6, y6, z6 = L, W, thickness
    x7, y7, z7 = 0.0, W, thickness
    
    blockmesh_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v7                                    |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Mesh parameters
// Length (x): {L} m, Width (y): {W} m, Thickness (z): {thickness} m
// Cells: {Nx} x {Ny} x 1
// Cell size: dx = {dx:.6f} m, dy = {dy:.6f} m

scale   1;

vertices
(
    ({x0} {y0} {z0})    // vertex 0: origin bottom
    ({x1} {y1} {z1})    // vertex 1: +x bottom  
    ({x2} {y2} {z2})    // vertex 2: +x+y bottom
    ({x3} {y3} {z3})    // vertex 3: +y bottom
    ({x4} {y4} {z4})    // vertex 4: origin top
    ({x5} {y5} {z5})    // vertex 5: +x top
    ({x6} {y6} {z6})    // vertex 6: +x+y top  
    ({x7} {y7} {z7})    // vertex 7: +y top
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({Nx} {Ny} 1) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    walls
    {{
        type wall;
        faces
        (
            (3 7 6 2)
            (1 5 4 0)
        );
    }}
    inlet
    {{
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }}
    outlet
    {{
        type patch;
        faces
        (
            (2 6 5 1)
        );
    }}
    emptyFaces
    {{
        type empty;
        faces
        (
            (0 3 2 1)
            (4 5 6 7)
        );
    }}
);

mergePatchPairs
(
);

// ************************************************************************* //
"""

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Write file
    output_path = os.path.join(output_dir, "blockMeshDict")
    with open(output_path, 'w') as f:
        f.write(blockmesh_content)
    
    print(f"blockMeshDict written to: {output_path}")
    print(f"Mesh: {Nx} x {Ny} cells, Domain: {L} x {W} m")
    print(f"Cell size: {dx:.6f} x {dy:.6f} m")
    
    return blockmesh_content


def create_topoSetDict(cell_ids, output_path="system/topoSetDict"):
    """
    Generates an OpenFOAM topoSetDict file to create a cellSet from a list of cell IDs.

    Args:
        cell_ids (list of int): A list of the cell IDs to include in the set.
        output_path (str): The full path for the output file.
                           Defaults to 'system/topoSetDict' for a standard case structure.
    """
    # Ensure the directory for the output file exists
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")

    # Convert the list of cell IDs into a single space-separated string
    cell_string = ' '.join(map(str, cell_ids))

    # Use an f-string to create the file content
    file_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

actions
(
    {{
        name    grainCells;
        type    cellSet;
        action  new;
        source  labelToCell;
        sourceInfo
        {{
            value ({cell_string});
        }}
    }}
);

// ************************************************************************* //
"""

    # Write the content to the specified file
    try:
        with open(output_path, 'w') as f:
            f.write(file_content)
        print(f"✅ Successfully created topoSetDict at: {output_path}")
    except IOError as e:
        print(f"❌ Error writing to file: {e}")

def to_1d(values_2d, index_map):
    flat_values = values_2d[index_map >= 0]
    flat_indices = index_map[index_map >= 0]
    
    max_index = flat_indices.max()
    array_1d = np.empty(max_index + 1, dtype=values_2d.dtype)
    array_1d[flat_indices] = flat_values
    return array_1d

def to_2d(array_1d, index_map):
    values_2d = np.full(index_map.shape, fill_value=np.nan, dtype=array_1d.dtype)
    valid = index_map >= 0
    values_2d[valid] = array_1d[index_map[valid]]
    return values_2d

def get_boundaries(binary_image):
    # 4-connectivity (only horizontal and vertical neighbors)
    boundaries = skimage.segmentation.find_boundaries(binary_image.astype(int), 
                                                   connectivity=1,  # 4-connectivity
                                                   mode='inner')
    # 
    # # 8-connectivity (includes diagonal neighbors)
    # boundaries = skimage.segmentation.find_boundaries(binary_image.astype(int), 
    #                                                connectivity=2,  # 8-connectivity  
    #                                                mode='inner')

    boundary_coords = np.column_stack(np.where(boundaries))
    
    return boundary_coords


