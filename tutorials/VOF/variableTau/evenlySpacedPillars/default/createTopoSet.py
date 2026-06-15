#!/usr/bin/env python3
"""
Script to generate OpenFOAM topoSetDict file for multiple cylinders in z-direction
"""

import os

def generate_toposet_dict(cylinders, z_min=0, z_max=1, output_file="system/topoSetDict", scale=1):
    """
    Generate topoSetDict file for multiple cylinders aligned in z-direction
    
    Parameters:
    -----------
    cylinders : list of tuples
        Each tuple contains (x_center, y_center, radius)
    z_min : float
        Bottom z-coordinate of cylinders (default: 0)
    z_max : float  
        Top z-coordinate of cylinders (default: 1)
    output_file : str
        Path to output topoSetDict file
    """
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Start writing the topoSetDict file
    with open(output_file, 'w') as f:
        # Write header
        f.write("""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2012                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      topoSetDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

actions
(
""")
        
        # Write first cylinder action (creates new cellSet)
        if cylinders:
            x, y, radius = cylinders[0]
            x *= scale
            y *= scale
            radius *= scale
            f.write(f"""    {{
        name    cylinderCells;
        type    cellSet;
        action  new;
        source  cylinderToCell;
        sourceInfo
        {{
            p1      ({x:.6f} {y:.6f} {z_min:.6f});
            p2      ({x:.6f} {y:.6f} {z_max:.6f});
            radius  {radius:.6f};
        }}
    }}
""")
        
        # Write remaining cylinders (add to existing cellSet)
        for i, (x, y, radius) in enumerate(cylinders[1:], 1):
            x *= scale
            y *= scale
            radius *= scale
            f.write(f"""    {{
        name    cylinderCells;
        type    cellSet;
        action  add;
        source  cylinderToCell;
        sourceInfo
        {{
            p1      ({x:.6f} {y:.6f} {z_min:.6f});
            p2      ({x:.6f} {y:.6f} {z_max:.6f});
            radius  {radius:.6f};
        }}
    }}
""")

        # write invert
        f.write("""    {
        name    cylinderCells;
        type    cellSet;
        action  invert;
    }
""")
        
        # Write footer
        f.write(""");

// ************************************************************************* //
""")
    
    print(f"Generated topoSetDict with {len(cylinders)} cylinders")
    print(f"Output file: {output_file}")
    print("\nNext steps:")
    print("1. Run: topoSet")
    print("2. Run: subsetMesh cylinderCells -overwrite")

def main():
    # Example usage - modify this section with your cylinder data
    
    # # Method 1: Define cylinders directly
    # cylinders = [
    #     (0.2, 0.3, 0.05),  # (x, y, radius)
    #     (0.7, 0.4, 0.08),
    #     (0.5, 0.8, 0.06),
    #     (0.9, 0.7, 0.04)
    # ]


    cylinders = [
        (-0.0025, -0.0025, 0.005), 
        (-0.0025, 0.0175, 0.005), 
        (-0.0025, 0.0375, 0.005), 
        (0.0148, 0.0075, 0.005), 
        (0.0148, 0.0275, 0.005), 
        (0.0148, 0.0475, 0.005), 
        (0.0321, -0.0025, 0.005), 
        (0.0321, 0.0175, 0.005), 
        (0.0321, 0.0375, 0.005), 
        (0.0495, 0.0075, 0.005), 
        (0.0495, 0.0275, 0.005), 
        (0.0495, 0.0475, 0.005), 
        (0.0668, -0.0025, 0.005), 
        (0.0668, 0.0175, 0.005), 
        (0.0668, 0.0375, 0.005),
    ]
    
    # Method 2: Load from file (uncomment if you have a data file)
    # cylinders = load_cylinders_from_file("cylinder_data.txt")
    
    # Method 3: Generate regular pattern (uncomment for testing)
    # cylinders = generate_regular_pattern(nx=3, ny=3, spacing=0.3, radius=0.05)
    
    # Set z-range for cylinders (modify as needed)
    z_min = 0.0
    z_max = 1.0
    scale = 1e-3
    
    # Generate the topoSetDict file
    generate_toposet_dict(cylinders, z_min, z_max, scale=scale)

def load_cylinders_from_file(filename):
    """
    Load cylinder data from a text file
    Expected format: x y radius (one per line)
    """
    cylinders = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 3:
                    x, y, radius = float(parts[0]), float(parts[1]), float(parts[2])
                    cylinders.append((x, y, radius))
    return cylinders

def generate_regular_pattern(nx, ny, spacing, radius, x_offset=0, y_offset=0):
    """
    Generate a regular grid pattern of cylinders
    """
    cylinders = []
    for i in range(nx):
        for j in range(ny):
            x = x_offset + i * spacing
            y = y_offset + j * spacing
            cylinders.append((x, y, radius))
    return cylinders

if __name__ == "__main__":
    main()