# File: src_delft3d4/delft3d4_dep.py
# Author: Shabnam
# Modified: May 2025
# Create a linear sloped .rgh file using .grd dimensions (for non-square grids)

import numpy as np
import os

def read_grid_shape_from_grd(grd_file):
    """
    Reads NX and NY (number of cells) from line 7 of a .grd file,
    and returns the number of nodes: NX+1, NY+1
    """
    with open(grd_file, 'r') as f:
        lines = f.readlines()
        line7 = lines[6].strip()
        parts = line7.split()
        nx = int(parts[0]) + 1  # convert cells → nodes
        ny = int(parts[1]) + 1
        return nx, ny

def create_sloped_dep_file(grd_file='45x45.grd',
                           slope_start=0.5,
                           slope_end=-0.5,
                           boundary_value=-999.0,
                           values_per_line=12):
    # Derive output DEP filename
    dep_filename = os.path.splitext(grd_file)[0] + '.dep'

    # Read node size from grid
    nx, ny = read_grid_shape_from_grd(grd_file)  # number of nodes

    # Initialize grid
    grid = np.zeros((ny, nx))  # ny rows, nx columns

    # Generate slope across columns (excluding last dummy)
    slope_values = np.linspace(slope_start, slope_end, nx - 1, endpoint=True)

    # Fill in elevation values
    for i in range(ny - 1):
        for j in range(nx - 1):
            grid[i, j] = slope_values[j]

    # Set boundary dummy values
    grid[:, -1] = boundary_value
    grid[-1, :] = boundary_value

    # Write file with Quickin-style formatting
    with open(dep_filename, 'w', encoding='utf-8') as f:
        for i in range(ny):
            row = grid[i]
            for j, val in enumerate(row):
                val_str = f'{val:.7E}'
                spacing = '  ' if val < 0 else '   '
                f.write(spacing + val_str)

                # Wrap line every N values
                if (j + 1) % values_per_line == 0 and j != nx - 1:
                    f.write('\n')
            f.write('\n')

    print(f"✅ DEP file '{dep_filename}' created using grid from '{grd_file}' ({nx}x{ny} nodes)")

# Run
create_sloped_dep_file('45x45.grd')
