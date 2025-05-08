# File: utils/delft3d4_dep.py
# Author: Shabnam
# Modified: April 2025

import numpy as np

def create_sloped_dep_file(filename='45x45.dep',
                           grid_size=47,
                           slope_start=0.5,
                           slope_end=-0.5,
                           boundary_value=-999.0,
                           values_per_line=12):
    # Initialize grid
    grid = np.zeros((grid_size, grid_size))

    # Generate 46 slope values (from col 0 to 45)
    slope_values = np.linspace(slope_start, slope_end, grid_size - 1, endpoint=True)

    for i in range(grid_size - 1):       # 0 to 45 (46 rows)
        for j in range(grid_size - 1):   # 0 to 45 (46 cols)
            grid[i, j] = slope_values[j]

    # Last column and last row
    grid[:, -1] = boundary_value
    grid[-1, :] = boundary_value

    # Write file with Quickin-style formatting
    with open(filename, 'w', encoding='utf-8') as f:
        for i in range(grid_size):
            row = grid[i]
            for j, val in enumerate(row):
                val_str = f'{val:.7E}'

                # Determine spacing before value
                if j == 0:
                    spacing = '  ' if val < 0 else '   '
                else:
                    spacing = '  ' if val < 0 else '   '

                f.write(spacing + val_str)

                # Wrap line every N values
                if (j + 1) % values_per_line == 0 and j != grid_size - 1:
                    f.write('\n')

            f.write('\n')

    print(f"✅ DEP file '{filename}' created correctly in Quickin style.")

# Run
create_sloped_dep_file()
