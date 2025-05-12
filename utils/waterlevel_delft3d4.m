% ------------------------------
% File: utils/waterlevel_delft3d4.m
% Purpose: Convert Delft3D trim output to NetCDF + CSV
% Author: Shabnam
% Modified: May 2025
% ------------------------------

clc; clear;

% --- 1. Add Delft3D MATLAB toolbox to the path
addpath(fullfile(fileparts(mfilename('fullpath')), 'delft3d_matlab'));

% --- 2. Set model output name (base filename)
name_output = '45x45';

% --- 3. Open the Delft3D trim file (NEFIS format)
fprintf('Opening Delft3D trim file: ../trim-%s.dat\n', name_output);
trim = vs_use(['../trim-', name_output, '.dat'], 'quiet');

% --- 4. Extract water level data (S1) [nt x M x N]
fprintf('Extracting water level (S1)...\n');
water = vs_let(trim, 'map-series', 'S1', 'quiet');
[nt, m, n] = size(water);
points = m * n;

% --- 5. Flatten water data to [nt x points]
water_flat = reshape(water, [nt, points]);

% --- 6. Create time vector (adjust dt_sec if needed)
dt_sec = 600;  % Time interval in seconds
time = (0:nt-1)' * dt_sec;

% --- 7. Extract real-world X and Y coordinates
x = squeeze(vs_let(trim, 'map-const', 'XCOR', 'quiet'));  % [M x N]
y = squeeze(vs_let(trim, 'map-const', 'YCOR', 'quiet'));  % [M x N]
x_coords = reshape(x, [1, points]);
y_coords = reshape(y, [1, points]);

% --- 8. Write to NetCDF
outfile = '../trim_45x45.nc';
fprintf('Saving to NetCDF: %s\n', outfile);
if exist(outfile, 'file'); delete(outfile); end

nccreate(outfile, 'time',        'Dimensions', {'time', nt},      'Datatype', 'double');
nccreate(outfile, 'x',           'Dimensions', {'point', points}, 'Datatype', 'double');
nccreate(outfile, 'y',           'Dimensions', {'point', points}, 'Datatype', 'double');
nccreate(outfile, 'water_level', 'Dimensions', {'time', nt, 'point', points}, 'Datatype', 'double');

ncwrite(outfile, 'time',        time);
ncwrite(outfile, 'x',           x_coords);
ncwrite(outfile, 'y',           y_coords);
ncwrite(outfile, 'water_level', water_flat);

fprintf('✅ NetCDF file written: %s\n', outfile);

% --- 9. Save as CSV
csv_out   = '../water_level_extracted.csv';
coord_out = '../grid_coords_extracted.csv';

% --- Time series CSV
csv_data = [time, water_flat];
header = ['time', arrayfun(@(i) sprintf('pt%d', i), 1:points, 'UniformOutput', false)];
header_line = strjoin(header, ',');

fid = fopen(csv_out, 'w');
fprintf(fid, '%s\n', header_line);
fclose(fid);
dlmwrite(csv_out, csv_data, '-append');

% --- Grid coordinate CSV
coord_data = [ (1:points)', x_coords', y_coords' ];
fid = fopen(coord_out, 'w');
fprintf(fid, 'pt,x,y\n');
fclose(fid);
dlmwrite(coord_out, coord_data, '-append');

fprintf('✅ CSV files written to: %s and %s\n', csv_out, coord_out);
