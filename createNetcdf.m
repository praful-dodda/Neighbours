filename = '../../analysis/1data/CAMx/modelLOC2LATLON_4K.csv';
data = readtable(filename);

% Extract the variables
x = data.X;
y = data.Y;
Lon = data.LON;
Lat = data.LAT;

% Define time range
time = 1:100;

% Determine the grid dimensions
unique_x = unique(x);
unique_y = unique(y);
num_x = length(unique_x);
num_y = length(unique_y);
num_time = length(time);

% Ensure that the data represents a regular grid
if num_x * num_y ~= length(x)
    error('The provided data does not form a regular grid.');
end

% Reshape the spatial coordinates into 2D grids
[X, Y] = meshgrid(unique_x, unique_y);

% Generate random data for the 'BENZ' field
% Random data with dimensions (num_x x num_y x num_time)
BENZ = rand(num_x, num_y, num_time);

% Create a NetCDF file
ncfile = 'SampleBENZ.nc';
ncid = netcdf.create(ncfile, 'CLOBBER');

% Define dimensions
dimid_x = netcdf.defDim(ncid, 'x', num_x);
dimid_y = netcdf.defDim(ncid, 'y', num_y);
dimid_time = netcdf.defDim(ncid, 'time', num_time);

% Define variables
varid_x = netcdf.defVar(ncid, 'x', 'double', dimid_x);
varid_y = netcdf.defVar(ncid, 'y', 'double', dimid_y);
varid_Lon = netcdf.defVar(ncid, 'Lon', 'double', [dimid_x, dimid_y]);
varid_Lat = netcdf.defVar(ncid, 'Lat', 'double', [dimid_x, dimid_y]);
varid_time = netcdf.defVar(ncid, 'time', 'double', dimid_time);
varid_BENZ = netcdf.defVar(ncid, 'BENZ', 'double', [dimid_x, dimid_y, dimid_time]);

% End definitions and leave define mode
netcdf.endDef(ncid);

% Write data to the NetCDF file
netcdf.putVar(ncid, varid_x, unique_x);
netcdf.putVar(ncid, varid_y, unique_y);
netcdf.putVar(ncid, varid_Lon, reshape(Lon, [num_x, num_y]));
netcdf.putVar(ncid, varid_Lat, reshape(Lat, [num_x, num_y]));
netcdf.putVar(ncid, varid_time, time);
netcdf.putVar(ncid, varid_BENZ, BENZ);

% Close the NetCDF file
netcdf.close(ncid);

disp('NetCDF file created successfully.');