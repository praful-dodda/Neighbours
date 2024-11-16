function [psub, zsub, dsub, nsub, index] = neighbours_stug(p0, netcdf_data, nmax, dmax)
    
    % Extract the relevant fields from the NetCDF data
    x = ncread(netcdf_data, 'x');
    y = ncread(netcdf_data, 'y');
    time = ncread(netcdf_data, 'time');
    Lon = ncread(netcdf_data, 'Lon');
    Lat = ncread(netcdf_data, 'Lat');
    BENZ = ncread(netcdf_data, 'BENZ');

    % Get the indices from p0
    idx_x = p0(4);
    idx_y = p0(5);
    idx_t = p0(6);

    % Define the search limits based on nmax
    x_min = max(1, idx_x - ceil(sqrt(nmax)));
    x_max = min(length(x), idx_x + ceil(sqrt(nmax)));
    y_min = max(1, idx_y - ceil(sqrt(nmax)));
    y_max = min(length(y), idx_y + ceil(sqrt(nmax)));
    
    t_min = max(1, idx_t - nmax);
    t_max = min(length(time), idx_t + nmax);

    % Extract the subset of coordinates and data within the nmax limits
    [X, Y, T] = ndgrid(x(x_min:x_max), y(y_min:y_max), time(t_min:t_max));
    psub_all = [Lon(x_min:x_max, y_min:y_max), Lat(x_min:x_max, y_min:y_max), T(:), ...
                X(:), Y(:), T(:)];
    
    % Calculate the spatial and temporal distances
    spatial_distances_all = sqrt((psub_all(:,1) - p0(1)).^2 + (psub_all(:,2) - p0(2)).^2);
    temporal_distances_all = abs(psub_all(:,3) - p0(3));
    space_time_distances_all = spatial_distances_all + dmax(3) * temporal_distances_all;

    % Sort by combined space-time distance and take the closest nmax points
    [~, sorted_idx_all] = sort(space_time_distances_all);
    closest_points_idx_all = sorted_idx_all(1:min(nmax, numel(sorted_idx_all)));

    % Apply dmax constraints on the closest nmax points
    valid_points_idx = closest_points_idx_all(spatial_distances_all(closest_points_idx_all) <= dmax(1) & ...
                                              temporal_distances_all(closest_points_idx_all) <= dmax(2));

    % Prepare the outputs
    psub = psub_all(valid_points_idx, :);
    linear_idx = sub2ind(size(BENZ), psub(:,4), psub(:,5), psub(:,6));
    zsub = BENZ(linear_idx);
    dsub = [spatial_distances_all(valid_points_idx), temporal_distances_all(valid_points_idx)];
    nsub = length(zsub);
    index = valid_points_idx;
end