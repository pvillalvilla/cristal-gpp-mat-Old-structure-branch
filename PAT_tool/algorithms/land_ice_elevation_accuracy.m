function [land_ice_elevation_values, land_ice_elevation_met] = ...
    land_ice_elevation_accuracy(glacier_start_lat_idx, glacier_end_lat_idx, land_ice_elevation_expected, ...
    elevation_accuracy_required, land_ice_elevation)

N_glaciers = length(glacier_start_lat_idx);
for n = 1:N_glaciers
    land_ice_elevation_values(n) = mean(land_ice_elevation_expected(n) - ...
        land_ice_elevation(glacier_start_lat_idx(n):glacier_end_lat_idx(n)));
end
land_ice_elevation_met = land_ice_elevation_values < elevation_accuracy_required;