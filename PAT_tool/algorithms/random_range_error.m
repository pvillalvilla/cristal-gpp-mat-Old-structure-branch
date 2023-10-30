function [range_rand_err_values, range_rand_err_met] = random_range_error(range,...
    altitude, height_ssh_expected, random_range_error_required, start_lat_idx, ...
    end_lat_idx)

N_zones = length(start_lat_idx);
error_total = [];
for n = 1:N_zones
    height_ssh_meas = altitude(start_lat_idx(n):end_lat_idx(n)) - range(start_lat_idx(n):end_lat_idx(n));
    error_zone = height_ssh_meas - height_ssh_expected;
    error_total = [error_total error_zone];
end
range_rand_err_values = stdev(error_total);
range_rand_err_met = range_rand_err_values < random_range_error_required;