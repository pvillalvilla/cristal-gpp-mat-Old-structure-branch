function [range_random_error_value, range_random_error_met] = ...
    range_random_error(range, point_target_idx, required_range_random_error, N_central_beams)
idx_range = point_target_idx - floor(N_central_beams/2):point_target_idx + floor(N_central_beams/2);
range_random_error_value = stdev(range(idx_range));
range_random_error_met = range_random_error_value < required_range_random_error;
end