function [range_slope_value, range_slope_met] = ...
    range_slope(range, time,  point_target_idx, required_range_slope, N_central_beams)

idx_range = point_target_idx - floor(N_central_beams/2):point_target_idx + floor(N_central_beams/2);
range_fitting = fit(time(idx_range), range(idx_range), 'poly1');
range_slope_value = range_fitting.p1;
range_slope_met = range_slope_value < required_range_slope;
    
end
