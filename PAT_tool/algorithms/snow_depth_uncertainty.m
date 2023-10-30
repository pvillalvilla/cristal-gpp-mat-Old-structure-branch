function [snow_depth_uncertainty_values, snow_depth_uncertainty_met] = snow_depth_uncertainty(range_ku, ...
    range_ka, time_ku, time_ka, snow_depth_uncertainty_req, start_lat_idx_ku, end_lat_idx_ku, ...
    start_lat_idx_ka, end_lat_idx_ka, mode_improved_approximate, cst)
    
    N_zones = length(start_lat_idx);
    
    for n = 1:N_zones
        range_ku_corrected = correct_slope(range_ku, time_ku, start_lat_idx_ku(n), end_lat_idx_ku(n)); % correct range_ku slope
        range_ka_corrected = correct_slope(range_ka, time_ka, start_lat_idx_ka(n), end_lat_idx_ka(n)); % correct range_ka slope
        stdev_range_ku = std(range_ku_corrected); % stdev range Ku
        stdev_range_ka = std(range_ka_corrected); % stdev range Ka
        
        if mode_improved_approximate==0
            range_ka_corrected = decimate(range_ka_corrected, 2.6481); 
        end
            
        correlation_ku_ka = corrcoef(range_ku_corrected, range_ka_corrected); % correlation between range Ku a & range Ka
        C=cst.c_snow/cst.c;
    
        snow_depth_uncertainty_values(n) = ...
            C*sqrt(stdev_range_ku^2 + stdev_range_ka^2 - 2*correlation_ku_ka); % propagation of uncertainty
        snow_depth_uncertainty_met(n) = snow_depth_uncertainty_values(n) ...
            < snow_depth_uncertainty_req; % 0.05;
    end
    
end