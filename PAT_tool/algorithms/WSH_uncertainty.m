function [WSH_uncertainty, WSH_uncertainty_met] = WSH_uncertainty(altitude, ...
    range, time, SSH_uncertainty_required, start_lat_idx, end_lat_idx)

N_zones = length(start_lat_idx);
for n = 1:N_zones
     
    range_corrected = correct_slope(range, time, start_lat_idx(n), end_lat_idx(n)); % correct slope range
    altitude_corrected = correct_slope(altitude, time, start_lat_idx(n), end_lat_idx(n)); %correct slope altitude 

    stdev_range = std(range_corrected); % stdev range 
    stdev_altitude = std(altitude_corrected); % stdev altitude 
    
    correlation_range_altitude = corrcoef(range_corrected, altitude_corrected); % correlation between range & altitude
    
    WSH_uncertainty(n) = C*sqrt(stdev_range^2 + stdev_altitude^2 - 2*correlation_range_altitude); % stdev WSH
    
    WSH_uncertainty_met(n) = WSH_uncertainty(n) < SSH_uncertainty_required; % 0.1

end

end