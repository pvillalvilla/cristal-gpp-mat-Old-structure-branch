function [values_bias_ku,values_bias_ka, bias_met_ku, bias_met_ka] = bias_per_zones(measured_values, ...
        expected_values, bias_requirement, ku_start_lat_idx, ku_end_lat_idx, ka_start_lat_idx, ...
        ka_end_lat_idx)
    N_zones = length(expected_values);
    for n = 1:N_zones
        mean_val_ku = mean(measured_values(ku_start_lat_idx(n):ku_end_lat_idx(n)));
        mean_val_ka = mean(measured_values(ka_start_lat_idx(n):ka_end_lat_idx(n)));
        
        values_bias_ku(n) = mean_val_ku - expected_values(n);
        bias_met_ku(n) = values_bias_ku(n) < bias_requirement;
        values_bias_ka(n) = mean_val_ka - expected_values(n);
        bias_met_ka(n) = values_bias_ka(n) < bias_requirement;
    end
end