function [requirements_values, requirements_met]= requirements_validation_TIW1(L1B_HR,chd)
    %% SIRS-793 Water Surface Height uncertainty < 10 cm
    % Ku Band
    [requirements_values.WSH_uncertainty_ku, requirements_met.WSH_uncertainty_ku] = WSH_uncertainty(L1B_HR.altitude_ku, ...
        L1B_HR.range_ku, L1B_HR.time_ku, chd.SSH_uncertainty, chd.ku_start_lat_idx, chd.ku_end_lat_idx);
    % Ka Band
    [requirements_values.WSH_uncertainty_ka, requirements_met.WSH_uncertainty_ka] = WSH_uncertainty(L1B_HR.altitude_ka, ...
        L1B_HR.range_ka, L1B_HR.time_ka, chd.SSH_uncertainty, chd.ka_start_lat_idx, chd.ka_end_lat_idx);
    
end