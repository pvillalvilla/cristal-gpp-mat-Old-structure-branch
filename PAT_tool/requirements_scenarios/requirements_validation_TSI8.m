function [requirements_values, requirements_met] = requirements_validation_TSI8 (GR,chd,cst,cnf)
    % SARIn OB Ku / SAR OB Ka GR HR
    %% SIRS-276: Snow depth uncertainty < 0.05 m
    [requirements_values.snow_depth_uncertainty, requirements_met.snow_depth_uncertainty] = ...
        snow_depth_uncertainty(GR.range_ku, GR.range_ka, GR.time_ku, GR.time_ka, ...
        chd.snow_depth_uncertainty, chd.ku_start_lat_idx, chd.ku_end_lat_idx, ...
        chd.ka_start_lat_idx, chd.ka_end_lat_idx, cnf.mode_improved_approximate, cst);
 
    %% bias (expected snow depth vs obtained snow depth)
    [requirements_values.bias_snow_depth_ku, requirements_values.bias_snow_depth_ka, requirements_met.bias_snow_depth_ku, ...
        requirements_met.bias_snow_depth_ka] = bias_per_zones(GR.snow_depth, ...
        chd.snow_depth, chd.snow_depth_bias, chd.ku_start_lat_idx, chd.ku_end_lat_idx, chd.ka_start_lat_idx, ...
        chd.ka_end_lat_idx);
    
end