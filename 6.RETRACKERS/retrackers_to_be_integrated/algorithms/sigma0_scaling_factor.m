function [data]   = sigma0_scaling_factor (data)


global N_total_pulses_b_chd
global power_tx_ant_ku_chd antenna_gain_ku_chd wv_length_ku
global c_cst pi_cst
global cnf_p
global earth_radius_cst
global PTR_width_chd prf_chd





%   > Surface area (of every beam)
norm_vel_sat = data.GEO.V;
range_sat_surf = c_cst/2.*data.MEA.win_delay;

azimuth_distance = (1+range_sat_surf./earth_radius_cst).*wv_length_ku .* range_sat_surf ./ (1./prf_chd) ./ (2 * (norm_vel_sat) .* N_total_pulses_b_chd);

range_distance = 2 * sqrt (c_cst * range_sat_surf * (PTR_width_chd).*...
    earth_radius_cst./(earth_radius_cst + range_sat_surf) );

switch lower(cnf_p.window_type_a)
    case 'hamming'
        wf_a=1.486*0.92;%widening factor as denoted by Dinardo
    case 'hanning'
        wf_a=1.0;% TBD
    otherwise
        wf_a=1.0;
end

switch lower(cnf_p.window_type_r)
    case 'hamming'
        wf_r=1.486*0.92;%widening factor as denoted by Dinardo
    case 'hanning'
        wf_r=1.0;% TBD
    otherwise
        wf_r=1.0;
end

surface_area = (wf_a.*azimuth_distance .* (wf_r.*range_distance)).*0.886;

%added by EM 30.11.2015: consider the surface and not the mean over
%the different beams & compensate for norm. in fft range
%compression & TBP or pulse compression gain
data.HRM.s0_sf = 10*log10(64) + 30*log10(pi_cst)...
    + 40*log10(range_sat_surf) - 10*log10(power_tx_ant_ku_chd) - 2*(antenna_gain_ku_chd)...
    - 20*log10(wv_length_ku) - 10*log10(surface_area);




    
    
    
    
    
    
    
end