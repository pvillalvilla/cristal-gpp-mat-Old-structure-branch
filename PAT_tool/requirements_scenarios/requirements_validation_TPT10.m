function [requirements_values, requirements_met]= requirements_validation_TPT10(L1B_HR, ~, L1B_SL,...
    L1A_ku, chd, cst, cnf) % ~ corresponds to L1BS
%% Create pwr wfms SL (FF)
% ku band

if (cnf.flag_waveform == 0) % Rx1
    wfm_AC_ku = abs(L1B_SL.i_samples_rx1_ku.*L1B_SL.iq_scale_factor_rx1_ku' ...
        + 1i*L1B_SL.q_samples_rx1_ku.*L1B_SL.iq_scale_factor_rx1_ku')^2;
    L1B_HR.power_waveform_ku = L1B_HR.power_waveform_rx1_ku;
    
elseif (cnf_L2.flag_waveform == 1) % Rx2
    wfm_AC_ku = abs(L1B_SL.i_samples_rx2_ku.*L1B_SL.iq_scale_factor_rx2_ku' ...
        + 1i*L1B_SL.q_samples_rx2_ku.*L1B_SL.iq_scale_factor_rx2_ku')^2;
    L1B_HR.power_waveform_ku = L1B_HR.power_waveform_rx2_ku;
    
elseif (cnf_L2.flag_waveform == 2) % Combined
    wfm_AC_ku_rx1 = abs(L1B_SL.i_samples_rx1_ku.*L1B_SL.iq_scale_factor_rx1_ku' ...
        + 1i*L1B_SL.q_samples_rx1_ku.*L1B_SL.iq_scale_factor_rx1_ku')^2;
   wfm_AC_ku_rx2 = abs(L1B_SL.i_samples_rx2_ku.*L1B_SL.iq_scale_factor_rx2_ku' ...
        + 1i*L1B_SL.q_samples_rx2_ku.*L1B_SL.iq_scale_factor_rx2_ku')^2;
   wfm_AC_ku =  (wfm_AC_ku_rx1 + wfm_AC_ku_rx2)/2;
   
end

% 2D Re-sampling of the complex matrix
wfm_AC_interp_ku=interpolsinc_2D(wfm_AC_ku,cnf.FFt.zp_ac_TRP,cnf.FFt.zp_al_TRP);

%% DATATION & WINDOW DELAY & PRI & T0 PER PULSE (FF)

% ku band
[time_pulse_ku,~,~,~,~,L1A_ku] = ...
    datation_win_delay_pulse(L1A_ku, cnf, chd, cst);

%% OSV & ATTITUDE SELECTION PER PULSE (MEASURED SIGNAL) (FF)
[~, ~, ~,x_vel_sat_sar_pulse_ku,  y_vel_sat_sar_pulse_ku,  z_vel_sat_sar_pulse_ku,...
    ~, ~, alt_sar_sat_pulse_ku,~,~, ~, ~] = ...
    osv_attitude_pulse(L1A_ku, time_pulse_ku, cnf, chd, cst);

%% SIRS-328 Range bias < 1cm (L1B_HR)
% Ku band
[requirements_values.range_bias_HR_ku, requirements_met.range_bias_HR_ku] = ...
    range_bias(L1B_HR.tracker_range_calibrated_ku, L1B_HR.slant_range_correction_applied_ku,...
    L1B_HR.power_waveform_ku, L1B_HR.Gap_flag_ku, chd.range_PT_ku,chd.range_bias_HR, chd.left_beams_ku,...
    chd.right_beams_ku, chd.point_target_idx_ku, chd, cst);

%% SIRS- 330 Datation bias < 100  (L1B_HR)
% Ku band
[requirements_values.datation_bias_HR_ku, requirements_met.datation_bias_HR_ku] = ...
    datation_bias_HR(L1B_HR.tracker_range_calibrated_ku, L1B_HR.slant_range_correction_applied_ku,...
    chd.point_target_idx_ku, chd, cst);

%% SIRS- 330 Datation bias < 100  (L1B_SL)
% ku band
[requirements_values.datation_bias_SL_ku, requirements_met.datation_bias_SL_ku, pos_max_along_focused_ku] = ...
    datation_bias_SL(wfm_AC_interp_ku, L1B_SL.lat_ku, L1B_SL.lon_ku, x_vel_sat_sar_pulse_ku, ...
    y_vel_sat_sar_pulse_ku, z_vel_sat_sar_pulse_ku, alt_sar_sat_pulse_ku, ...
    cnf, chd, cst);

%% SIRS-332 Range bias < 1 mm (L1BS)
% ku band  
% [requirements_values.range_bias_L1BS_ku, requirements_met.range_bias_L1BS_ku] = ...
%     range_bias(L1BS.tracker_range_calibrated_ku, L1BS.slant_range_correction_applied_ku, ...
%     L1BS.power_waveform_ku, L1BS.Gap_flag_ku, chd.range_PT_ku, chd.range_bias_L1BS, chd.left_beams_ku,...
%     chd.right_beams_ku, chd.point_target_idx_ku, chd, cst);

%% SIRS-334 Range random error M central beams < 1 mm (L1BS)
% ku band
% [requirements_values.range_random_error_ku, requirements_met.range_random_error_ku] = ...
%     range_random_error(L1BS.range_ku, chd.point_target_idx_ku, chd.range_random_error, chd.N_central_beams);

%% SIRS-336 Range slope M central beams < 0.01 mm/beam (L1BS)
% ku band
% [requirements_values.range_slope_L1BS_ku, requirements_met.range_slope_L1BS_ku] = ...
%     range_slope(L1BS.range_ku, L1BS.time_ku, chd.point_target_idx_ku, chd.range_slope, chd.N_central_beams);


%% SIRS-338 Burst Air PSLR < -13 dB (L1B_SL)
% ku band
[requirements_values.burst_air_pslr_SL_ku, requirements_met.burst_air_pslr_SL_ku] = ...
    burst_air_PSLR(wfm_AC_interp_ku, pos_max_along_focused_ku, chd);

%% SIRS-339 Peak power side lobes ratio < 0.01 dB (L1B_SL)
% S'ha de mirar el ratio entre els dos side lobes? o side lobes i main lobe
end