function [L1AP] = create_L1AP_struct(cnf, chd, N_bursts)


L1AP         = [];

N_pulses_total = chd.N_pulses_burst * N_bursts;

L1AP = setfield_new(L1AP,'pri_sar',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'ambiguity_order_sar_isp',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'burst',zeros(1,N_pulses_total));
% L0 = setfield_new(L0,'burst_fbr',0);
L1AP = setfield_new(L1AP,'lat_sar_sat',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'lon_sar_sat',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'alt_sar_sat',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'alt_rate_sar_sat',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'x_vel_sat_sar',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'y_vel_sat_sar',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'z_vel_sat_sar',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'roll_sar',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'pitch_sar',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'yaw_sar',zeros(1,N_pulses_total));


L1AP = setfield_new(L1AP,'win_delay',zeros(1,N_pulses_total));
if(strcmp(cnf.mode,'SIN'))
    L1AP = setfield_new(L1AP,'win_delay_2',zeros(1,N_pulses_total));
end
%L1A = setfield_new(L1A,'h0_sar_isp',0);
%L1A = setfield_new(L1A,'cai_sar_isp',0);
%L1A = setfield_new(L1A,'cor2_sar_isp',0);
%L1A = setfield_new(L1A,'fai_sar_isp',0);
%L1A = setfield_new(L1A,'ATT1_science',0);
%L1A = setfield_new(L1A,'ATT2_science',0);
%L1A = setfield_new(L1A,'att_isp',0);
%L1A = setfield_new(L1A,'tot_fixed_gain_1',0);
%L1A = setfield_new(L1A,'tot_fixed_gain_2',0);
%L1A = setfield_new(L1A,'transmit_power',0);
L1AP = setfield_new(L1AP,'doppler_range_correction',zeros(1,N_pulses_total));
% L1A = setfield_new(L1A,'instrument_range_correction_tx_rx',0);
%L1A = setfield_new(L1A,'instrument_range_correction_rx',0);
%L1A = setfield_new(L1A,'instrument_sigma0_correction_tx_rx',0);
%L1A = setfield_new(L1A,'instrument_sigma0_correction_rx',0);
L1AP = setfield_new(L1AP,'internal_phase_correction',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'external_phase_correction',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'noise_power',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'phase_slope_correction',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'gain_corr_instr_sar_1',zeros(1,N_pulses_total));
if(strcmp(cnf.mode,'SIN'))
    L1AP = setfield_new(L1AP,'gain_corr_instr_sar_2',zeros(1,N_pulses_total));
end


L1AP = setfield_new(L1AP,'wfm_cal_gain_corrected',zeros(N_pulses_total,chd.N_samples_sar));
if(strcmp(cnf.mode,'SIN'))
    L1AP = setfield_new(L1AP,'wfm_cal_gain_corrected_2',zeros(N_pulses_total,chd.N_samples_sar));
end

L1AP = setfield_new(L1AP,'time',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'x_sar_sat',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'y_sar_sat',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'z_sar_sat',zeros(1,N_pulses_total));
L1AP = setfield_new(L1AP,'lat_sar_surf',0);
L1AP = setfield_new(L1AP,'lon_sar_surf',0);
L1AP = setfield_new(L1AP,'alt_sar_surf',0);
L1AP = setfield_new(L1AP,'x_sar_surf',0);
L1AP = setfield_new(L1AP,'y_sar_surf',0);
L1AP = setfield_new(L1AP,'z_sar_surf',0);
L1AP = setfield_new(L1AP,'doppler_ang_sar_sat',0);

L1AP = setfield_new(L1AP,'T0_sar',0);
L1AP = setfield_new(L1AP,'bri_sar',0);

%beam angles
L1AP = setfield_new(L1AP,'N_beams',0);
L1AP = setfield_new(L1AP,'surf_loc_index',0);
L1AP = setfield_new(L1AP,'beam_ang',0);
L1AP = setfield_new(L1AP,'beam_ang_nadir_index',0);
L1AP = setfield_new(L1AP,'beam_ang_index',0);
L1AP = setfield_new(L1AP,'start_beam',0);
L1AP = setfield_new(L1AP,'end_beam',0);
L1AP = setfield_new(L1AP,'beam_index',0);

%focussing
L1AP = setfield_new(L1AP,'beams_focused_shifted',zeros(N_pulses_total,chd.N_samples_sar));
if(strcmp(cnf.processing_mode,'SIN'))
    L1AP = setfield_new(L1AP,'beams_focused_shifted_2',zeros(N_pulses_total,chd.N_samples_sar));
end



end

