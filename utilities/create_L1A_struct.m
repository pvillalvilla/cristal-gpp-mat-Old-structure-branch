function [L1A] = create_L1A_struct(cnf, chd)


L1A         = [];
% L0 = setfield_new(L0,'source_seq_count_sar_isp',0);
% L0 = setfield_new(L0,'days',0);
% L0 = setfield_new(L0,'seconds',0);
% L0 = setfield_new(L0,'microseconds',0);
% L0 = setfield_new(L0,'USO_correction',0);
L1A = setfield_new(L1A,'ProcessID','');
% L0 = setfield_new(L0,'source_seq_count_fbr',0);
L1A = setfield_new(L1A,'ins_id',0);
L1A = setfield_new(L1A,'ins_tracking_mode',0);
% L0 = setfield_new(L0,'ins_loop_stat',0);
L1A = setfield_new(L1A,'inst_id_sar_isp',0);
L1A = setfield_new(L1A,'pri_sar',0);
L1A = setfield_new(L1A,'ambiguity_order_sar_isp',0);
L1A = setfield_new(L1A,'burst',0);
%flag and mask for OB processing indicating if CAL pulses inside beam 
L1A = setfield_new(L1A,'OB_cal_flag',4); %1=1st burst,4 holes; 2=2nd burst,3 holes; 3=3rd burst,1 hole; 4=4-12 burst, 0 holes
L1A = setfield_new(L1A,'OB_cal_mask',ones(1,64));
% L0 = setfield_new(L0,'burst_fbr',0);
L1A = setfield_new(L1A,'lat_sar_sat',0);
L1A = setfield_new(L1A,'lon_sar_sat',0);
L1A = setfield_new(L1A,'alt_sar_sat',0);
L1A = setfield_new(L1A,'alt_rate_sar_sat',0);
L1A = setfield_new(L1A,'x_vel_sat_sar',0);
L1A = setfield_new(L1A,'y_vel_sat_sar',0);
L1A = setfield_new(L1A,'z_vel_sat_sar',0);
L1A = setfield_new(L1A,'roll_sar',0);
L1A = setfield_new(L1A,'pitch_sar',0);
L1A = setfield_new(L1A,'yaw_sar',0);
L1A = setfield_new(L1A,'confi_block_degraded',[]);
L1A = setfield_new(L1A,'mea_conf_data',0);
L1A = setfield_new(L1A,'win_delay',0);
if(strcmp(cnf.mode,'SIN'))
    L1A = setfield_new(L1A,'win_delay_2',0);
end
L1A = setfield_new(L1A,'h0_comp_sar_isp',0);
L1A = setfield_new(L1A,'cor2_comp_sar_isp',0);
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
L1A = setfield_new(L1A,'doppler_range_correction',0);
% L1A = setfield_new(L1A,'instrument_range_correction_tx_rx',0);
%L1A = setfield_new(L1A,'instrument_range_correction_rx',0);
%L1A = setfield_new(L1A,'instrument_sigma0_correction_tx_rx',0);
%L1A = setfield_new(L1A,'instrument_sigma0_correction_rx',0);
L1A = setfield_new(L1A,'internal_phase_correction',0);
L1A = setfield_new(L1A,'external_phase_correction',0);
L1A = setfield_new(L1A,'noise_power',0);
L1A = setfield_new(L1A,'phase_slope_correction',0);
L1A = setfield_new(L1A,'gain_corr_instr_sar_1',0);
if(strcmp(cnf.mode,'SIN'))
    L1A = setfield_new(L1A,'gain_corr_instr_sar_2',0);
end

L1A = setfield_new(L1A,'dry_tropo_correction_bursts',0);
L1A = setfield_new(L1A,'wet_tropo_correction_bursts',0);
L1A = setfield_new(L1A,'inverse_baro_correction_bursts',0);
L1A = setfield_new(L1A,'Dynamic_atmospheric_correction_bursts',0);
L1A = setfield_new(L1A,'GIM_iono_correction_bursts',0);
L1A = setfield_new(L1A,'model_iono_correction_bursts',0);
L1A = setfield_new(L1A,'ocean_equilibrium_tide_bursts',0);
L1A = setfield_new(L1A,'long_period_tide_height_bursts',0);
L1A = setfield_new(L1A,'ocean_loading_tide_bursts',0);
L1A = setfield_new(L1A,'solid_earth_tide_bursts',0);
L1A = setfield_new(L1A,'geocentric_polar_tide_bursts',0);
L1A = setfield_new(L1A,'surface_type_flag_bursts',0);
L1A = setfield_new(L1A,'nimp_sar_isp',0);
L1A = setfield_new(L1A,'wfm_cal_gain_corrected',zeros(chd.N_pulses_burst,chd.N_samples_sar));
if(strcmp(cnf.mode,'SIN'))
    L1A = setfield_new(L1A,'wfm_cal_gain_corrected_2',zeros(chd.N_pulses_burst,chd.N_samples_sar));
end
L1A = setfield_new(L1A,'time_rx_1st',0);
L1A = setfield_new(L1A,'time',0);
L1A = setfield_new(L1A,'x_sar_sat',0);
L1A = setfield_new(L1A,'y_sar_sat',0);
L1A = setfield_new(L1A,'z_sar_sat',0);
L1A = setfield_new(L1A,'lat_sar_surf',0);
L1A = setfield_new(L1A,'lon_sar_surf',0);
L1A = setfield_new(L1A,'alt_sar_surf',0);
L1A = setfield_new(L1A,'x_sar_surf',0);
L1A = setfield_new(L1A,'y_sar_surf',0);
L1A = setfield_new(L1A,'z_sar_surf',0);
L1A = setfield_new(L1A,'doppler_ang_sar_sat',0);

L1A = setfield_new(L1A,'T0_sar',0);
L1A = setfield_new(L1A,'bri_sar',0);

%beam angles
L1A = setfield_new(L1A,'N_beams',0);
L1A = setfield_new(L1A,'surf_loc_index',0);
L1A = setfield_new(L1A,'beam_ang',0);
L1A = setfield_new(L1A,'beam_ang_nadir_index',0);
L1A = setfield_new(L1A,'beam_ang_index',0);
L1A = setfield_new(L1A,'start_beam',0);
L1A = setfield_new(L1A,'end_beam',0);
L1A = setfield_new(L1A,'beam_index',0);

%azimuth processing
L1A = setfield_new(L1A,'beams_focused_shifted',zeros(chd.N_pulses_burst,chd.N_samples_sar));
if(strcmp(cnf.processing_mode,'SIN'))
    L1A = setfield_new(L1A,'beams_focused_shifted_2',zeros(chd.N_pulses_burst,chd.N_samples_sar));
end


end

