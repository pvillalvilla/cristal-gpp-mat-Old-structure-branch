function [L1B_FF]= reading_L1B (filesBulk, chd, cnf, cst)
L1B_FF.time     = double(ncread(filesBulk.filename_L1BFF_SL, 'time_l1b_echo'));
L1B_FF.lat     = double(ncread(filesBulk.filename_L1BFF_SL, 'lat_l1b_echo'));
L1B_FF.lon     = double(ncread(filesBulk.filename_L1BFF_SL, 'lon_l1b_echo'));
L1B_FF.alt     = double(ncread(filesBulk.filename_L1BFF_SL, 'alt_l1b_echo'));
L1B_FF.x_vel   = double(ncread(filesBulk.filename_L1BFF_SL, 'x_vel_l1b_echo'));
L1B_FF.y_vel   = double(ncread(filesBulk.filename_L1BFF_SL, 'y_vel_l1b_echo'));
L1B_FF.z_vel   =  double(ncread(filesBulk.filename_L1BFF_SL, 'z_vel_l1b_echo'));
L1B_FF.alt_rate        = double(ncread(filesBulk.filename_L1BFF_SL, 'orb_alt_rate_l1b_echo'));
L1B_FF.pri        = double(ncread(filesBulk.filename_L1BFF_SL, 'pri_l1b_echo'));
L1B_FF.satellite_mispointing_l1b_echo  = double(ncread(filesBulk.filename_L1BFF_SL, 'satellite_mispointing_l1b_echo'));
L1B_FF.roll    = L1B_FF.satellite_mispointing_l1b_echo(1,:)*0;
L1B_FF.pitch   = L1B_FF.satellite_mispointing_l1b_echo(2,:)*0;
L1B_FF.yaw     = L1B_FF.satellite_mispointing_l1b_echo(3,:)*0;
L1B_FF.x_vel   = double(ncread(filesBulk.filename_L1BFF_SL, 'x_vel_l1b_echo'));
L1B_FF.y_vel   = double(ncread(filesBulk.filename_L1BFF_SL, 'y_vel_l1b_echo'));
L1B_FF.z_vel   = double(ncread(filesBulk.filename_L1BFF_SL, 'z_vel_l1b_echo'));
L1B_FF.vel     = sqrt(L1B_FF.x_vel.^2 + L1B_FF.y_vel.^ 2+ L1B_FF.z_vel.^2);
L1B_FF.nb_stack_l1b_echo               = double(ncread(filesBulk.filename_L1BFF_SL, 'nb_stack_l1b_echo'));
L1B_FF.sigma0           = double(ncread(filesBulk.filename_L1BFF_SL, 'scale_factor_ku_l1b_echo'));
L1B_FF.i2q2_meas_ku_l1b_echo           = double(ncread(filesBulk.filename_L1BFF_SL, 'i2q2_meas_ku_l1b_echo')).';
L1B_FF.waveform_scale_factor           = double(ncread(filesBulk.filename_L1BFF_SL, 'waveform_scale_factor_l1b_echo'));
if(strcmp(cnf.processing_mode,'SIN'))
    L1B_FF.i2q2_meas_ku_l1b_echo_2          = double(ncread(filesBulk.filename_L1BFF_SL, 'i2q2_meas_ku_l1b_echo_2')).';
    L1B_FF.phase_diff_meas_ku_l1b_echo = double(ncread(filesBulk.filename_L1BFF_SL, 'phase_diff_meas_ku_l1b_echo')).';
    L1B_FF.coherence_meas_ku_l1b_echo  = double(ncread(filesBulk.filename_L1BFF_SL, 'coherence_meas_ku_l1b_echo')).';
    L1B_FF.scaled_waveforms_2 = (((L1B_FF.waveform_scale_factor) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(L1B_FF.i2q2_meas_ku_l1b_echo_2).').'; 
end
L1B_FF.scaled_waveforms = (((L1B_FF.waveform_scale_factor) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(L1B_FF.i2q2_meas_ku_l1b_echo).').';
L1B_FF.range_ku_l1b_echo   = double(ncread(filesBulk.filename_L1BFF_SL, 'range_ku_l1b_echo'));
L1B_FF.look_ang_start      = double(ncread(filesBulk.filename_L1BFF_SL, 'look_angle_start_l1b_echo'));
L1B_FF.look_ang_stop       = double(ncread(filesBulk.filename_L1BFF_SL, 'look_angle_stop_l1b_echo'));
L1B_FF.doppler_angle_start = double(ncread(filesBulk.filename_L1BFF_SL, 'doppler_angle_start_l1b_echo'));
L1B_FF.doppler_angle_stop  = double(ncread(filesBulk.filename_L1BFF_SL, 'doppler_angle_stop_l1b_echo'));
L1B_FF.Doppler_mask        = double(ncread(filesBulk.filename_L1BFF_SL, 'stack_mask_range_bin_l1b_echo'));
%% FF-ML variables
%time variables
L1B_FF.time_surf_ML     = double(ncread(filesBulk.filename_L1BFF_ML, 'time_surf_ML'));
%orbit/attitude variables
L1B_FF.lat_surf_ML     = double(ncread(filesBulk.filename_L1BFF_ML, 'lat_surf_ML'));
L1B_FF.lon_surf_ML     = double(ncread(filesBulk.filename_L1BFF_ML, 'lon_surf_ML'));
L1B_FF.alt_surf_ML     = double(ncread(filesBulk.filename_L1BFF_ML, 'alt_surf_ML'));
L1B_FF.alt_rate_surf_ML        = double(ncread(filesBulk.filename_L1BFF_ML, 'alt_rate_surf_ML'));
%position/velocity variables
L1B_FF.x_vel_sat_surf_ML   = double(ncread(filesBulk.filename_L1BFF_ML, 'x_vel_sat_surf_ML'));
L1B_FF.y_vel_sat_surf_ML   = double(ncread(filesBulk.filename_L1BFF_ML, 'y_vel_sat_surf_ML'));
L1B_FF.z_vel_sat_surf_ML   =  double(ncread(filesBulk.filename_L1BFF_ML, 'z_vel_sat_surf_ML'));
%range variables
L1B_FF.range_ML   = double(ncread(filesBulk.filename_L1BFF_ML, 'range_ML'));
%waveforms variables
L1B_FF.power_waveform_ML   = double(ncread(filesBulk.filename_L1BFF_ML, 'power_waveform')).';
L1B_FF.waveform_scale_factor_ML   = double(ncread(filesBulk.filename_L1BFF_ML, 'waveform_scale_factor'));
L1B_FF.power_scaled_waveforms_ML = (((L1B_FF.waveform_scale_factor_ML) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(L1B_FF.power_waveform_ML).').';
if(strcmp(cnf.processing_mode,'SIN'))
     L1B_FF.power_waveform_ML_2   = double(ncread(filesBulk.filename_L1BFF_ML, 'power_waveform_2')).';
     L1B_FF.waveform_scale_factor_ML_2   = double(ncread(filesBulk.filename_L1BFF_ML, 'waveform_scale_factor_2'));
     L1B_FF.power_scaled_waveforms_ML_2 = (((L1B_FF.waveform_scale_factor_ML_2) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(L1B_FF.power_waveform_ML_2).').';
     L1B_FF.phase_difference_ML   = double(ncread(filesBulk.filename_L1BFF_ML, 'phase_diff_ML')).';
     L1B_FF.coherence_ML   = double(ncread(filesBulk.filename_L1BFF_ML, 'wvfm_coh_ML')).';
     L1B_FF.power_waveform_ML_comb   = double(ncread(filesBulk.filename_L1BFF_ML, 'power_waveform_comb')).';
     L1B_FF.waveform_scale_factor_ML_comb   = double(ncread(filesBulk.filename_L1BFF_ML, 'waveform_scale_factor_comb'));
     L1B_FF.power_scaled_waveforms_ML_comb = (((L1B_FF.waveform_scale_factor_ML_comb) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(L1B_FF.power_waveform_ML_comb).').';
end

end