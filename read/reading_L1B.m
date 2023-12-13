function [L1B]= reading_L1B (filesBulk, chd, cnf, cst)
L1B.time     = double(ncread(filesBulk.filename_L1B, 'time_l1b_echo'));
L1B.lat     = double(ncread(filesBulk.filename_L1B, 'lat_l1b_echo'));
L1B.lon     = double(ncread(filesBulk.filename_L1B, 'lon_l1b_echo'));
L1B.alt     = double(ncread(filesBulk.filename_L1B, 'alt_l1b_echo'));
L1B.x_vel   = double(ncread(filesBulk.filename_L1B, 'x_vel_l1b_echo'));
L1B.y_vel   = double(ncread(filesBulk.filename_L1B, 'y_vel_l1b_echo'));
L1B.z_vel   =  double(ncread(filesBulk.filename_L1B, 'z_vel_l1b_echo'));
L1B.alt_rate        = double(ncread(filesBulk.filename_L1B, 'orb_alt_rate_l1b_echo'));
L1B.pri        = double(ncread(filesBulk.filename_L1B, 'pri_l1b_echo'));
% L1B.satellite_mispointing_l1b_echo  = double(ncread(filesBulk.filename_L1B, 'satellite_mispointing_l1b_echo'));
L1B.roll    = 0; % L1B.satellite_mispointing_l1b_echo(1,:)*0;
L1B.pitch   = 0; % L1B.satellite_mispointing_l1b_echo(2,:)*0;
L1B.yaw     = 0; % L1B.satellite_mispointing_l1b_echo(3,:)*0;
L1B.x_vel   = double(ncread(filesBulk.filename_L1B, 'x_vel_l1b_echo'));
L1B.y_vel   = double(ncread(filesBulk.filename_L1B, 'y_vel_l1b_echo'));
L1B.z_vel   = double(ncread(filesBulk.filename_L1B, 'z_vel_l1b_echo'));
L1B.vel     = sqrt(L1B.x_vel.^2 + L1B.y_vel.^ 2+ L1B.z_vel.^2);
L1B.nb_stack_l1b_echo               = double(ncread(filesBulk.filename_L1B, 'nb_stack_l1b_echo'));
L1B.sigma0           = double(ncread(filesBulk.filename_L1B, 'scale_factor_ku_l1b_echo'));
L1B.i2q2_meas_ku_l1b_echo           = double(ncread(filesBulk.filename_L1B, 'i2q2_meas_ku_l1b_echo')).';
L1B.waveform_scale_factor           = double(ncread(filesBulk.filename_L1B, 'waveform_scale_factor_l1b_echo'));
if(strcmp(cnf.processing_mode,'SIN'))
    L1B.phase_diff_meas_ku_l1b_echo = double(ncread(filesBulk.filename_L1B, 'phase_diff_meas_ku_l1b_echo')).';
    L1B.coherence_meas_ku_l1b_echo  = double(ncread(filesBulk.filename_L1B, 'coherence_meas_ku_l1b_echo')).';
end
L1B.scaled_waveforms = (((L1B.waveform_scale_factor) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(L1B.i2q2_meas_ku_l1b_echo).').';
L1B.range_ku_l1b_echo   = double(ncread(filesBulk.filename_L1B, 'range_ku_l1b_echo'));
try
L1B.win_delay_l1b_echo   = double(ncread(filesBulk.filename_L1B, 'win_delay_l1b_echo'));
end
L1B.look_ang_start      = double(ncread(filesBulk.filename_L1B, 'look_angle_start_l1b_echo'));
L1B.look_ang_stop       = double(ncread(filesBulk.filename_L1B, 'look_angle_stop_l1b_echo'));
L1B.doppler_angle_start = double(ncread(filesBulk.filename_L1B, 'doppler_angle_start_l1b_echo'));
L1B.doppler_angle_stop  = double(ncread(filesBulk.filename_L1B, 'doppler_angle_stop_l1b_echo'));
L1B.Doppler_mask        = double(ncread(filesBulk.filename_L1B, 'stack_mask_range_bin_l1b_echo'));
L1B.stdev_stack_l1b_echo= double(ncread(filesBulk.filename_L1B, 'stdev_stack_l1b_echo'));
end