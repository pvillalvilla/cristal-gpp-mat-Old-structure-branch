function [L1B]= reading_L1B_SIN (filesBulk, chd, cnf, cst)
    L1B.time                = double(ncread(filesBulk.filename_L1B, 'time_20_ku'));
    L1B.lat                 = double(ncread(filesBulk.filename_L1B, 'lat_20_ku'));
    L1B.lon                 = double(ncread(filesBulk.filename_L1B, 'lon_20_ku'));
    L1B.alt                 = double(ncread(filesBulk.filename_L1B, 'alt_20_ku'));
    L1B.vel                 = double(ncread(filesBulk.filename_L1B, 'sat_vel_vec_20_ku'));
    L1B.alt_rate            = double(ncread(filesBulk.filename_L1B, 'orb_alt_rate_20_ku'));
    L1B.roll                = double(ncread(filesBulk.filename_L1B, 'off_nadir_roll_angle_str_20_ku'));
    L1B.pitch               = double(ncread(filesBulk.filename_L1B, 'off_nadir_pitch_angle_str_20_ku'));
    L1B.yaw                 = double(ncread(filesBulk.filename_L1B, 'off_nadir_yaw_angle_str_20_ku'));
    L1B.echo_scale_factor   = double(ncread(filesBulk.filename_L1B, 'echo_scale_factor_20_ku'));
    L1B.echo_scale_power    = double(ncread(filename, 'echo_scale_pwr_20_ku'));
    L1B.waveform_scale_factor = double(ncread(filesBulk.filename_L1B, 'waveform_scale_factor_20_ku'));
    L1B.phase_diff = double(ncread(filesBulk.filename_L1B, 'ph_diff_waveform_20_ku')).';
    L1B.coherence  = double(ncread(filesBulk.filename_L1B, 'coherence_waveform_20_ku')).';
    L1B.look_ang_start      = double(ncread(filesBulk.filename_L1B, 'look_angle_start_20_ku'));
    L1B.look_ang_stop       = double(ncread(filesBulk.filename_L1B, 'look_angle_stop_20_ku'));
    L1B.doppler_angle_start = double(ncread(filesBulk.filename_L1B, 'dop_angle_start_20_ku'));
    L1B.doppler_angle_stop  = double(ncread(filesBulk.filename_L1B, 'dop_angle_stop_20_ku'));
    L1B.Doppler_mask        = double(ncread(filesBulk.filename_L1B, 'stack_mask_range_bin_20_ku'));
    L1B.Records_num         = double(ncread(filename, 'rec_count_20_ku'));
    L1B.window_delay        = double(ncread(filename, 'window_del_20_ku'))*1e-12;
    L1B.averaged_power_echo_waveform = double(ncread(filename, 'pwr_waveform_20_ku'));

end