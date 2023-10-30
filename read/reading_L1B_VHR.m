function [L1B]= reading_L1B_VHR(filesBulk, chd, cnf, cst)

L1B.time_ku = double(ncread(filesBulk.filename_L1B, 'time_ku'));

L1B.i2q2_meas_ku_l1b_echo       = double(ncread(filesBulk.filename_L1B, 'power_waveform'));
L1B.phase_diff_meas_ku_l1b_echo = double(ncread(filesBulk.filename_L1B, 'phase_waveform'));
L1B.coherence_meas_ku_l1b_echo  = double(ncread(filesBulk.filename_L1B, 'coherence_waveform'));

L1B.lat     = double(ncread(filesBulk.filename_L1B, 'latitude'));
L1B.lon     = double(ncread(filesBulk.filename_L1B, 'longitude'));
L1B.alt     = double(ncread(filesBulk.filename_L1B, 'altitude'));

L1B.roll    = double(ncread(filesBulk.filename_L1B, 'attitude_roll'));
L1B.pitch   = double(ncread(filesBulk.filename_L1B, 'attitude_pitch'));
L1B.yaw     = double(ncread(filesBulk.filename_L1B, 'attitude_yaw'));

L1B.range_ku_l1b_echo = double(ncread(filesBulk.filename_L1B, 'window_delay'));

L1B.number_of_looks = double(ncread(filesBulk.filename_L1B, 'number_of_looks'));
L1B.rg_sampling_step = double(ncread(filesBulk.filename_L1B, 'rg_sampling_step'));

L1B.scaled_waveforms = L1B.i2q2_meas_ku_l1b_echo;

end