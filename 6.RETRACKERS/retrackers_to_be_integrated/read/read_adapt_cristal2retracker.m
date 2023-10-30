function [data, cst_p, chd_p] = read_adapt_cristal2retracker(L1B, cnf_L2, cst, chd) 





data.N_records = size(L1B.scaled_waveforms,1);
data.N_samples = chd.N_samples_sar;

data.GEO.LAT = L1B.lat;
data.GEO.LON = L1B.lon;
data.GEO.H   = L1B.alt;
data.GEO.H_rate = L1B.alt_rate;

data.GEO.roll  = L1B.roll;
data.GEO.pitch = L1B.pitch;
data.GEO.yaw   = L1B.yaw;

data.GEO.V = L1B.vel;

data.HRM.pri_surf = L1B.pri;
data.HRM.fs_clock_ku_surf = ones(size(L1B.scaled_waveforms,1),1).*chd.bw; %JPLZ: is this one correct?
data.HRM.Neff = L1B.nb_stack_l1b_echo;

% data.HRM.look_ang_stack = ;
data.HRM.look_ang_start_surf = L1B.look_ang_start;
data.HRM.look_ang_stop_surf = L1B.look_ang_stop;
data.HRM.doppler_angle_start_surf = L1B.doppler_angle_start;
data.HRM.doppler_angle_stop_surf = L1B.doppler_angle_stop;
%data.HRM.Doppler_mask = L1B.Doppler_mask.';
%data.HRM.Doppler_mask = L1B.Doppler_mask;
data.HRM.s0_sf = L1B.sigma0;




data.HRM.power_wav = L1B.scaled_waveforms.';

data.MEA.seed(1)= cnf_L2.seed;

data.MEA.win_delay = -0.1270.*ones(1,data.N_records);
try
data.MEA.win_delay = L1B.win_delay_l1b_echo;
end

% Convert cst to retracker input format
cst_p.semi_minor_axis_cst = cst.semi_minor_axis;
cst_p.mean_sat_alt_cst = cst.mean_sat_alt;
cst_p.sec_in_day_cst = cst.sec_in_day;
cst_p.c_cst  = cst.c;
cst_p.pi_cst = cst.pi;
cst_p.flat_coeff_cst = cst.flat_coeff;
cst_p.earth_radius_cst = cst.earth_radius;
cst_p.semi_major_axis_cst = cst.semi_major_axis;

% Convert chd to retracker input format
chd_p.N_samples_chd = chd.N_samples_sar;
chd_p.N_samples_sar_rmc_chd = chd.N_samples_sar_rmc;
chd_p.uso_freq_nom_chd = chd.uso_freq_nom;
chd_p.alt_freq_multiplier_chd = chd.alt_freq_multiplier;
chd_p.antenna_gain_chd = chd.antenna_gain;
chd_p.tx1_sar_chd = chd.tx1_sar;
chd_p.ext_delay_ground_chd = chd.ext_delay_ground;
chd_p.int_delay_ground_chd = chd.int_delay_ground;
chd_p.power_tx_ant_chd = chd.power_tx_ant;
chd_p.simulation_chd = chd.simulation;
chd_p.antenna_beamwidth_alt_ku_chd = chd.antenna_beamwidth_along_track;
chd_p.antenna_beamwidth_act_ku_chd = chd.antenna_beamwidth_across_track;
chd_p.N_cal_pulses_burst_chd = chd.N_cal_pulses_burst;
chd_p.N_c_pulses_burst_chd = chd.N_c_pulses_burst;
chd_p.pulse_length_chd = chd.pulse_length;
chd_p.bw_rx_ku_chd = chd.bw;
chd_p.alt_freq_chd = chd.alt_freq;
chd_p.T0_nom_chd = chd.T0_nom;
chd_p.start_freq_chd = chd.start_freq;
chd_p.phase_departure_chd = chd.phase_departure;
chd_p.meas_mode_chd = chd.meas_mode;
chd_p.N_total_pulses_b_chd = chd.N_ku_pulses_burst;
chd_p.N_bursts_cycle_chd = chd.N_bursts_cycle_sar;
chd_p.SWST_chd = chd.SWST;
chd_p.prf_chd = chd.prf;
chd_p.pri_nom_chd = chd.pri_nom;
chd_p.brf_chd = chd.brf;
chd_p.burst_duration_sar_chd = chd.burst_duration_sar;
chd_p.freq_ku_chd = chd.freq;
chd_p.wv_length_chd = chd.wv_length;
chd_p.wv_length_ku =  chd.wv_length; %JPLZ: has to be written this way?
chd_p.band_chd = chd.band;
chd_p.num_rx_chains_chd = chd.num_rx_chains;
chd_p.x_cog_ant_chd = chd.x_cog_ant;
chd_p.y_cog_ant_chd = chd.y_cog_ant;
chd_p.z_cog_ant_chd = chd.z_cog_ant;
chd_p.x_cog_ant_2_chd = chd.x_cog_ant_2;
chd_p.y_cog_ant_2_chd = chd.y_cog_ant_2;
chd_p. z_cog_ant_2_chd = chd. z_cog_ant_2;
try
chd_p.snow_height_chd = chd.snow_height;
chd_p.snow_depth_chd = chd.snow_depth;
end
chd_p.chirp_slope_chd = chd.chirp_slope;
chd_p.N_pulses_burst_chd = chd.N_pulses_burst;
chd_p.N_pri_chd = chd.N_pri;
chd_p.pri_T0_unit_conv_chd = chd.pri_T0_unit_conv;
chd_p.h0_cor2_unit_conv_chd = chd.h0_cor2_unit_conv;
chd_p.T0_h0_unit_conv_chd = chd.T0_h0_unit_conv;
chd_p.cai_cor2_unit_conv_chd = chd.cai_cor2_unit_conv;
chd_p.i_sample_start_chd = chd.i_sample_start;
chd_p.fai_shift_number_chd = chd.fai_shift_number;
chd_p.N_max_beams_stack_chd = chd.N_max_beams_stack;
chd_p.RMC_start_sample_chd = chd.RMC_start_sample;
chd_p.PTR_width_chd = chd.PTR_width;
chd_p.mask_look_angles_CR2_max_chd = chd.mask_look_angles_CR2_max;
chd_p.mask_look_angles_CR2_min_chd = chd.mask_look_angles_CR2_min;
chd_p.A_s2Ga_chd = chd.A_s2Ga;
chd_p.alpha_ga_chd = chd.alpha_ga;
chd_p.A_s2Gr_chd = chd.A_s2Gr;
chd_p.alpha_gr_chd = chd.alpha_gr;

chd_p.fs_clock_ku_chd = 1.*chd.bw; %JPLZ: is this one correct?

end