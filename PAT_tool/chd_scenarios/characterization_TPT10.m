%% Point Target Location
chd.ECEF_X = 5.219993764835501e+06;
chd.ECEF_y = -2.104394012926102e+06;
chd.ECEF_Z = -2.990576451934797e+06;

chd.point_target_idx_ku = 62;

chd.right_beams_ku = 50;
chd.left_beams_ku = 50;

% Passar lat/lon/alt a ECEF
ecef_coord_ku = lla2ecef([L1B.lat_ku(chd.point_target_idx_ku), L1B.lon_ku(chd.point_target_idx_ku), ...
    L1B.alt_ku(chd.point_target_idx_ku)]);
ecef_x_sat_ku = ecef_coord_ku(1);
ecef_y_sat_ku = ecef_coord_ku(2);
ecef_z_sat_ku = ecef_coord_ku(3);

%% Expected values
% calcular distanies a cada punt del PT
chd.range_PT_ku = sqrt((chd.ECEF_X - ecef_x_sat_ku)^2 + (chd.ECEF_Y - ecef_y_sat_ku)^2 ...
    + (chd.ECEF_Z - ecef_z_sat_ku)^2);

%% Cst
chd.uso_freq_nom = 10e6;
chd.alt_freq_multiplier = 39.5;
chd.T0_nom = 1/(chd.uso_freq_nom*chd.alt_freq_multiplier);
chd.N_central_beams_ku = 20;
chd.N_central_beams_ka = 20;

%% Required Values
chd.range_bias_HR = 0.01;
chd.datation_bias = 100e-6; 
chd.range_bias_L1BS = 0.001;
chd.burst_air_PSLR = -13; %[dB]
chd.range_random_error = 1e-3;
chd.range_slope = 0.01e-3;
chd.side_lobes_ratio = 0.01; %db

