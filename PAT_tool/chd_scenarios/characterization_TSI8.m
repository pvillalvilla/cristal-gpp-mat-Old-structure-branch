
%% TSI.8 cut of the 4 different snow depth areas (according to Alessio Izzo mail 11/01/2022)
%% Ku band latitudes/longitudes indexs (L1B_HR)
% ZONE 1
[chd.ku_start_lat_HR(1), chd.ku_start_lat_idx_HR(1)] = min(abs(-0.494381*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat_HR(1), chd.ku_end_lat_idx_HR(1)]     = min(abs(-0.491598*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon_HR(1), chd.ku_start_lon_idx_HR(1)] = min(abs(-0.383444*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon_HR(1), chd.ku_end_lon_idx_HR(1)]     = min(abs(-0.383793*180/pi - L1B_HR.lon_ku));

% ZONE 2
[chd.ku_start_lat_HR(2), chd.ku_start_lat_idx_HR(2)] = min(abs(-0.491649*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat_HR(2), chd.ku_end_lat_idx_HR(2)]     = min(abs(-0.490536*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon_HR(2), chd.ku_start_lon_idx_HR(2)] = min(abs(-0.383787*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon_HR(2), chd.ku_end_lon_idx_HR(2)]     = min(abs(-0.383914*180/pi - L1B_HR.lon_ku));

% ZONE 3   
[chd.ku_start_lat_HR(3), chd.ku_start_lat_HR_idx_HR(3)] = min(abs(-0.490605*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat_HR(3), chd.ku_end_lat_HR_idx_HR(3)]     = min(abs(-0.489498*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon_HR(3), chd.ku_start_lon_idx_HR(3)] = min(abs(-0.383906*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon_HR(3), chd.ku_end_lon_idx_HR(3)]     = min(abs(-0.384033*180/pi - L1B_HR.lon_ku));

% ZONE 4   
[chd.ku_start_lat_HR(4), chd.ku_start_lat_idx_HR(4)] = min(abs(-0.489546*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat_HR(4), chd.ku_end_lat_idx_HR(4)]     = min(abs(-0.486803*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon_HR(4), chd.ku_start_lon_idx_HR(4)] = min(abs(-0.384027*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon_HR(4), chd.ku_end_lon_idx_HR(4)]     = min(abs(-0.384339*180/pi - L1B_HR.lon_ku));

%% Ka band latitudes/longitudes indexs (L1B_HR)
    % ZONE 1
[chd.ka_start_lat_HR(1), chd.ka_start_lat_idx_HR(1)] = min(abs(-0.494381*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat_HR(1), chd.ka_end_lat_idx_HR(1)]     = min(abs(-0.491598*180/pi - L1B_HR.lat_ka));

[chd.ka_start_lon_HR(1), chd.ka_start_lon_idx_HR(1)] = min(abs(-0.383444*180/pi - L1B_HR.lon_ka));
[chd.ka_end_lon_HR(1), chd.ka_end_lon_idx_HR(1)]     = min(abs(-0.383793*180/pi - L1B_HR.lon_ka));

% ZONE 2
[chd.ka_start_lat_HR(2), chd.ka_start_lat_idx_HR(2)] = min(abs(-0.491649*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat_HR(2), chd.ka_end_lat_idx_HR(2)]     = min(abs(-0.490536*180/pi - L1B_HR.lat_ka));

[chd.ka_start_lon_HR(2), chd.ka_start_lon_idx_HR(2)] = min(abs(-0.383787*180/pi - L1B_HR.lon_ka));
[chd.ka_end_lon_HR(2), chd.ka_end_lon_idx_HR(2)]     = min(abs(-0.383914*180/pi - L1B_HR.lon_ka));

% ZONE 3   
[chd.ka_start_lat_HR(3), chd.ka_start_lat_idx_HR(3)] = min(abs(-0.490605*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat_HR(3), chd.ka_end_lat_idx_HR(3)]     = min(abs(-0.489498*180/pi - L1B_HR.lat_ka));

[chd.ka_start_lon_HR(3), chd.ka_start_lon_idx_HR(3)] = min(abs(-0.383906*180/pi - L1B_HR.lon_ka));
[chd.ka_end_lon_HR(3), chd.ka_end_lon_idx_HR(3)]     = min(abs(-0.384033*180/pi - L1B_HR.lon_ka));

% ZONE 4   
[chd.ka_start_lat_HR(4), chd.ka_start_lat_idx_HR(4)] = min(abs(-0.489546*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat_HR(4), chd.ka_end_lat_idx_HR(4)]     = min(abs(-0.486803*180/pi - L1B_HR.lat_ka));

[chd.ka_start_lon_HR(4), chd.ka_start_lon_idx_HR(4)] = min(abs(-0.384027*180/pi - L1B_HR.lon_ka));
[chd.ka_end_lon_HR(4), chd.ka_end_lon_idx_HR(4)]     = min(abs(-0.384339*180/pi - L1B_HR.lon_ka));

 %% Ku band latitudes/longitudes indexs (L1B_ML)
% ZONE 1
[chd.ku_start_lat_ML_FF(1), chd.ku_start_lat_idx_ML_FF(1)] = min(abs(-0.494381*180/pi - L1B_ML.lat_ku));
[chd.ku_end_lat_ML_FF(1), chd.ku_end_lat_idx_ML_FF(1)]     = min(abs(-0.491598*180/pi - L1B_ML.lat_ku));

[chd.ku_start_lon_ML_FF(1), chd.ku_start_lon_idx_ML_FF(1)] = min(abs(-0.383444*180/pi - L1B_ML.lon_ku));
[chd.ku_end_lon_ML_FF(1), chd.ku_end_lon_idx_ML_FF(1)]     = min(abs(-0.383793*180/pi - L1B_ML.lon_ku));

% ZONE 2
[chd.ku_start_lat_ML_FF(2), chd.ku_start_lat_idx_ML_FF(2)] = min(abs(-0.491649*180/pi - L1B_ML.lat_ku));
[chd.ku_end_lat_ML_FF(2), chd.ku_end_lat_idx_ML_FF(2)]     = min(abs(-0.490536*180/pi - L1B_ML.lat_ku));

[chd.ku_start_lon_ML_FF(2), chd.ku_start_lon_idx_ML_FF(2)] = min(abs(-0.383787*180/pi - L1B_ML.lon_ku));
[chd.ku_end_lon_ML_FF(2), chd.ku_end_lon_idx_ML_FF(2)]     = min(abs(-0.383914*180/pi - L1B_ML.lon_ku));

% ZONE 3   
[chd.ku_start_lat_ML_FF(3), chd.ku_start_lat_ML_FF_idx_ML_FF(3)] = min(abs(-0.490605*180/pi - L1B_ML.lat_ku));
[chd.ku_end_lat_ML_FF(3), chd.ku_end_lat_ML_FF_idx_ML_FF(3)]     = min(abs(-0.489498*180/pi - L1B_ML.lat_ku));

[chd.ku_start_lon_ML_FF(3), chd.ku_start_lon_idx_ML_FF(3)] = min(abs(-0.383906*180/pi - L1B_ML.lon_ku));
[chd.ku_end_lon_ML_FF(3), chd.ku_end_lon_idx_ML_FF(3)]     = min(abs(-0.384033*180/pi - L1B_ML.lon_ku));

% ZONE 4   
[chd.ku_start_lat_ML_FF(4), chd.ku_start_lat_idx_ML_FF(4)] = min(abs(-0.489546*180/pi - L1B_ML.lat_ku));
[chd.ku_end_lat_ML_FF(4), chd.ku_end_lat_idx_ML_FF(4)]     = min(abs(-0.486803*180/pi - L1B_ML.lat_ku));

[chd.ku_start_lon_ML_FF(4), chd.ku_start_lon_idx_ML_FF(4)] = min(abs(-0.384027*180/pi - L1B_ML.lon_ku));
[chd.ku_end_lon_ML_FF(4), chd.ku_end_lon_idx_ML_FF(4)]     = min(abs(-0.384339*180/pi - L1B_ML.lon_ku));

%% Ka band latitudes/longitudes indexs (L1B_ML)
    % ZONE 1
[chd.ka_start_lat_ML_FF(1), chd.ka_start_lat_idx_ML_FF(1)] = min(abs(-0.494381*180/pi - L1B_ML.lat_ka));
[chd.ka_end_lat_ML_FF(1), chd.ka_end_lat_idx_ML_FF(1)]     = min(abs(-0.491598*180/pi - L1B_ML.lat_ka));

[chd.ka_start_lon_ML_FF(1), chd.ka_start_lon_idx_ML_FF(1)] = min(abs(-0.383444*180/pi - L1B_ML.lon_ka));
[chd.ka_end_lon_ML_FF(1), chd.ka_end_lon_idx_ML_FF(1)]     = min(abs(-0.383793*180/pi - L1B_ML.lon_ka));

% ZONE 2
[chd.ka_start_lat_ML_FF(2), chd.ka_start_lat_idx_ML_FF(2)] = min(abs(-0.491649*180/pi - L1B_ML.lat_ka));
[chd.ka_end_lat_ML_FF(2), chd.ka_end_lat_idx_ML_FF(2)]     = min(abs(-0.490536*180/pi - L1B_ML.lat_ka));

[chd.ka_start_lon_ML_FF(2), chd.ka_start_lon_idx_ML_FF(2)] = min(abs(-0.383787*180/pi - L1B_ML.lon_ka));
[chd.ka_end_lon_ML_FF(2), chd.ka_end_lon_idx_ML_FF(2)]     = min(abs(-0.383914*180/pi - L1B_ML.lon_ka));

% ZONE 3   
[chd.ka_start_lat_ML_FF(3), chd.ka_start_lat_idx_ML_FF(3)] = min(abs(-0.490605*180/pi - L1B_ML.lat_ka));
[chd.ka_end_lat_ML_FF(3), chd.ka_end_lat_idx_ML_FF(3)]     = min(abs(-0.489498*180/pi - L1B_ML.lat_ka));

[chd.ka_start_lon_ML_FF(3), chd.ka_start_lon_idx_ML_FF(3)] = min(abs(-0.383906*180/pi - L1B_ML.lon_ka));
[chd.ka_end_lon_ML_FF(3), chd.ka_end_lon_idx_ML_FF(3)]     = min(abs(-0.384033*180/pi - L1B_ML.lon_ka));

% ZONE 4   
[chd.ka_start_lat_ML_FF(4), chd.ka_start_lat_idx_ML_FF(4)] = min(abs(-0.489546*180/pi - L1B_ML.lat_ka));
[chd.ka_end_lat_ML_FF(4), chd.ka_end_lat_idx_ML_FF(4)]     = min(abs(-0.486803*180/pi - L1B_ML.lat_ka));

[chd.ka_start_lon_ML_FF(4), chd.ka_start_lon_idx_ML_FF(4)] = min(abs(-0.384027*180/pi - L1B_ML.lon_ka));
[chd.ka_end_lon_ML_FF(4), chd.ka_end_lon_idx_ML_FF(4)]     = min(abs(-0.384339*180/pi - L1B_ML.lon_ka));
%% Expected values

chd.snow_depth(1) = 5e-2;
chd.snow_depth(2) = 2e-1;
chd.snow_depth(3) = 3.5e-1;
chd.snow_depth(4) = 5e-1;

%% Requirements
chd.snow_depth_uncertainty = 0.05;
chd.snow_depth_bias = 0; % Needs to be changed



