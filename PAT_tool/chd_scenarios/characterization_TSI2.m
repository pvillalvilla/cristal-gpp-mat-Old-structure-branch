%% Ku band latitudes/longitudes indexs (L1B_HR)
% ZONE 1
[chd.ku_start_lat_HR(1), chd.ku_start_lat_idx_HR(1)] = min(abs(-0.49489*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat_HR(1), chd.ku_end_lat_idx_HR(1)]     = min(abs(-0.492378*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon_HR(1), chd.ku_start_lon_idx_HR(1)] = min(abs(-0.38339*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon_HR(1), chd.ku_end_lon_idx_HR(1)]     = min(abs(-0.383704*180/pi - L1B_HR.lon_ku));

% ZONE 2
[chd.ku_start_lat_HR(2), chd.ku_start_lat_idx_HR(2)] = min(abs(-0.492378*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat_HR(2), chd.ku_end_lat_idx_HR(2)]     = min(abs(-0.492376*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon_HR(2), chd.ku_start_lon_idx_HR(2)] = min(abs(-0.383704*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon_HR(2), chd.ku_end_lon_idx_HR(2)]     = min(abs(-0.3837044*180/pi - L1B_HR.lon_ku));

% ZONE 3   
[chd.ku_start_lat_HR(3), chd.ku_start_lat_idx_HR(3)] = min(abs(-0.492376*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat_HR(3), chd.ku_end_lat_idx_HR(3)]     = min(abs(-0.48986*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon_HR(3), chd.ku_start_lon_idx_HR(3)] = min(abs(-0.3837044*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon_HR(3), chd.ku_end_lon_idx_HR(3)]     = min(abs(-0.383988*180/pi - L1B_HR.lon_ku));

%% Ka band latitudes/longitudes indexs (L1B_HR)
% ZONE 1
[chd.ka_start_lat_HR(1), chd.ka_start_lat_idx_HR(1)] = min(abs(-0.49489*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat_HR(1), chd.ka_end_lat_idx_HR(1)]     = min(abs(-0.492378*180/pi - L1B_HR.lat_ka));

[chd.ka_start_lon_HR(1), chd.ka_start_lon_idx_HR(1)] = min(abs(-0.38339*180/pi - L1B_HR.lon_ka));
[chd.ka_end_lon_HR(1), chd.ka_end_lon_idx_HR(1)]     = min(abs(-0.383704*180/pi - L1B_HR.lon_ka));

% ZONE 2
[chd.ka_start_lat_HR(2), chd.ka_start_lat_idx_HR(2)] = min(abs(-0.492378*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat_HR(2), chd.ka_end_lat_idx_HR(2)]     = min(abs(-0.492376*180/pi - L1B_HR.lat_ka));

[chd.ka_start_lon_HR(2), chd.ka_start_lon_idx_HR(2)] = min(abs(-0.383704*180/pi - L1B_HR.lon_ka));
[chd.ka_end_lon_HR(2), chd.ka_end_lon_idx_HR(2)]     = min(abs(-0.3837044*180/pi - L1B_HR.lon_ka));

% ZONE 3   
[chd.ka_start_lat_HR(3), chd.ka_start_lat_idx_HR(3)] = min(abs(-0.492376*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat_HR(3), chd.ka_end_lat_idx_HR(3)]     = min(abs(-0.48986*180/pi - L1B_HR.lat_ka));

[chd.ka_start_lon_HR(3), chd.ka_start_lon_idx_HR(3)] = min(abs(-0.3837044*180/pi - L1B_HR.lon_ka));
[chd.ka_end_lon_HR(3), chd.ka_end_lon_idx_HR(3)]     = min(abs(-0.383988*180/pi - L1B_HR.lon_ka));


%% Ku band latitudes/longitudes indexs (L1B_ML)
% ZONE 1
[chd.ku_start_lat_ML_FF(1), chd.ku_start_lat_idx_ML_FF(1)] = min(abs(-0.49489*180/pi - L1B_ML.lat_ku));
[chd.ku_end_lat_ML_FF(1), chd.ku_end_lat_idx_ML_FF(1)]     = min(abs(-0.492378*180/pi - L1B_ML.lat_ku));

[chd.ku_start_lon_ML_FF(1), chd.ku_start_lon_idx_ML_FF(1)] = min(abs(-0.38339*180/pi - L1B_ML.lon_ku));
[chd.ku_end_lon_ML_FF(1), chd.ku_end_lon_idx_ML_FF(1)]     = min(abs(-0.383704*180/pi - L1B_ML.lon_ku));

% ZONE 2
[chd.ku_start_lat_ML_FF(2), chd.ku_start_lat_idx_ML_FF(2)] = min(abs(-0.492378*180/pi - L1B_ML.lat_ku));
[chd.ku_end_lat_ML_FF(2), chd.ku_end_lat_idx_ML_FF(2)]     = min(abs(-0.492376*180/pi - L1B_ML.lat_ku));

[chd.ku_start_lon_ML_FF(2), chd.ku_start_lon_idx_ML_FF(2)] = min(abs(-0.383704*180/pi - L1B_ML.lon_ku));
[chd.ku_end_lon_ML_FF(2), chd.ku_end_lon_idx_ML_FF(2)]     = min(abs(-0.3837044*180/pi - L1B_ML.lon_ku));

% ZONE 3   
[chd.ku_start_lat_ML_FF(3), chd.ku_start_lat_idx_ML_FF(3)] = min(abs(-0.492376*180/pi - L1B_ML.lat_ku));
[chd.ku_end_lat_ML_FF(3), chd.ku_end_lat_idx_ML_FF(3)]     = min(abs(-0.48986*180/pi - L1B_ML.lat_ku));

[chd.ku_start_lon_ML_FF(3), chd.ku_start_lon_idx_ML_FF(3)] = min(abs(-0.3837044*180/pi - L1B_ML.lon_ku));
[chd.ku_end_lon_ML_FF(3), chd.ku_end_lon_idx_ML_FF(3)]     = min(abs(-0.383988*180/pi - L1B_ML.lon_ku));

%% Ka band latitudes/longitudes indexs (L1B_ML)
% ZONE 1
[chd.ka_start_lat_ML_FF(1), chd.ka_start_lat_idx_ML_FF(1)] = min(abs(-0.49489*180/pi - L1B_ML.lat_ka));
[chd.ka_end_lat_ML_FF(1), chd.ka_end_lat_idx_ML_FF(1)]     = min(abs(-0.492378*180/pi - L1B_ML.lat_ka));

[chd.ka_start_lon_ML_FF(1), chd.ka_start_lon_idx_ML_FF(1)] = min(abs(-0.38339*180/pi - L1B_ML.lon_ka));
[chd.ka_end_lon_ML_FF(1), chd.ka_end_lon_idx_ML_FF(1)]     = min(abs(-0.383704*180/pi - L1B_ML.lon_ka));

% ZONE 2
[chd.ka_start_lat_ML_FF(2), chd.ka_start_lat_idx_ML_FF(2)] = min(abs(-0.492378*180/pi - L1B_ML.lat_ka));
[chd.ka_end_lat_ML_FF(2), chd.ka_end_lat_idx_ML_FF(2)]     = min(abs(-0.492376*180/pi - L1B_ML.lat_ka));

[chd.ka_start_lon_ML_FF(2), chd.ka_start_lon_idx_ML_FF(2)] = min(abs(-0.383704*180/pi - L1B_ML.lon_ka));
[chd.ka_end_lon_ML_FF(2), chd.ka_end_lon_idx_ML_FF(2)]     = min(abs(-0.3837044*180/pi - L1B_ML.lon_ka));

% ZONE 3   
[chd.ka_start_lat_ML_FF(3), chd.ka_start_lat_idx_ML_FF(3)] = min(abs(-0.492376*180/pi - L1B_ML.lat_ka));
[chd.ka_end_lat_ML_FF(3), chd.ka_end_lat_idx_ML_FF(3)]     = min(abs(-0.48986*180/pi - L1B_ML.lat_ka));

[chd.ka_start_lon_ML_FF(3), chd.ka_start_lon_idx_ML_FF(3)] = min(abs(-0.3837044*180/pi - L1B_ML.lon_ka));
[chd.ka_end_lon_ML_FF(3), chd.ka_end_lon_idx_ML_FF(3)]     = min(abs(-0.383988*180/pi - L1B_ML.lon_ka));
%% Expected values 
chd.height_ssh_ku(1) = 1.5 + 0.2; % Ice
chd.height_ssh_ka(1) = 1.5 + 0.2 + 0.5; %Ice + Snow
chd.height_ssh_ku(2) = 1.5 + 0.2; % Ice
chd.height_ssh_ka(2) = 1.5 + 0.2 + 0.5; %Ice + Snow
%% Requirements
chd.random_range_error = 0.08;
