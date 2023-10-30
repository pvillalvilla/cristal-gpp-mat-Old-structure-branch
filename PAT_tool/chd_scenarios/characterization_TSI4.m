%% Ku band latitudes/longitudes indexs (L1B_HR)--------------fets amb coord correctes
% ZONE 1
[chd.ku_start_lat_HR(1), chd.ku_start_lat_idx_HR(1)] = min(abs(-0.4936318*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat_HR(1), chd.ku_end_lat_idx_HR(1)]     = min(abs(-0.49237151*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon_HR(1), chd.ku_start_lon_idx_HR(1)] = min(abs(-0.38356121*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon_HR(1), chd.ku_end_lon_idx_HR(1)]     = min(abs(-0.3837050*180/pi - L1B_HR.lon_ku));

% ZONE 2 LEAD
[chd.ku_start_lat_HR(2), chd.ku_start_lat_idx_HR(2)] = min(abs(-0.49237151*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat_HR(2), chd.ku_end_lat_idx_HR(2)]     = min(abs(-0.49236776*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon_HR(2), chd.ku_start_lon_idx_HR(2)] = min(abs(-0.3837050*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon_HR(2), chd.ku_end_lon_idx_HR(2)]     = min(abs(-0.3837054*180/pi - L1B_HR.lon_ku));

% ZONE 3   
[chd.ku_start_lat_HR(3), chd.ku_start_lat_idx_HR(3)] = min(abs(-0.49236776*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat_HR(3), chd.ku_end_lat_idx_HR(3)]     = min(abs(-0.4909818*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon_HR(3), chd.ku_start_lon_idx_HR(3)] = min(abs(-0.3837054*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon_HR(3), chd.ku_end_lon_idx_HR(3)]     = min(abs(-0.38386352*180/pi - L1B_HR.lon_ku));

%% Ka band latitudes/longitudes indexs (L1B_HR)----------fet
% ZONE 1
[chd.ka_start_lat_HR(1), chd.ka_start_lat_idx_HR(1)] = min(abs(-0.4936318*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat_HR(1), chd.ka_end_lat_idx_HR(1)]     = min(abs(-0.49237151*180/pi - L1B_HR.lat_ka));

[chd.ka_start_lon_HR(1), chd.ka_start_lon_idx_HR(1)] = min(abs(-0.38356121*180/pi - L1B_HR.lon_ka));
[chd.ka_end_lon_HR(1), chd.ka_end_lon_idx_HR(1)]     = min(abs(-0.3837050*180/pi - L1B_HR.lon_ka));

% ZONE 2 LEAD
[chd.ka_start_lat_HR(2), chd.ka_start_lat_idx_HR(2)] = min(abs(-0.49237151*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat_HR(2), chd.ka_end_lat_idx_HR(2)]     = min(abs(-0.49236776*180/pi - L1B_HR.lat_ka));

[chd.ka_start_lon_HR(2), chd.ka_start_lon_idx_HR(2)] = min(abs(-0.3837050*180/pi - L1B_HR.lon_ka));
[chd.ka_end_lon_HR(2), chd.ka_end_lon_idx_HR(2)]     = min(abs(-0.3837054*180/pi - L1B_HR.lon_ka));

% ZONE 3   
[chd.ka_start_lat_HR(3), chd.ka_start_lat_idx_HR(3)] = min(abs(-0.49236776*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat_HR(3), chd.ka_end_lat_idx_HR(3)]     = min(abs(-0.4909818*180/pi - L1B_HR.lat_ka));

[chd.ku_start_lon_HR(3), chd.ku_start_lon_idx_HR(3)] = min(abs(-0.3837054*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon_HR(3), chd.ku_end_lon_idx_HR(3)]     = min(abs(-0.38386352*180/pi - L1B_HR.lon_ku));
%% Ku band latitudes/longitudes indexs (L1B_ML)--------------fets amb coord correctes
% ZONE 1
[chd.ku_start_lat_ML_FF(1), chd.ku_start_lat_idx_ML_FF(1)] = min(abs(-0.4936318*180/pi - L1B_ML.lat_ku));
[chd.ku_end_lat_ML_FF(1), chd.ku_end_lat_idx_ML_FF(1)]     = min(abs(-0.49237151*180/pi - L1B_ML.lat_ku));

[chd.ku_start_lon_ML_FF(1), chd.ku_start_lon_idx_ML_FF(1)] = min(abs(-0.38356121*180/pi - L1B_ML.lon_ku));
[chd.ku_end_lon_ML_FF(1), chd.ku_end_lon_idx_ML_FF(1)]     = min(abs(-0.3837050*180/pi - L1B_ML.lon_ku));

% ZONE 2
[chd.ku_start_lat_ML_FF(2), chd.ku_start_lat_idx_ML_FF(2)] = min(abs(-0.49237151*180/pi - L1B_ML.lat_ku));
[chd.ku_end_lat_ML_FF(2), chd.ku_end_lat_idx_ML_FF(2)]     = min(abs(-0.49236776*180/pi - L1B_ML.lat_ku));

[chd.ku_start_lon_ML_FF(2), chd.ku_start_lon_idx_ML_FF(2)] = min(abs(-0.3837050*180/pi - L1B_ML.lon_ku));
[chd.ku_end_lon_ML_FF(2), chd.ku_end_lon_idx_ML_FF(2)]     = min(abs(-0.3837054*180/pi - L1B_ML.lon_ku));

% ZONE 3   
[chd.ku_start_lat_ML_FF(3), chd.ku_start_lat_idx_ML_FF(3)] = min(abs(-0.49236776*180/pi - L1B_ML.lat_ku));
[chd.ku_end_lat_ML_FF(3), chd.ku_end_lat_idx_ML_FF(3)]     = min(abs(-0.4909818*180/pi - L1B_ML.lat_ku));

[chd.ku_start_lon_ML_FF(3), chd.ku_start_lon_idx_ML_FF(3)] = min(abs(-0.3837054*180/pi - L1B_ML.lon_ku));
[chd.ku_end_lon_ML_FF(3), chd.ku_end_lon_idx_ML_FF(3)]     = min(abs(-0.38386352*180/pi - L1B_ML.lon_ku));

%% Ka band latitudes/longitudes indexs (L1B_ML)
% ZONE 1
[chd.ka_start_lat_ML_FF(1), chd.ka_start_lat_idx_ML_FF(1)] = min(abs(-0.4936318*180/pi - L1B_ML.lat_ka));
[chd.ka_end_lat_ML_FF(1), chd.ka_end_lat_idx_ML_FF(1)]     = min(abs(-0.49237151*180/pi - L1B_ML.lat_ka));

[chd.ka_start_lon_ML_FF(1), chd.ka_start_lon_idx_ML_FF(1)] = min(abs(-0.38356121*180/pi - L1B_ML.lon_ka));
[chd.ka_end_lon_ML_FF(1), chd.ka_end_lon_idx_ML_FF(1)]     = min(abs(-0.3837050*180/pi - L1B_ML.lon_ka));

% ZONE 2
[chd.ka_start_lat_ML_FF(2), chd.ka_start_lat_idx_ML_FF(2)] = min(abs(-0.49237151*180/pi - L1B_ML.lat_ka));
[chd.ka_end_lat_ML_FF(2), chd.ka_end_lat_idx_ML_FF(2)]     = min(abs(-0.49236776*180/pi - L1B_ML.lat_ka));

[chd.ka_start_lon_ML_FF(2), chd.ka_start_lon_idx_ML_FF(2)] = min(abs(-0.3837050*180/pi - L1B_ML.lon_ka));
[chd.ka_end_lon_ML_FF(2), chd.ka_end_lon_idx_ML_FF(2)]     = min(abs(-0.3837054*180/pi - L1B_ML.lon_ka));

% ZONE 3   
[chd.ka_start_lat_ML_FF(3), chd.ka_start_lat_idx_ML_FF(3)] = min(abs(-0.49236776*180/pi - L1B_ML.lat_ka));
[chd.ka_end_lat_ML_FF(3), chd.ka_end_lat_idx_ML_FF(3)]     = min(abs(-0.4909818*180/pi - L1B_ML.lat_ka));

[chd.ka_start_lon_ML_FF(3), chd.ka_start_lon_idx_ML_FF(3)] = min(abs(-0.3837054*180/pi - L1B_ML.lon_ka));
[chd.ka_end_lon_ML_FF(3), chd.ka_end_lon_idx_ML_FF(3)]     = min(abs(-0.38386352*180/pi - L1B_ML.lon_ka));

%% Expected values 
chd.height_ssh_ku(1) = 1.5; % Reference height
chd.height_ssh_ka(1) = 1.5 + 0.5; % Reference height + Snow
chd.height_ssh_ku(2) = 1.5; 
chd.height_ssh_ka(2) = 1.5 + 0.5; % Snow
%% Requirements
chd.random_range_error = 0.08;
