%% TIW.1 cut of the 3 different areas land-lake-land (according to Alessio Izzo mail 04/04/2022)
%% Ku band 
% ZONE 1
[chd.ku_start_lat(1), chd.ku_start_lat_idx(1)] = min(abs(-0.492928*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat(1), chd.ku_end_lat_idx(1)]     = min(abs(-0.492009*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon(1), chd.ku_start_lon_idx(1)] = min(abs(-0.383641*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon(1), chd.ku_end_lon_idx(1)]     = min(abs(-0.383746*180/pi - L1B_HR.lon_ku));

% ZONE 2
[chd.ku_start_lat(2), chd.ku_start_lat_idx(2)] = min(abs(-0.492009*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat(2), chd.ku_end_lat_idx(2)]     = min(abs(-0.485482*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon(2), chd.ku_start_lon_idx(2)] = min(abs(-0.383746*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon(2), chd.ku_end_lon_idx(2)]     = min(abs(-0.38449*180/pi - L1B_HR.lon_ku));

% ZONE 3   
[chd.ku_start_lat(3), chd.ku_start_lat_idx(3)] = min(abs(-0.485492*180/pi - L1B_HR.lat_ku));
[chd.ku_end_lat(3), chd.ku_end_lat_idx(3)]     = min(abs(-0.484559*180/pi - L1B_HR.lat_ku));

[chd.ku_start_lon(3), chd.ku_start_lon_idx(3)] = min(abs(-0.384489*180/pi - L1B_HR.lon_ku));
[chd.ku_end_lon(3), chd.ku_end_lon_idx(3)]     = min(abs(-0.384595*180/pi - L1B_HR.lon_ku));

%% Ka band 
% ZONE 1
[chd.ka_start_lat(1), chd.ka_start_lat_idx(1)] = min(abs(-0.492928*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat(1), chd.ka_end_lat_idx(1)]     = min(abs(-0.492009*180/pi - L1B_HR.lat_ka));

[chd.ka_start_lon(1), chd.ka_start_lon_idx(1)] = min(abs(-0.383641*180/pi - L1B_HR.lon_ka));
[chd.ka_end_lon(1), chd.ka_end_lon_idx(1)]     = min(abs(-0.383746*180/pi - L1B_HR.lon_ka));

% ZONE 2
[chd.ka_start_lat(2), chd.ka_start_lat_idx(2)] = min(abs(-0.492009*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat(2), chd.ka_end_lat_idx(2)]     = min(abs(-0.485482*180/pi - L1B_HR.lat_ka));

[chd.ka_start_lon(2), chd.ka_start_lon_idx(2)] = min(abs(-0.383746*180/pi - L1B_HR.lon_ka));
[chd.ka_end_lon(2), chd.ka_end_lon_idx(2)]     = min(abs(-0.38449*180/pi - L1B_HR.lon_ka));

% ZONE 3   
[chd.ka_start_lat(3), chd.ka_start_lat_idx(3)] = min(abs(-0.485492*180/pi - L1B_HR.lat_ka));
[chd.ka_end_lat(3), chd.ka_end_lat_idx(3)]     = min(abs(-0.484559*180/pi - L1B_HR.lat_ka));

[chd.ka_start_lon(3), chd.ka_start_lon_idx(3)] = min(abs(-0.384489*180/pi - L1B_HR.lon_ka));
[chd.ka_end_lon(3), chd.ka_end_lon_idx(3)]     = min(abs(-0.384595*180/pi - L1B_HR.lon_ka));

%% Expected values
chd.SSH(1) = 1;
chd.SSH(2) = 0;
chd.SSH(3) = 1;
%% Requiremets
chd.SSH_uncertainty = 0.1;
