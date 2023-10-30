%% Ku band 
    % ZONE 1
    [chd.ku_start_lat(1), chd.ku_start_lat_idx(1)] = min(abs(1.44477*180/pi - GR.lat_ku));
    [chd.ku_end_lat(1), chd.ku_end_lat_idx(1)]     = min(abs(1.45066*180/pi - GR.lat_ku));
    
    [chd.ku_start_lon(1), chd.ku_start_lon_idx(1)] = min(abs(-0.37787*180/pi - GR.lon_ku));
    [chd.ku_end_lon(1), chd.ku_end_lon_HR_idx(1)]     = min(abs(-0.39278*180/pi - GR.lon_ku));

    % ZONE 2
    [chd.ku_start_lat(2), chd.ku_start_lat_idx(2)] = min(abs(1.45066*180/pi - GR.lat_ku));
    [chd.ku_end_lat(2), chd.ku_end_lat_idx(2)]     = min(abs(1.45076*180/pi - GR.lat_ku));

    [chd.ku_start_lon(2), chd.ku_start_lon_idx(2)] = min(abs(-0.39278*180/pi - GR.lon_ku));
    [chd.ku_end_lon(2), chd.ku_end_lon_idx(2)]     = min(abs(-0.39304*180/pi - GR.lon_ku));

    % ZONE 3   
    [chd.ku_start_lat(3), chd.ku_start_lat_idx(3)] = min(abs(1.45076*180/pi - GR.lat_ku));
    [chd.ku_end_lat(3), chd.ku_end_lat_idx(3)]     = min(abs(1.456715*180/pi - GR.lat_ku));

    [chd.ku_start_lon(3), chd.ku_start_lon_idx(3)] = min(abs(-0.39304*180/pi - GR.lon_ku));
    [chd.ku_end_lon(3), chd.ku_end_lon_idx(3)]     = min(abs(-0.40856*180/pi - GR.lon_ku));
    
    %% Ka band 
    % ZONE 1
    [chd.ka_start_lat(1), chd.ka_start_lat_idx(1)] = min(abs(1.44477*180/pi - GR.lat_ka));
    [chd.ka_end_lat(1), chd.ka_end_lat_idx(1)]     = min(abs(1.45066*180/pi - GR.lat_ka));
    
    [chd.ka_start_lon(1), chd.ka_start_lon_idx(1)] = min(abs(-0.37787*180/pi - GR.lon_ka));
    [chd.ka_end_lon(1), chd.ka_end_lon_HR_idx(1)]     = min(abs(-0.39278*180/pi - GR.lon_ka));

    % ZONE 2
    [chd.ka_start_lat(2), chd.ka_start_lat_idx(2)] = min(abs(1.45066*180/pi - GR.lat_ka));
    [chd.ka_end_lat(2), chd.ka_end_lat_idx(2)]     = min(abs(1.45076*180/pi - GR.lat_ka));

    [chd.ka_start_lon(2), chd.ka_start_lon_idx(2)] = min(abs(-0.39278*180/pi - GR.lon_ka));
    [chd.ka_end_lon(2), chd.ka_end_lon_idx(2)]     = min(abs(-0.39304*180/pi - GR.lon_ka));

    % ZONE 3   
    [chd.ka_start_lat(3), chd.ka_start_lat_idx(3)] = min(abs(1.45076*180/pi - GR.lat_ka));
    [chd.ka_end_lat(3), chd.ka_end_lat_idx(3)]     = min(abs(1.456715*180/pi - GR.lat_ka));

    [chd.ka_start_lon(3), chd.ka_start_lon_idx(3)] = min(abs(-0.39304*180/pi - GR.lon_ka));
    [chd.ka_end_lon(3), chd.ka_end_lon_idx(3)]     = min(abs(-0.40856*180/pi - GR.lon_ka));
    
    %% Expected values
    chd.land_ice_elevation = 2.05 + 2.168e-1; % = Reference_Height(Rough Surface) + Thickness(Ice)

    %% Requirements 
    chd.elevation_accuracy = 2;