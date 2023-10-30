function GR = reading_GR(filesBulk)
GR.lat_ku     = double(ncread(filesBulk.filename_L2_GR, 'lat_ku'));
GR.lat_ka     = double(ncread(filesBulk.filename_L2_GR, 'lat_ka'));
GR.lon_ku     = double(ncread(filesBulk.filename_L2_GR, 'lon_ku'));
GR.lon_ka     = double(ncread(filesBulk.filename_L2_GR, 'lon_ka'));
GR.land_ice_elevation_ku     = double(ncread(filesBulk.filename_L2_GR, '.land_ice_elevation_ku'));
end

