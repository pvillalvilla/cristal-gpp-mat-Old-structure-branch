function L1B_ML = reading_L1B_ML(filesBulk)
L1B_ML.lat_ku            = double(ncread(filesBulk.filename_L1B_ML, 'lat_ku'));
L1B_ML.lat_ka            = double(ncread(filesBulk.filename_L1B_ML, 'lat_ka'));
L1B_ML.lon_ku            = double(ncread(filesBulk.filename_L1B_ML, 'lon_ku'));
L1B_ML.lon_ka            = double(ncread(filesBulk.filename_L1B_ML, 'lon_ka'));
L1B_ML.altitude_ku       = double(ncread(filesBulk.filename_L1B_ML, 'altitude_ku'));
L1B_ML.altitude_ka       = double(ncread(filesBulk.filename_L1B_ML, 'altitude_ka'));
L1B_ML.range_ku          = double(ncread(filesBulk.filename_L1B_ML, 'range_ku'));
L1B_ML.range_ka          = double(ncread(filesBulk.filename_L1B_ML, 'range_ka'));
end