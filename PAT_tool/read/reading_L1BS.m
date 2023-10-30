function L1BS = reading_L1BS(filesBulk)
L1BS.lat_ku            = double(ncread(filesBulk.filename_L1BS, 'lat_ku'));
L1BS.lat_ka            = double(ncread(filesBulk.filename_L1BS, 'lat_ka'));
L1BS.lon_ku            = double(ncread(filesBulk.filename_L1BS, 'lon_ku'));
L1BS.lon_ka            = double(ncread(filesBulk.filename_L1BS, 'lon_ka'));
L1BS.time_ku            = double(ncread(filesBulk.filename_L1BS, 'time_ku'));
L1BS.time_ka            = double(ncread(filesBulk.filename_L1BS, 'time_ka'));
L1BS.range_ku            = double(ncread(filesBulk.filename_L1BS, 'range_ku'));
L1BS.tracker_range_calibrated_ku            = double(ncread(filesBulk.filename_L1BS, 'tracker_range_calibrated_ku'));
L1BS.tracker_range_calibrated_ka            = double(ncread(filesBulk.filename_L1BS, 'tracker_range_calibrated_ka'));
L1BS.slant_range_correction_applied_ku      = double(ncread(filesBulk.filename_L1BS, 'slant_range_correction_applied_ku'));
L1BS.slant_range_correction_applied_ka      = double(ncread(filesBulk.filename_L1BS, 'slant_range_correction_applied_ka'));
L1BS.power_waveform_ku      = double(ncread(filesBulk.filename_L1BS, 'power_waveform_ku'));
L1BS.power_waveform_ka      = double(ncread(filesBulk.filename_L1BS, 'power_waveform_ka'));
L1BS.Gap_flag_ku            = double(ncread(filesBulk.filename_L1BS, 'Gap_flag_ku'));
L1BS.Gap_flag_ka            = double(ncread(filesBulk.filename_L1BS, 'Gap_flag_ka'));
end