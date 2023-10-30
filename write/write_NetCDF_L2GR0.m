%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the CODING & PACKING 
% algorithm for Level-1B product using Sentinel-3 like format
% Ref: Product Data Format Specification - SRAL/MWR Level 1 & 2 Instrument
% Products issue 2.0
%
% ---------------------------------------------------------
% Objective: Pack variables and write the NETCDF
% 
% INPUTs : Workspace
% OUTPUTs: TM Structure as defined on isardSAT_JasonCS_DPM
%
% ----------------------------------------------------------
% Author:    Eduard Makhoul/ isardSAT
%            Gorka Moyano  / isardSAT
%            Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT ()
% 
% Versions
% 1.0 
% 1.1 Updated time conversion for data between 2010 and 2016 (CR2)
% 2.0 Transformed to a function. Writting one record
% 2.1 changed chd.N_samples_sar_sar_chd by chd.N_samples_sar
% 2.2 (04.10.2016, EM) Adding ACDC output related variables
% 2.3 nb_stack_l1b_echo is L1B.N_beams_contributing instead of L1BS.N_beams_stack
% 2.4 added cnf flag cnf.processing_mode
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils

function write_NetCDF_L2GR0(files,L2,N_rec,cnf,chd,cst, ncid)


% time_utc    = L2.ISD_time_surf(L2.idx_int_ISD); %leap seconds are already taken into account in L1B time values
% UTC_day_L2_echo = L2.ISD_UTC_day(L2.idx_int_ISD);
% UTC_sec_L2_echo = L2.ISD_UTC_sec(L2.idx_int_ISD);
%L2_record_counter = ;

% UTC_day_L2_echo = int16(floor(L1BS.time_surf./ sec_in_day_cst));
% UTC_sec_L2_echo = double((time_utc - double(UTC_day_L2_echo) * sec_in_day_cst));

%---------- Time and counter ------------------------------
time = double(L2.time);
% time_tai = double(L2.time_tai);
% l2_record_counter = int32(L2.l2_record_counter);
% l1b_record_counter = int32(L2.l1b_record_counter);
% tm_source_sequence_counter = int32(L2.tm_source_sequence_counter);

%---------- Orbit variables ------------------------------
altitude = int32((L2.alt-700000) * 1e4);
altitude_rate = int32(L2.alt_rate * 1e4); 
latitude = int32(L2.lat * 1e6);
longitude = int32(L2.lon * 1e6); 
% latitude_surf = int32(L2.lat_surf * 1e6);
% longitude_surf = int32(L2.lon_surf * 1e6);

position_vector = [double(L2.x_pos),double(L2.y_pos),double(L2.z_pos)];
velocity_vector = [int32(L2.x_vel * 1e4),double(L2.y_vel * 1e4),double(L2.z_vel * 1e4)];

off_nadir_pitch_angle = int32(L2.off_nadir_pitch_angle * 1e7); %degrees?
off_nadir_roll_angle = int32(L2.off_nadir_roll_angle * 1e7); %degrees? 
off_nadir_yaw_angle = int32(L2.off_nadir_yaw_angle * 1e7); %degrees?

%---------- Configuration variables ------------------------------
% iris_mode_flag = int8(L2.iris_mode_flag);
% telemetry_type_flag = int8(L2.telemetry_type_flag);

%---------- Altimeter range variables ------------------------------
tracker_range_calibrated = int32((L2.tracker_range_calibrated - 700000) * 1e4);

%---------- Altimeter power variables ------------------------------
sig0_scaling_factor = int16(L2.sig0_scaling_factor * 1e2);

%---------- Altimeter flag variables ------------------------------
surface_classification_flag = int8(L2.surface_classification_flag);

%---------- Geophysical corrections variables --------------------------------------
rad_wet_tropo_cor = int16(L2.rad_wet_tropo_cor * 1e4);
rad_atm_cor_sig0 = int16(L2.rad_atm_cor_sig0 * 1e2);
ocean_geo_corrections = int32(L2.ocean_geo_corrections * 1e4);
seaice_geo_corrections = int32(L2.seaice_geo_corrections * 1e4);
lead_geo_corrections = int32(L2.lead_geo_corrections * 1e4); 
iceshelve_geo_corrections = int32(L2.iceshelve_geo_corrections * 1e4); 
landice_geo_corrections = int32(L2.landice_geo_corrections * 1e4); 
inland_geo_corrections = int32(L2.inland_geo_corrections * 1e4);

%---------- Reference surfaces variables --------------------------------------
% geoid = int32(L2.geoid * 1e4);
mean_sea_surface = int32(L2.mean_sea_surface * 1e4);

%---------- Waveform characteristics variables --------------------------------------
peakiness = int16(L2.peakiness * 1e2);
noise_floor = int32(L2.noise_floor * 1e6);
% number_of_peaks = int8(L2.number_of_peaks);
% prominence = int16(L2.prominence * 1e2);
% slope_trailing_edge = int32(L2.slope_trailing_edge);
peak_noise_ratio = int32(L2.peak_noise_ratio * 1e6);

%---------- Retracking estimates variables --------------------------------------
epoch_ocog = int32(L2.epoch_ocog * 1e15);
epoch_tcog = int32(L2.epoch_tcog * 1e15);
epoch_tfmra = int32(L2.epoch_tfmra * 1e15);
epoch = int32(L2.epoch * 1e15);

amplitude_ocog = int32(L2.amplitude_ocog * 1e3);
amplitude_tcog = int32(L2.amplitude_tcog * 1e3);
amplitude_tfmra = int32(L2.amplitude_tfmra * 1e3);
amplitude = int32(L2.amplitude * 1e3);

%---------- Retracking retrievals variables --------------------------------------
range_ocog = int32((L2.range_ocog - 700000)* 1e4);
range_tcog = int32((L2.range_tcog - 700000)* 1e4);
range_tfmra = int32((L2.range_tfmra - 700000)* 1e4);
range = int32((L2.range - 700000)* 1e4);

sig0_ocog = int16(L2.sig0_ocog * 1e2);
sig0_tcog = int16(L2.sig0_tcog * 1e2);
sig0_tfmra = int16(L2.sig0_tfmra * 1e2);
sig0 = int16(L2.sig0 * 1e2);

%---------- Ocean geophysical retrievals variables --------------------------------------
ssha = int32(L2.ssha * 1e4);
swh = int16(L2.swh * 1e3);
% wind_speed_alt = int16(L2.wind_speed_alt * 1e3);

%---------- Snow geophysical retrievals variables --------------------------------------
snow_depth = int16 (L2.snow_depth * 1e4);
snow_depth_correction = int16 (L2.snow_depth_correction * 1e4);
snow_depth_correction = int16 (L2.snow_depth_correction * 1e4);

%---------- Sea-ice geophysical retrievals variables --------------------------------------
sea_ice_concentration = int16(L2.sea_ice_concentration * 1e2);
sea_ice_type = int8(L2.sea_ice_type);
waveform_classification_flag = int8(L2.waveform_classification_flag);
sea_ice_freeboard = int32(L2.sea_ice_freeboard * 1e4);
sea_ice_thickness = int32(L2.sea_ice_thickness * 1e4);
ssha_lead_interpolated = int32(L2.ssha_lead_interpolated * 1e4);
% iceberg_detection_flag = int8(L2.iceberg_detection_flag);
% iceberg_freeboard = int32(L2.iceberg_freeboard * 1e4);
% iceberg_sig0 = int16(L2.iceberg_sig0 * 1e2);

%---------- Land Ice geophysical retrievals variables --------------------------------------
land_ice_elevation = int32(L2.land_ice_elevation * 1e4);

%---------- Inland waters retrievals variables--------------------------------------
water_level_height = int32(L2.water_level_height * 1e4);


%% PACKING L2
ncid = netcdf.open(files.filename_L2GR);
% data_ncid = netcdf.inqNcid(ncid,'data');
% data_ku_ncid = netcdf.inqNcid(data_ncid,'ku');
% data_ka_ncid = netcdf.inqNcid(data_ncid,'ku');

%---------- Time and counter ------------------------------
var_id = netcdf.inqVarID(ncid,'time');
netcdf.putVar(ncid,var_id,time);

% var_id = netcdf.inqVarID(ncid,'l2_record_counter');
% netcdf.putVar(ncid,var_id,l2_record_counter);
% 
% var_id = netcdf.inqVarID(ncid,'l1b_record_counter');
% netcdf.putVar(ncid,var_id,l1b_record_counter);
% 
% var_id = netcdf.inqVarID(data_ku_ncid,'tm_source_sequence_counter');
% netcdf.putVar(data_ku_ncid,var_id,tm_source_sequence_counter);
%------------------------ Orbit variables ---------------------------------
var_id = netcdf.inqVarID(ncid,'altitude');
netcdf.putVar(ncid,var_id,altitude);

var_id = netcdf.inqVarID(ncid,'altitude_rate');
netcdf.putVar(ncid,var_id,altitude_rate);

var_id = netcdf.inqVarID(ncid,'latitude');
netcdf.putVar(ncid,var_id,latitude);

var_id = netcdf.inqVarID(ncid,'longitude');
netcdf.putVar(ncid,var_id,longitude);

var_id = netcdf.inqVarID(ncid,'latitude_surf');
netcdf.putVar(ncid,var_id,latitude_surf);

var_id = netcdf.inqVarID(ncid,'longitude_surf');
netcdf.putVar(ncid,var_id,longitude_surf);

var_id = netcdf.inqVarID(ncid,'position_vector');
netcdf.putVar(ncid,var_id,position_vector);

var_id = netcdf.inqVarID(ncid,'velocity_vector');
netcdf.putVar(ncid,var_id,velocity_vector);

var_id = netcdf.inqVarID(ncid,'off_nadir_pitch_angle');
netcdf.putVar(ncid,var_id,off_nadir_pitch_angle);

var_id = netcdf.inqVarID(ncid,'off_nadir_roll_angle');
netcdf.putVar(ncid,var_id,off_nadir_roll_angle);

var_id = netcdf.inqVarID(ncid,'off_nadir_yaw_angle');
netcdf.putVar(ncid,var_id,off_nadir_roll_angle);
%--------- Configuration Variables --------------------------------------
var_id = netcdf.inqVarID(ncid,'iris_mode_flag');
netcdf.putVar(ncid,var_id,iris_mode_flag);

var_id = netcdf.inqVarID(ncid,'telemetry_type_flag');
netcdf.putVar(ncid,var_id,telemetry_type_flag);

%---------- Altimeter range variables --------------------------------------
var_id = netcdf.inqVarID(ncid,'tracker_range_calibrated');
netcdf.putVar(ncid,tracker_range_calibrated);
%---------- Altimeter power variables --------------------------------------

var_id = netcdf.inqVarID(ncid,'sig0_scaling_factor');
netcdf.putVar(ncid,var_id,tracker_range_calibrated);

%---------- Flag Variables --------------------------------------

var_id = netcdf.inqVarID(ncid,'surface_classification_flag');
netcdf.putVar(ncid,var_id,surface_classification_flag);

%---------- Geophysical corrections variables --------------------------------------

var_id = netcdf.inqVarID(ncid,'rad_wet_tropo_cor');
netcdf.putVar(ncid,var_id,rad_wet_tropo_cor);

var_id = netcdf.inqVarID(ncid,'rad_atm_cor_sig0');
netcdf.putVar(ncid,var_id,rad_atm_cor_sig0);

var_id = netcdf.inqVarID(ncid,'ocean_geo_corrections');
netcdf.putVar(ncid,var_id,ocean_geo_corrections);

var_id = netcdf.inqVarID(ncid,'seaice_geo_corrections');
netcdf.putVar(ncid,var_id,seaice_geo_corrections);

var_id = netcdf.inqVarID(ncid,'lead_geo_corrections');
netcdf.putVar(ncid,var_id,lead_geo_corrections);

var_id = netcdf.inqVarID(ncid,'iceshelve_geo_corrections');
netcdf.putVar(ncid,var_id,iceshelve_geo_corrections);

var_id = netcdf.inqVarID(ncid,'landice_geo_corrections');
netcdf.putVar(ncid,var_id,landice_geo_corrections);

var_id = netcdf.inqVarID(ncid,'inland_geo_corrections');
netcdf.putVar(ncid,var_id,inland_geo_corrections);

%---------- Reference surfaces variables --------------------------------------
% var_id = netcdf.inqVarID(data_ncid,'geoid');
% netcdf.putVar(data_ncid,var_id,geoid);
% 
% var_id = netcdf.inqVarID(data_ncid,'mean_sea_surface');
% netcdf.putVar(data_ncid,var_id,mean_sea_surface);

%---------- Waveform characteristics variables --------------------------------------

var_id = netcdf.inqVarID(ncid,'peakiness');
netcdf.putVar(ncid,var_id,peakiness);

var_id = netcdf.inqVarID(ncid,'noise_floor');
netcdf.putVar(ncid,var_id,noise_floor);

var_id = netcdf.inqVarID(ncid,'number_of_peaks');
netcdf.putVar(ncid,var_id,number_of_peaks);

var_id = netcdf.inqVarID(ncid,'prominence');
netcdf.putVar(ncid,var_id,prominence);

var_id = netcdf.inqVarID(ncid,'slope_trailing_edge');
netcdf.putVar(ncid,var_id,slope_trailing_edge);

var_id = netcdf.inqVarID(ncid,'peak_noise_ratio');
netcdf.putVar(ncid,var_id,peak_noise_ratio);

%---------- Retracking estimates variables --------------------------------------
var_id = netcdf.inqVarID(ncid,'epoch_ocog');
netcdf.putVar(ncid,var_id,epoch_ocog);

var_id = netcdf.inqVarID(ncid,'epoch_tcog');
netcdf.putVar(ncid,var_id,epoch_tcog);

var_id = netcdf.inqVarID(ncid,'epoch_tfmra');
netcdf.putVar(ncid,var_id,epoch_tfmra);

var_id = netcdf.inqVarID(ncid,'epoch');
netcdf.putVar(ncid,var_id,epoch);

var_id = netcdf.inqVarID(ncid,'amplitude_ocog');
netcdf.putVar(ncid,var_id,amplitude_ocog);

var_id = netcdf.inqVarID(ncid,'amplitude_tcog');
netcdf.putVar(ncid,var_id,amplitude_tcog);

var_id = netcdf.inqVarID(ncid,'amplitude_tfmra');
netcdf.putVar(ncid,var_id,amplitude_tfmra);

var_id = netcdf.inqVarID(ncid,'amplitude');
netcdf.putVar(ncid,var_id,amplitude);

%---------- Retracking retrievals variables --------------------------------------
var_id = netcdf.inqVarID(ncid,'range_ocog');
netcdf.putVar(ncid,var_id,range_ocog);

var_id = netcdf.inqVarID(ncid,'range_tcog');
netcdf.putVar(ncid,var_id,range_tcog);

var_id = netcdf.inqVarID(ncid,'range_tfmra');
netcdf.putVar(ncid,var_id,range_tfmra);

var_id = netcdf.inqVarID(ncid,'range');
netcdf.putVar(ncid,var_id,range);

var_id = netcdf.inqVarID(ncid,'sig0_ocog');
netcdf.putVar(ncid,var_id,sig0_ocog);

var_id = netcdf.inqVarID(ncid,'sig0_tcog');
netcdf.putVar(ncid,var_id,sig0_tcog);

var_id = netcdf.inqVarID(ncid,'sig0_tfmra');
netcdf.putVar(ncid,var_id,sig0_tfmra);

var_id = netcdf.inqVarID(ncid,'sig0');
netcdf.putVar(ncid,var_id,sig0);

%---------- Ocean geophysical retrievals variables --------------------------------------

var_id = netcdf.inqVarID(ncid,'ssha');
netcdf.putVar(ncid,var_id,ssha);

var_id = netcdf.inqVarID(ncid,'swh');
netcdf.putVar(ncid,var_id,swh);

var_id = netcdf.inqVarID(ncid,'wind_speed_alt');
netcdf.putVar(ncid,var_id,wind_speed_alt);

%---------- Snow geophysical retrievals variables --------------------------------------
% var_id = netcdf.inqVarID(data_ncid,'snow_depth');
% netcdf.putVar(data_ncid,var_id,snow_depth);

var_id = netcdf.inqVarID(ncid,'snow_depth_correction');
netcdf.putVar(ncid,var_id,snow_depth_correction);

var_id = netcdf.inqVarID(ncid,'snow_interface_alignement');
netcdf.putVar(ncid,var_id,snow_interface_alignement);

%---------- Sea-ice geophysical retrievals variables --------------------------------------

var_id = netcdf.inqVarID(ncid,'sea_ice_concentration');
netcdf.putVar(ncid,var_id,sea_ice_concentration);

var_id = netcdf.inqVarID(ncid,'sea_ice_type');
netcdf.putVar(ncid,var_id,sea_ice_type);

var_id = netcdf.inqVarID(ncid,'waveform_classification_flag');
netcdf.putVar(ncid,var_id,waveform_classification_flag);

var_id = netcdf.inqVarID(ncid,'sea_ice_freeboard');
netcdf.putVar(ncid,var_id,sea_ice_freeboard);

var_id = netcdf.inqVarID(ncid,'sea_ice_thickness');
netcdf.putVar(ncid,var_id,sea_ice_thickness);

var_id = netcdf.inqVarID(ncid,'ssha_lead_interpolated');
netcdf.putVar(ncid,var_id,ssha_lead_interpolated);

var_id = netcdf.inqVarID(ncid,'iceberg_detection_flag');
netcdf.putVar(ncid,var_id,iceberg_detection_flag);

var_id = netcdf.inqVarID(ncid,'iceberg_freeboard');
netcdf.putVar(ncid,var_id,iceberg_freeboard);

var_id = netcdf.inqVarID(ncid,'iceberg_sig0');
netcdf.putVar(ncid,var_id,iceberg_sig0);

%---------- Land Ice geophysical retrievals variables --------------------------------------

var_id = netcdf.inqVarID(ncid,'land_ice_elevation');
netcdf.putVar(ncid,var_id,land_ice_elevation);

%---------- Inland waters retrievals variables --------------------------------------

var_id = netcdf.inqVarID(ncid,'water_level_height');
netcdf.putVar(ncid,var_id,water_level_height);

netcdf.close(L2.ncid);

end