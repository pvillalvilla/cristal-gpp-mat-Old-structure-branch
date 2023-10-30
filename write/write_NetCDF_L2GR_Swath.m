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

function write_NetCDF_L2GR(files,SWATH,N_rec, max_swath_length,cnf,chd,cst)


% time_utc    = L2.ISD_time_surf(L2.idx_int_ISD); %leap seconds are already taken into account in L1B time values
% UTC_day_L2_echo = L2.ISD_UTC_day(L2.idx_int_ISD);
% UTC_sec_L2_echo = L2.ISD_UTC_sec(L2.idx_int_ISD);
%L2_record_counter = ;

% UTC_day_L2_echo = int16(floor(L1BS.time_surf./ sec_in_day_cst));
% UTC_sec_L2_echo = double((time_utc - double(UTC_day_L2_echo) * sec_in_day_cst));

%---------- Time ------------------------------
time = double(SWATH.time);
time_tai = double(SWATH.time_tai);
%---------- Swath retrievals variables ------------------------------
length_swath = int8(length(SWATH.lat));
bin_index = int16(SWATH.sample);
height_swath = int32((SWATH.elev - 700000)* 1e4);
latitude_swath = int32(SWATH.lat* 1e6);
longitude_swath = int32(SWATH.lon* 1e6);
angle_of_arrival_unwrapped = int32(SWATH.AoA* 1e6);
coherence = int32(SWATH.coherence* 1e3);
dem_diff = int32(SWATH.individual_dem_diff* 1e4);
dem_diff_mad = int32(SWATH.MAD* 1e4);

%---------- Land Ice geophysical retrievals variables ------------------------------
land_ice_elevation = int32(SWATH.land_ice_elevation* 1e4);

%---------- Inland waters retrievals variables ------------------------------
water_level_height = int32(SWATH.water_level_height* 1e4);


%% PACKING L2

data_ncid = netcdf.inqNcid(ncid,'data');
data_ku_ncid = netcdf.inqNcid(data_ncid,'ku');

%---------- Time  ------------------------------
var_id = netcdf.inqVarID(data_ku_ncid,'time');
netcdf.putVar(data_ku_ncid,var_id,time);

var_id = netcdf.inqVarID(data_ku_ncid,'time_tai');
netcdf.putVar(data_ku_ncid,var_id,time_tai);

%---------- Swath retrievals variables ------------------------------
var_id = netcdf.inqVarID(data_ku_ncid,'length_swath');
netcdf.putVar(data_ku_ncid,var_id,length_swath);

var_id = netcdf.inqVarID(data_ku_ncid,'bin_index');
netcdf.putVar(data_ku_ncid,var_id,bin_index);

var_id = netcdf.inqVarID(data_ku_ncid,'height_swath');
netcdf.putVar(data_ku_ncid,var_id,height_swath);

var_id = netcdf.inqVarID(data_ku_ncid,'latitude_swath');
netcdf.putVar(data_ku_ncid,var_id,latitude_swath);

var_id = netcdf.inqVarID(data_ku_ncid,'longitude_swath');
netcdf.putVar(data_ku_ncid,var_id,longitude_swath);

var_id = netcdf.inqVarID(data_ku_ncid,'angle_of_arrival_unwrapped');
netcdf.putVar(data_ku_ncid,var_id,angle_of_arrival_unwrapped);

var_id = netcdf.inqVarID(data_ku_ncid,'coherence');
netcdf.putVar(data_ku_ncid,var_id,coherence);

var_id = netcdf.inqVarID(data_ku_ncid,'dem_diff');
netcdf.putVar(data_ku_ncid,var_id,dem_diff);

var_id = netcdf.inqVarID(data_ku_ncid,'dem_diff_mad');
netcdf.putVar(data_ku_ncid,var_id,dem_diff_mad);

%---------- Land Ice geophysical retrievals variables --------------------------------------
var_id = netcdf.inqVarID(data_ku_ncid,'land_ice_elevation');
netcdf.putVar(data_ku_ncid,var_id,land_ice_elevation);



%---------- Inland waters retrievals variables --------------------------------------
var_id = netcdf.inqVarID(data_ku_ncid,'water_level_height');
netcdf.putVar(data_ku_ncid,var_id,water_level_height);



% %----------  Global Attributes definition -----------------------------------
% %---- attributes inherited from Sentinel-3 product description-------------
% 
% ncwriteatt(L2.filename_netCDF,'/','Conventions',netcdf_v4_format);
% ncwriteatt(L2.filename_netCDF,'/','altimeter_sensor_name',altimeter_sensor_name);
% ncwriteatt(L2.filename_netCDF,'/','first_meas_time',first_meas_time);
% ncwriteatt(L2.filename_netCDF,'/','last_meas_time',last_meas_time);
% ncwriteatt(L2.filename_netCDF,'/','first_meas_lat',first_meas_lat);
% ncwriteatt(L2.filename_netCDF,'/','last_meas_lat',last_meas_lat);
% ncwriteatt(L2.filename_netCDF,'/','first_meas_lon',first_meas_lon);
% ncwriteatt(L2.filename_netCDF,'/','last_meas_lon',last_meas_lon);
% ncwriteatt(L2.filename_netCDF,'/','semi_major_ellipsoid_axis',semi_major_ellipsoid_axis);
% ncwriteatt(L2.filename_netCDF,'/','ellipsoid_flattening',ellipsoid_flattening);
% %--------------- add the attributes related to intermediate product--------
% ncwriteatt(L2.filename_netCDF,'/','orbit_cycle_num',orbit_cycle_num);
% ncwriteatt(L2.filename_netCDF,'/','orbit_REL_Orbit',orbit_REL_Orbit);
% ncwriteatt(L2.filename_netCDF,'/','orbit_ABS_Orbit_Start',orbit_ABS_Orbit_Start);


netcdf.close(L2.ncid);

end