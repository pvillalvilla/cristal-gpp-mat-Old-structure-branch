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

function write_NetCDF_L2(files,L1A_buffer,L1BS, L1B,i_surf_stacked,cnf,chd,cst)
%{
i_surf = i_surf_stacked-1;
global cst.c chd.bw_ku_chd  cnf.zp_fact_range cst.sec_in_day cst.pi
global mission cnf.mode chd.N_samples_sar N_max_beams_stack_chd cnf.processing_mode
global cnf.compute_L1_stadistics
global cnf.include_wfms_aligned
%added EM: 04.10.2016
global ACDC_application_cnf

%}
%global netcdf_type
% ncid = netcdf.open(files.filename_netCDF,'WRITE');

% ku_rec_dimension = netcdf.defDim(ncid,'time_l1b_echo',N_surfs_loc_estimated);
% nl_dimension = netcdf.defDim(ncid,'Nl',N_max_beams_stack_chd);
% ns_dimension = netcdf.defDim(ncid,'Ns',chd.N_samples_sar*cnf.zp_fact_range);
% space_3D_dimension = netcdf.defDim(ncid,'space_3D',3);

% t6 = tic;

% date_creation = datestr(now, '_yyyymmdd_');
% switch mission
%     case 'CR2'        
%         switch cnf.mode
%             case 'SAR'                
%                 files.filename_netCDF = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA____',...
%                                 files.sph.product_info.product_id(20:20+30),...
%                                 date_creation,...
%                                 'isd','.nc');
%         end
%     case {'S3_','S3A','S3B'} 
%     case {'S6_'} 
% end


% dimensions_key = 'Dimensions';
% format_key = 'Format';
% data_type_key = 'DataType';
% fill_value_key='FillValue';
% 
% netcdf_v4_format = 'netcdf4';
% 
% ku_rec_dimension = 'time_l1b_echo');
% nl_dimension = 'Nl';
% nl_dimension_size = N_max_beams_stack_chd;
% ns_dimension = 'Ns';
% space_3D_dimension = 'space_3D';
% space_3D_dimension_size = 3;



% int8_type = 'int8';
% uint8_type = 'uint8';
% int16_type = 'int16';
% uint16_type = 'uint16';
% int32_type = 'int32';
% uint32_type = 'uint32';
% uint64_type = 'uint64';
% float_type = 'single';
% double_type='double';


time_utc    = L2.ISD_time_surf(L2.idx_int_ISD); %leap seconds are already taken into account in L1B time values
UTC_day_L2_echo = L2.ISD_UTC_day(L2.idx_int_ISD);
UTC_sec_L2_echo = L2.ISD_UTC_sec(L2.idx_int_ISD);
%L2_record_counter = ;

% UTC_day_L2_echo = int16(floor(L1BS.time_surf./ sec_in_day_cst));
% UTC_sec_L2_echo = double((time_utc - double(UTC_day_L2_echo) * sec_in_day_cst));


%---------- Orbit variables ------------------------------
latitude_sat = int32(L2.ISD_lat_sat(L2.idx_int_ISD)         * 1e6);
Longitude_sat = int32(L2.ISD_lon_sat(L2.idx_int_ISD)         * 1e6); 
com_altitude_sat = int32((L2.ISD_alt_sat(L2.idx_int_ISD)-700000) * 1e4);

com_altitude_rate = int32(L2.ISD_alt_rate(L2.idx_int_ISD)); 

com_position_vector = double(L2.ISD_pos_vector(L2.idx_int_ISD)); 
com_position_velocity = double(L2.ISD_pos_velocity(L2.idx_int_ISD)); 

off_nadir_roll_angle_pf = double(L2.ISD_roll_angle(L2.idx_int_ISD)); 
off_nadir_pitch_angle_pf = double(L2.ISD_pitch_angle(L2.idx_int_ISD)); 
off_nadir_yaw_angle_pf = double(L2.ISD_yaw_angle(L2.idx_int_ISD)); 

%---------- Geophysical corrections --------------------------------------
wet_tropo = double(L2.wet_tropo(L2.idx_int_ISD)); 
geoid = double(L2.geoid(L2.idx_int_ISD)); 
mss = double(L2.mss(L2.idx_int_ISD)); 

%---------- Measurement Parameters --------------------------------------
latitude_surf = int32(L2.ISD_lat_surf(L2.idx_int_ISD)         * 1e6);
Longitude_surf = int32(L2.ISD_lon_surf(L2.idx_int_ISD)         * 1e6); 

%---------- Empirical retrievals --------------------------------------
height_ku_xxxx = double(L2.height_ku(L2.idx_int_ISD));
height_ka_xxxx = double(L2.height_ka(L2.idx_int_ISD));
backs_ku_xxxx = int32(L2.backs_ku(L2.idx_int_ISD));
backs_ka_xxxx = int32(L2.backs_ka(L2.idx_int_ISD));

%---------- Physical retrievals --------------------------------------
ssh_ku = double(L2.ssh_ku(L2.idx_int_ISD));
swh_ku = double(L2.swh_ku(L2.idx_int_ISD));
wind_ku = double(L2.wind_ku(L2.idx_int_ISD));
ssh_ka = double(L2.ssh_ka(L2.idx_int_ISD));
swh_ka = double(L2.swh_ka(L2.idx_int_ISD));
wind_ka = double(L2.wind_ka(L2.idx_int_ISD));

%---------- Sea ice retrievals --------------------------------------
freeboard = double(L2.freeboard(L2.idx_int_ISD));
snow_depth = double(L2.snow_depth(L2.idx_int_ISD));
ice_thickness = double(L2.ice_thickness(L2.idx_int_ISD));

%----------  Swath Processing Retrievals --------------------------------
bin_index = int32(L2.bin_index(L2.idx_int_ISD));
height_swath = int32(L2.height_swath(L2.idx_int_ISD));
Latitude_swath = int32(L2.Latitude_swath(L2.idx_int_ISD)         * 1e6);
Longitude_swath = int32(L2.Longitude_swath(L2.idx_int_ISD)         * 1e6);
Unwrapped_AoA = double(L2.Unwrapped_AoA(L2.idx_int_ISD));
coherence = double(L2.coherence(L2.idx_int_ISD));
dem_diff = double(L2.dem_diff(L2.idx_int_ISD));
dem_diff_mad = double(L2.dem_diff_mad(L2.idx_int_ISD));
%{
%---------- Altimeter range -------------------------
range_ice_sheet_L2_echo                 = int32((L2.ISD_range-700000)   * 1e4);
range_ice_sheet_L2_echo_mean_error      = int32(L2.range_mean_error     * 1e4);
range_ice_sheet_L2_echo_rmse            = int32(L2.range_rmse           * 1e4);

%---------- Sigma0 --------------------------------
sig0_ice_sheet_L2_echo                  = int16(L2.ISD_sigma0           * 1e2);
sig0_ice_sheet_L2_echo_mean_error       = int16(L2.sig0_mean_error      * 1e2);
sig0_ice_sheet_L2_echo_rmse             = int16(L2.sig0_rmse            * 1e2);

%---------- Elevation of echoing points -------------------------
elevation_ice_sheet_L2_echo             = int32((L2.ISD_SSH)            * 1e4);
elevation_ice_sheet_L2_echo_mean_error  = int32(L2.elev_mean_error      * 1e4);
elevation_ice_sheet_L2_echo_rmse        = int32(L2.elev_rmse            * 1e4);
%}


%% PACKING L2
%----------A. Time variables ----------------------------------------------
var_id = netcdf.inqVarID(L2.ncid,'time_utc');
netcdf.putVar(L2.ncid,var_id,time_utc);

var_id = netcdf.inqVarID(L2.ncid,'UTC_day_L2_echo');
netcdf.putVar(L2.ncid,var_id,UTC_day_L2_echo);

var_id = netcdf.inqVarID(L2.ncid,'UTC_sec_L2_echo');
netcdf.putVar(L2.ncid,var_id,UTC_sec_L2_echo);
%{
var_id = netcdf.inqVarID(L2.ncid,'L2_record_counter');
netcdf.putVar(L2.ncid,var_id,L2_record_counter);
%}

%------------------------ Orbit variables ---------------------------------
var_id = netcdf.inqVarID(L2.ncid,'latitude_sat');
netcdf.putVar(L2.ncid,var_id,latitude_sat);

var_id = netcdf.inqVarID(L2.ncid,'Longitude_sat');
netcdf.putVar(L2.ncid,var_id,Longitude_sat);

var_id = netcdf.inqVarID(L2.ncid,'com_altitude_sat');
netcdf.putVar(L2.ncid,var_id,com_altitude_sat);


%------------------------ Altimeter range ---------------------------------
var_id = netcdf.inqVarID(L2.ncid,'range_ice_sheet_L2_echo');
netcdf.putVar(L2.ncid,var_id,range_ice_sheet_L2_echo);

var_id = netcdf.inqVarID(L2.ncid,'range_ice_sheet_L2_echo_mean_error');
netcdf.putVar(L2.ncid,var_id,range_ice_sheet_L2_echo_mean_error);

var_id = netcdf.inqVarID(L2.ncid,'range_ice_sheet_L2_echo_rmse');
netcdf.putVar(L2.ncid,var_id,range_ice_sheet_L2_echo_rmse);


%--------------------- Sigma0 scaling factor ------------------------------
var_id = netcdf.inqVarID(L2.ncid,'sig0_ice_sheet_L2_echo');
netcdf.putVar(L2.ncid,var_id,sig0_ice_sheet_L2_echo);

var_id = netcdf.inqVarID(L2.ncid,'sig0_ice_sheet_L2_echo_mean_error');
netcdf.putVar(L2.ncid,var_id,sig0_ice_sheet_L2_echo_mean_error);

var_id = netcdf.inqVarID(L2.ncid,'sig0_ice_sheet_L2_echo_rmse');
netcdf.putVar(L2.ncid,var_id,sig0_ice_sheet_L2_echo_rmse);


%----------------  Elevation of echoing points ----------------------------
var_id = netcdf.inqVarID(L2.ncid,'elevation_ice_sheet_L2_echo');
netcdf.putVar(L2.ncid,var_id,elevation_ice_sheet_L2_echo);

var_id = netcdf.inqVarID(L2.ncid,'elevation_ice_sheet_L2_echo_mean_error');
netcdf.putVar(L2.ncid,var_id,elevation_ice_sheet_L2_echo_mean_error);

var_id = netcdf.inqVarID(L2.ncid,'elevation_ice_sheet_L2_echo_rmse');
netcdf.putVar(L2.ncid,var_id,elevation_ice_sheet_L2_echo_rmse);



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