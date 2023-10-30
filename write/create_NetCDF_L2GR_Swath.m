%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the CODING & PACKING 
% algorithm for Level-2 product using Sentinel-3 like format
% Ref: S3-TN-ESA-SR-0433 SRAL L1A-1BS IODD V1.4
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
% Last rev.: Monica Roca   / isardSAT 
% 
% Versions
% 1.0 
% 1.1 Updated time conversion for data between 2010 and 2016 (CR2)
% 2.0 Transformed to a function. Writting one record
% 2.1 zp changed to int16 as int8 did not work for zp > 128
% 2.2 case S3_ filename
% 2.3 case 'SIN' added 
% 2.4 added cnf flag processing_mode_cnf
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function [files] = create_NetCDF_L2GR(files, N_rec, max_swath_length,cnf, chd, cst)

%global N_samples   semi_major_axis_cst flat_coeff_cst;
%global bw_ku_chd 
%zp_fact_range_cnf  sec_in_day_cst pi_cst
%global mission mode N_max_beams_stack_chd
%global compute_L1_stadistics_cnf include_wfms_aligned
%global optional_ext_file_flag file_ext_string
%added by EM: 04.10.2016 
%global ACDC_application_cnf cnf_p_ACDC
%global netcdf_type
%global processing_mode_cnf
% t6 = tic;

% date_creation = datestr(now, '_yyyymmddTHHMMSS_');


if strcmp(files.product_type, 'L1B_HR')
    L1B_date_find = strfind(files.filename_L1B_HR,string('_1B_'));
    L1B_mode_find = strfind(files.filename_L1B_HR, 'PICE_SIRS_');
    files.filename_L2GR = strcat(files.outputPath,'PICE_SIRS_',...
        files.filename_L1B_HR((L1B_mode_find(1)+9):(L1B_mode_find(1)+15)),'_L2_GR_S_',...
        files.filename_L1B_HR((L1B_date_find(1)+4):(L1B_date_find(1)+34)),...
        '_isd','.nc');
elseif strcmp(files.product_type, 'L1B_FF')
    L1B_date_find = strfind(files.filename_L1BFF_ML,string('_1B_'));
    L1B_mode_find = strfind(files.filename_L1BFF_ML, 'PICE_SIRS_');
    files.filename_L2GR = strcat(files.outputPath,'PICE_SIRS_',...
        files.filename_L1BFF_ML((L1B_mode_find(1)+9):(L1B_mode_find(1)+21)),'_L2_GR_S_',...
        files.filename_L1BFF_ML((L1B_date_find(1)+4):(L1B_date_find(1)+34)),...
        '_isd','.nc');
end

% L1B_date_find = strfind(files.filename_L1B,'_201');
% 
% files.filename_L2GR_Swath = strcat(files.outputPath,'SR_2_SRA____',...
%         files.filename_L1A((L1A_date_find(1)-8):(L1A_date_find(1)-3)),'L2_',...
%         files.filename_L1A((L1A_date_find(1)+1):(L1A_date_find(1)+30)),...
%         '_isd','.nc');


files.filename_L2GR = 'prova1';
N_rec = 5;
N_peaks = 5;
ncid = netcdf.create(files.filename_L2GR_Swath,'NETCDF4');
data_id = netcdf.defGrp(ncid,'data');
data_ku_id = netcdf.defGrp(data_id,'ku');

long_name_att = 'long_name';
std_name_att = 'standard_name';
calendar_name_att='calendar';
comment_att = 'comment';
units_att = 'units';
scale_factor_att = 'scale_factor';
add_offset_att = 'add_offset';
flag_values_att='flag_values';
flag_desc_att='flag_meanings';

nr_dimension = netcdf.defDim(data_ku_id,'nr',N_rec);
space_3D_dimension = netcdf.defDim(data_ku_id,'space_3D',3);
max_swath_length_dimension = netcdf.defDim(data_ku_id,'max_swath_length',max_swath_length);


day_units = 'day';
seconds_units = 'seconds';
meters_units = 'meters';
meters_per_second_units = 'm/s';
degrees_north_units = 'degrees_north';
degrees_east_units = 'degrees_east';
degrees_units = 'degrees';
dB_units = 'dB';
per_second_units = '/s';

int16_type = 'NC_SHORT';
int32_type = 'NC_INT';
double_type= 'NC_DOUBLE';
int8_type = 'NC_BYTE';

%% PACKING L2 GR Swath


%---------- Time ----------------------------------------------
time = 'time';
id_aux = netcdf.defVar(data_ku_id,time,double_type,nr_dimension);
netcdf.putAtt(data_ku_id,id_aux,std_name_att,'time');
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'Time in UTC');
netcdf.putAtt(data_ku_id,id_aux,calendar_name_att,'Gregorian');
netcdf.putAtt(data_ku_id,id_aux,units_att,seconds_units);
netcdf.putAtt(data_ku_id,id_aux,comment_att,'It contains the seconds since 1 Jan 2000 00:00:00. Time is in UTC.');

time_tai = 'time_tai';
id_aux = netcdf.defVar(data_ku_id,time_tai,double_type,nr_dimension);
netcdf.putAtt(data_ku_id,id_aux,std_name_att,'time');
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'Time in TAI');
netcdf.putAtt(data_ku_id,id_aux,calendar_name_att,'Gregorian');
netcdf.putAtt(data_ku_id,id_aux,units_att,seconds_units);
netcdf.putAtt(data_ku_id,id_aux,comment_att,'Time of measurement in seconds in the TAI time scale since 1 Jan 2000 00:00:00 TAI.');

%---------- Swath retrievals variables ------------------------------
length_swath_name = 'length_swath';
id_aux = netcdf.defVar(data_ku_id, length_swath_name,int8_type,nr_dimension);
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'number of points across track');
netcdf.putAtt(data_ku_id,id_aux,comment_att,'Number of across track points in the strip.');

bin_index_name = 'bin_index';
id_aux = netcdf.defVar(data_ku_id,bin_index_name,int16_type,[nr_dimension, max_swath_length_dimension]);
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'sample index from the L1B product of the retrieved measurement');
netcdf.putAtt(data_ku_id,id_aux,comment_att,'Sample index from the L1B product of the retrieved measurement, taken from the l1b_record_counter variable.');

height_swath_name = 'height_swath';
id_aux = netcdf.defVar(data_ku_id,height_swath_name,int32_type,[nr_dimension, max_swath_length_dimension]);
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'height of the surface at the measurement point w.r.t. the reference ellipsoid');
netcdf.putAtt(data_ku_id,id_aux,units_att,meters_units);
netcdf.putAtt(data_ku_id,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_ku_id,id_aux,add_offset_att,700000);
netcdf.putAtt(data_ku_id,id_aux,comment_att,'Measured height of the surface above the reference ellipsoid (WGS84) at each across track point location');

latitude_swath_name = 'latitude_swath';
id_aux = netcdf.defVar(data_ku_id,latitude_swath_name,int32_type,[nr_dimension, max_swath_length_dimension]);
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'latitude of each strip point');
netcdf.putAtt(data_ku_id,id_aux,std_name_att,'latitude');
netcdf.putAtt(data_ku_id,id_aux,units_att,degrees_north_units);
netcdf.putAtt(data_ku_id,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_ku_id,id_aux,comment_att,'Latitude of each across track points [-90, +90]. Positive latitude is North latitude, negative latitude is South latitude.');

longitude_swath_name = 'longitude_swath';
id_aux = netcdf.defVar(data_ku_id,longitude_swath_name,int32_type,[nr_dimension, max_swath_length_dimension]);
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'longitude of each strip point');
netcdf.putAtt(data_ku_id,id_aux,std_name_att,'longitude');
netcdf.putAtt(data_ku_id,id_aux,units_att,degrees_east_units);
netcdf.putAtt(data_ku_id,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_ku_id,id_aux,comment_att,'Longitude of each across track points [0, 360). East longitude relative to Greenwich meridian.');

angle_of_arrival_unwrapped_name = 'angle_of_arrival_unwrapped';
id_aux = netcdf.defVar(data_ku_id,angle_of_arrival_unwrapped_name,int32_type,[nr_dimension, max_swath_length_dimension]);
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'angle of arrival retrieved after the unwrapping process');
netcdf.putAtt(data_ku_id,id_aux,units_att,degrees_east_units);
netcdf.putAtt(data_ku_id,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_ku_id,id_aux,comment_att,'Angle of arrival at each across track point, computed from the DEM model in SAR mode and from the phase difference in SARin mode.');

coherence_name = 'coherence';
id_aux = netcdf.defVar(data_ku_id,coherence_name,int32_type,[nr_dimension, max_swath_length_dimension]);
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'coherence of the SARIn echoes at the retracking point');
netcdf.putAtt(data_ku_id,id_aux,scale_factor_att,1.e-3);
netcdf.putAtt(data_ku_id,id_aux,comment_att,'Not used in SAR mode.');

dem_diff_name = 'dem_diff';
id_aux = netcdf.defVar(data_ku_id,dem_diff_name,int32_type,[nr_dimension, max_swath_length_dimension]);
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'difference w.r.t. the reference Digital Elevation Modelt');
netcdf.putAtt(data_ku_id,id_aux,units_att,meters_units);
netcdf.putAtt(data_ku_id,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_ku_id,id_aux,comment_att,'Difference w.r.t. the reference DEM at each across track point.');

dem_diff_mad_name = 'dem_diff_mad';
id_aux = netcdf.defVar(data_ku_id,dem_diff_mad_name,int32_type,[nr_dimension, max_swath_length_dimension]);
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'median absolute deviation w.r.t. the reference digital elevation model');
netcdf.putAtt(data_ku_id,id_aux,units_att,meters_units);
netcdf.putAtt(data_ku_id,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_ku_id,id_aux,comment_att,'Median absolute deviation w.r.t. the reference digital elevation model at each across track point.');

%---------- Land Ice geophysical retrievals variables ------------------------------
land_ice_elevation_name = 'land_ice_elevation';
id_aux = netcdf.defVar(data_ku_id,land_ice_elevation_name,int32_type,nr_dimension);
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'land ice elevation');
netcdf.putAtt(data_ku_id,id_aux,units_att,meters_units);
netcdf.putAtt(data_ku_id,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_ku_id,id_aux,comment_att,'Measured height of the surface above the reference ellipsoid (WGS84) at each measurement location. All instrumental and land ice geophysical corrections included.');

%---------- Inland waters retrievals variables ------------------------------
water_level_height_name = 'water_level_height';
id_aux = netcdf.defVar(data_ku_id,water_level_height_name,int32_type,nr_dimension);
netcdf.putAtt(data_ku_id,id_aux,long_name_att,'water level height');
netcdf.putAtt(data_ku_id,id_aux,units_att,meters_units);
netcdf.putAtt(data_ku_id,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_ku_id,id_aux,comment_att,'Measured height of the surface above the reference ellipsoid (WGS84) at each measurement location. All instrumental and inland water geophysical corrections included.');

% %----------  Global Attributes definition -----------------------------------
% %---- attributes inherited from Sentinel-3 product description-------------
% id_aux = netcdf.getConstant('NC_GLOBAL');
% netcdf.putAtt(data_ku_id,id_aux,'creation_time',date_creation);
% netcdf.putAtt(data_ku_id,id_aux,'Conventions',netcdf_v4_format);
% netcdf.putAtt(data_ku_id,id_aux,'altimeter_sensor_name',altimeter_sensor_name);
% netcdf.putAtt(data_ku_id,id_aux,'first_meas_time',first_meas_time);
% netcdf.putAtt(data_ku_id,id_aux,'last_meas_time',last_meas_time);
% netcdf.putAtt(data_ku_id,id_aux,'first_meas_lat',first_meas_lat);
% netcdf.putAtt(data_ku_id,id_aux,'last_meas_lat',last_meas_lat);
% netcdf.putAtt(data_ku_id,id_aux,'first_meas_lon',first_meas_lon);
% netcdf.putAtt(data_ku_id,id_aux,'last_meas_lon',last_meas_lon);
% netcdf.putAtt(data_ku_id,id_aux,'semi_major_ellipsoid_axis',semi_major_ellipsoid_axis);
% netcdf.putAtt(data_ku_id,id_aux,'ellipsoid_flattening',ellipsoid_flattening);
% %--------------- add the attributes related to intermediate product--------
% netcdf.putAtt(data_ku_id,id_aux,'orbit_cycle_num',orbit_cycle_num);
% netcdf.putAtt(data_ku_id,id_aux,'orbit_REL_Orbit',orbit_REL_Orbit);
% netcdf.putAtt(data_ku_id,id_aux,'orbit_ABS_Orbit_Start',orbit_ABS_Orbit_Start);

netcdf.endDef(ncid);
netcdf.close(ncid);
% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);

end
