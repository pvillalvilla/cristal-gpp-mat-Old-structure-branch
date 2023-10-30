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
function [files,ncid] = create_NetCDF_L2GR0(files, N_rec, cnf, chd, cst)

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

if strcmp(files.product_type, 'L1B_LR')

elseif strcmp(files.product_type, 'L1B_LRRMC')

elseif strcmp(files.product_type, 'L1B_LROS')

elseif strcmp(files.product_type, 'L1B_HR')
    L1B_date_find = strfind(files.filename_L1B_HR,string('_1B_'));
    L1B_mode_find = strfind(files.filename_L1B_HR, 'PICE_SIRS_');
    files.filename_L2GR = strcat(files.outputPath,'PICE_SIRS_',...
        files.filename_L1B_HR((L1B_mode_find(1)+9):(L1B_mode_find(1)+15)),'_L2_GR_',...
        files.filename_L1B_HR((L1B_date_find(1)+4):(L1B_date_find(1)+34)),...
        '_isd','.nc');
elseif strcmp(files.product_type, 'L1B_FF')
    L1B_date_find = strfind(files.filename_L1BFF_ML,string('_1B_'));
    L1B_mode_find = strfind(files.filename_L1BFF_ML, 'PICE_SIRS_');
    files.filename_L2GR = strcat(files.outputPath,'PICE_SIRS_',...
        files.filename_L1BFF_ML((L1B_mode_find(1)+9):(L1B_mode_find(1)+21)),'_L2_GR_',...
        files.filename_L1BFF_ML((L1B_date_find(1)+4):(L1B_date_find(1)+34)),...
        '_isd','.nc');
end
L1B_date_find = strfind(files.filename_L1BFF_ML,string('_1B_'));
L1B_mode_find = strfind(files.filename_L1BFF_ML, 'PICE_SIRS_');
files.filename_L2GR = strcat(files.outputPath,'PICE_SIRS_',...
        files.filename_L1BFF_ML((L1B_mode_find(1)+9):(L1B_mode_find(1)+21)),'_L2_GR_',...
        files.filename_L1BFF_ML((L1B_date_find(1)+4):(L1B_date_find(1)+34)),...
        '_isd','.nc');

files.filename_L2GR = 'prova3';
N_rec = 5;
N_peaks = 5;
ncid = netcdf.create(files.filename_L2GR,'NETCDF4');
% data_id = netcdf.defGrp(ncid,'data');
% ncid = netcdf.defGrp(data_id,'ku');
% data_ka_id = netcdf.defGrp(data_id,'ka');

long_name_att = 'long_name';
std_name_att = 'standard_name';
calendar_name_att='calendar';
comment_att = 'comment';
units_att = 'units';
scale_factor_att = 'scale_factor';
add_offset_att = 'add_offset';
flag_values_att='flag_values';
flag_desc_att='flag_meanings';

nr_dimension = netcdf.defDim(ncid,'nr',N_rec);
space_3D_dimension = netcdf.defDim(ncid,'space_3D',3);

nr_dimension_ku = netcdf.defDim(ncid,'nr_ku',N_rec);
np_dimension_ku = netcdf.defDim(ncid,'np_ku',N_peaks);
space_3D_dimension_ku = netcdf.defDim(ncid,'space_3D_ku',3);

% nr_dimension_ka = netcdf.defDim(data_ka_id,'nr',N_rec);
% np_dimension_ka = netcdf.defDim(data_ka_id,'np',N_peaks);

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

%% PACKING L2


%---------- Time and counter variables ----------------------------------------------
time_l2_echo_name = 'time';
id_aux = netcdf.defVar(ncid,time_l2_echo_name,double_type,nr_dimension);
netcdf.putAtt(ncid,id_aux,std_name_att,'time');
netcdf.putAtt(ncid,id_aux,long_name_att,'Time in UTC');
netcdf.putAtt(ncid,id_aux,calendar_name_att,'Gregorian');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,comment_att,'It contains the seconds since 1 Jan 2000 00:00:00. Time is in UTC.');

time_tai_l2_echo_name = 'time_tai';
id_aux = netcdf.defVar(ncid,time_tai_l2_echo_name,double_type,nr_dimension);
netcdf.putAtt(ncid,id_aux,std_name_att,'time');
netcdf.putAtt(ncid,id_aux,long_name_att,'Time in TAI');
netcdf.putAtt(ncid,id_aux,calendar_name_att,'Gregorian');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Time of measurement in seconds in the TAI time scale since 1 Jan 2000 00:00:00 TAI.');


L2_record_counter_name = 'l2_record_counter';
id_aux = netcdf.defVar(ncid,L2_record_counter_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,std_name_att,'records');
netcdf.putAtt(ncid,id_aux,long_name_att,'L2 record counter');
netcdf.putAtt(ncid,id_aux,comment_att,'L2 record counter');

L1B_record_counter_name = 'l1b_record_counter';
id_aux = netcdf.defVar(ncid,L1B_record_counter_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,std_name_att,'records');
netcdf.putAtt(ncid,id_aux,long_name_att,'L1B record counter');
netcdf.putAtt(ncid,id_aux,comment_att,'L1B record counter');

tm_source_sequence_counter_name = 'tm_source_sequence_counter';
id_aux = netcdf.defVar(ncid,tm_source_sequence_counter_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,std_name_att,'records');
netcdf.putAtt(ncid,id_aux,long_name_att,'instrument sequence counter of the Closest Burst');
netcdf.putAtt(ncid,id_aux,comment_att,'Instrument record counter');

%---------- Orbit and attitude variables ------------------------------
altitude_name = 'altitude';
id_aux = netcdf.defVar(ncid, altitude_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'center of mass altitude above the reference ellipsoid');
netcdf.putAtt(ncid,id_aux,std_name_att,'height_above_reference_ellipsoid');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Altitude of the satellite center of mass above the reference ellipsoid (WGS84)');

altitude_rate_name = 'altitude_rate';
id_aux = netcdf.defVar(ncid,altitude_rate_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'center of mass altitude rate with respect to the reference ellipsoid');
netcdf.putAtt(ncid,id_aux,std_name_att,'derivative_of_height_above_reference_ellipsoid_wrt_time');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Instantaneous altitude rate of the satellite center of mass with respect to the reference ellipsoid (WGS84)');

Latitude_sat_name = 'latitude';
id_aux = netcdf.defVar(ncid,Latitude_sat_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'latitude of the satellite center of mass at time stamp');
netcdf.putAtt(ncid,id_aux,std_name_att,'latitude');
netcdf.putAtt(ncid,id_aux,units_att,degrees_north_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Latitude of measurement [-90, +90]. Positive latitude is North latitude, negative latitude is South latitude');

longitude_name = 'longitude';
id_aux = netcdf.defVar(ncid,longitude_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'longitude of the satellite center of mass at time stamp');
netcdf.putAtt(ncid,id_aux,std_name_att,'longitude');
netcdf.putAtt(ncid,id_aux,units_att,degrees_east_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Longitude of measurement [0, 360). East longitude relative to Greenwich meridian');

latitude_surf_name = 'latitude_surf';
id_aux = netcdf.defVar(ncid,latitude_surf_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'latitude of the measurement at the surface');
netcdf.putAtt(ncid,id_aux,std_name_att,'latitude');
netcdf.putAtt(ncid,id_aux,units_att,degrees_north_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Latitude of point of closest approach [-90, +90]. Positive latitude is North latitude, negative latitude is South latitude');

longitude_surf_name = 'longitude_surf';
id_aux = netcdf.defVar(ncid,longitude_surf_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'longitude of the measurement at the surface');
netcdf.putAtt(ncid,id_aux,std_name_att,'longitude');
netcdf.putAtt(ncid,id_aux,units_att,degrees_east_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Longitude of point of closest approach [0, 360). East longitude relative to Greenwich meridian');

position_vector_name = 'position_vector';
id_aux = netcdf.defVar(ncid,position_vector_name,double_type,[space_3D_dimension_ku ,nr_dimension_ku]);
netcdf.putAtt(ncid,id_aux,long_name_att,'center of mass position vector in ITRF, components: [1] x, [2] y, [3] z');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Position vector at the center of mass in ITRF');

velocity_vector_name = 'velocity_vector';
id_aux = netcdf.defVar(ncid,velocity_vector_name,int32_type,[space_3D_dimension_ku ,nr_dimension_ku]);
netcdf.putAtt(ncid,id_aux,long_name_att,'center of mass velocity vector in ITRF, components: [1] x, [2] y, [3] z');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Velocity vector at the center of mass in ITRF');

off_nadir_pitch_angle_name = 'off_nadir_pitch_angle';
id_aux = netcdf.defVar(ncid,off_nadir_pitch_angle_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'off-nadir pitch angle corresponding to the waveform used in the retracking');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-7);
netcdf.putAtt(ncid,id_aux,comment_att,'Pitch angle, corresponding to the waveform used in the retracking, with respect to the nadir pointing, measured by the STRs and post-processed by AOCS or by ground facility');

off_nadir_roll_angle_name = 'off_nadir_roll_angle';
id_aux = netcdf.defVar(ncid,off_nadir_roll_angle_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'off-nadir roll angle corresponding to the waveform used in the retracking');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-7);
netcdf.putAtt(ncid,id_aux,comment_att,'Roll angle, corresponding to the waveform used in the retracking, with respect to the nadir pointing, measured by the STRs and post-processed by AOCS or by ground facility');

off_nadir_yaw_angle_name = 'off_nadir_yaw_angle';
id_aux = netcdf.defVar(ncid,off_nadir_yaw_angle_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'off-nadir yaw angle corresponding to the waveform used in the retracking');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-7);
netcdf.putAtt(ncid,id_aux,comment_att,'Yaw angle, corresponding to the waveform used in the retracking, with respect to the nadir pointing, measured by the STRs and post-processed by AOCS or by ground facility');

%---------- Configuration Variables --------------------------------------
iris_mode_flag_name = 'iris_mode_flag';
id_aux = netcdf.defVar(ncid,iris_mode_flag_name,int8_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'mode identifier from the IRIS instrument source packets');
netcdf.putAtt(ncid,id_aux,std_name_att,'status_flag');
netcdf.putAtt(ncid,id_aux, flag_values_att,'0, 1, 2, 3, 4, 5');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'SAR_CB_RMC SAR_CB_RAW  SAR_CB_LRM SARin_CB SARin_CB_OL SARin_OB');
netcdf.putAtt(ncid,id_aux,comment_att,'Flag containing the mode identifier from the IRIS instrument source packets');

telemetry_type_flag_name = 'telemetry_type_flag';
id_aux = netcdf.defVar(ncid,telemetry_type_flag_name,int8_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'telemetry type: RAW or RMC');
netcdf.putAtt(ncid,id_aux,std_name_att,'status_flag');
netcdf.putAtt(ncid,id_aux, flag_values_att,'0, 1');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'RAW RMC');
netcdf.putAtt(ncid,id_aux,comment_att,'Telemetry packet type can be: RAW or RMC.');

%---------- Altimeter range variables --------------------------------------
tracker_range_calibrated_name = 'tracker_range_calibrated';
id_aux = netcdf.defVar(ncid,tracker_range_calibrated_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'calibrated 1-way tracker range');
netcdf.putAtt(ncid,id_aux,std_name_att,'altimeter_range');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'This is the 1-way distance from the satellite center of mass to the middle of the range window (sample ns/2 from 0)');


%---------- Altimeter power variables --------------------------------------

sig0_scaling_factor_name = 'sig0_scaling_factor';
id_aux = netcdf.defVar(ncid,sig0_scaling_factor_name,int16_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'scaling factor for the sig0 evaluation (from the antenna flange)');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,comment_att,'The scaling factor in order to retrieve sigma-0 from L1B waveform');

%---------- Flag Variables --------------------------------------
surface_classification_flag_name = 'surface_classification_flag';
id_aux = netcdf.defVar(ncid,surface_classification_flag_name,int8_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'surface type classification');
netcdf.putAtt(ncid,id_aux,std_name_att,'status_flag');
netcdf.putAtt(ncid,id_aux, flag_values_att,'0, 1, 2, 3, 4, 5, 6');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'open_ocean land continental_water aquatic_vegetation continental_ice_snow floating_ice salted_basin');
netcdf.putAtt(ncid,id_aux,comment_att,'7-state surface type classification from MODIS and GlobCover.');

%---------- Geophysical corrections variables --------------------------------------
rad_wet_tropo_cor_name = 'rad_wet_tropo_cor';
id_aux = netcdf.defVar(ncid,rad_wet_tropo_cor_name,int16_type,nr_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'radiometer wet troposphere correction');
netcdf.putAtt(ncid,id_aux,std_name_att,'altimeter_range_correction_due_to_wet_troposphere');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Wet tropospheric correction derived from the on-board radiometer and valid over ocean surfaces only.');

rad_atm_cor_sig0_name = 'rad_atm_cor_sig0';
id_aux = netcdf.defVar(ncid,rad_atm_cor_sig0_name,int16_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'two-way atmospheric attenuation correction on backscatter coefficient');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,comment_att,'Two-way atmospheric attenuation to backscatter coefficient (sig0) based on the radiometer');

ocean_geo_corrections_name = 'ocean_geo_corrections';
id_aux = netcdf.defVar(ncid,ocean_geo_corrections_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'sum of meteo and geo corrections to be applied over ocean');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Sum of meteorological and geophysical corrections to be applied to correct the range over ocean, semi-enclosed sea and extended lakes.');

seaice_geo_corrections_name = 'seaice_geo_corrections';
id_aux = netcdf.defVar(ncid,seaice_geo_corrections_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'sum of meteo and geo corrections to be applied over sea-ice');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Sum of meteorological and geophysical corrections to be applied to correct the range over sea-ice areas.');

lead_geo_corrections_name = 'lead_geo_corrections';
id_aux = netcdf.defVar(ncid,lead_geo_corrections_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'sum of meteo and geo corrections to be applied over lead');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Sum of meteorological and geophysical corrections to be applied to correct the range over leads areas.');

iceshelve_geo_corrections_name = 'iceshelve_geo_corrections';
id_aux = netcdf.defVar(ncid,iceshelve_geo_corrections_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'sum of meteo and geo corrections to be applied over ice shelve');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Sum of meteorological and geophysical corrections to be applied to correct the range over ice shelve areas.');

landice_geo_corrections_name = 'landice_geo_corrections';
id_aux = netcdf.defVar(ncid,landice_geo_corrections_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'sum of meteo and geo corrections to be applied over land ice');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Sum of meteorological and geophysical corrections to be applied to correct the range over land ice areas or ice covered lakes.');

inland_geo_corrections_name = 'inland_geo_corrections';
id_aux = netcdf.defVar(ncid,inland_geo_corrections_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'sum of meteo and geo corrections to be applied over inland water');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Sum of meteorological and geophysical corrections to be applied to correct the range over inland water areas.');

%---------- Reference surfaces variables --------------------------------------

geoid_name = 'geoid';
id_aux = netcdf.defVar(ncid,geoid_name,int32_type,nr_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'geoid height');
netcdf.putAtt(ncid,id_aux,std_name_att,'geoid_height_above_reference_ellipsoid');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Height of the geoid above the reference ellipsoid (WGS84).');

mean_sea_surface_name = 'mean_sea_surface';
id_aux = netcdf.defVar(ncid,mean_sea_surface_name,int32_type,nr_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'mean sea surface from model');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Mean sea surface height above reference ellipsoid (WGS84).');

%---------- Waveform characteristics variables --------------------------------------

peakiness_name = 'peakiness';
id_aux = netcdf.defVar(ncid,peakiness_name,int16_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'pulse peakiness');
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,comment_att,'Peakiness of the waveform.');

noise_floor_name = 'noise_floor';
id_aux = netcdf.defVar(ncid,noise_floor_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'noise floor');
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Noise floor estimate.');

number_of_peaks_name = 'number_of_peaks';
id_aux = netcdf.defVar(ncid,number_of_peaks_name,int8_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'number of peaks');
netcdf.putAtt(ncid,id_aux,comment_att,'Number of peaks within the waveform.');

prominence_name = 'prominence';
id_aux = netcdf.defVar(ncid,prominence_name,int16_type,[nr_dimension_ku, np_dimension_ku]);
netcdf.putAtt(ncid,id_aux,long_name_att,'prominence of each peak');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,comment_att,'Prominence of each peak within the waveform.');

slope_trailing_edge_name = 'slope_trailing_edge';
id_aux = netcdf.defVar(ncid,slope_trailing_edge_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'slope of the trailing edge');
netcdf.putAtt(ncid,id_aux,units_att,per_second_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Slope of the trailing edge of the waveform.');

peak_noise_ratio_name = 'peak_noise_ratio';
id_aux = netcdf.defVar(ncid,peak_noise_ratio_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'peak to noise ratio');
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Ratio between the peak of the waveform and the noise floor.');

%---------- Retracking estimates variables --------------------------------------

epoch_ocog_name = 'epoch_ocog';
id_aux = netcdf.defVar(ncid,epoch_ocog_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'epoch derived from the OCOG retracker');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
netcdf.putAtt(ncid,id_aux,comment_att,'Retracking offset (in terms of 2-way time delay) from the middle of the range window derived from the OCOG retracker.');

epoch_tcog_name = 'epoch_tcog';
id_aux = netcdf.defVar(ncid,epoch_tcog_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'epoch derived from the TCOG retracker');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
netcdf.putAtt(ncid,id_aux,comment_att,'Retracking offset (in terms of 2-way time delay) from the middle of the range window derived from the TCOG retracker.');

epoch_tfmra_name = 'epoch_tfmra';
id_aux = netcdf.defVar(ncid,epoch_tfmra_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'epoch derived from the TFMRA retracker');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
netcdf.putAtt(ncid,id_aux,comment_att,'Retracking offset (in terms of 2-way time delay) from the middle of the range window derived from the TFMRA retracker.');

epoch_name = 'epoch';
id_aux = netcdf.defVar(ncid,epoch_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'epoch derived from the two-step retracker');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
netcdf.putAtt(ncid,id_aux,comment_att,'Retracking offset (in terms of 2-way time delay) from the middle of the range window derived from the physical retracker.');

amplitude_ocog_name = 'amplitude_ocog';
id_aux = netcdf.defVar(ncid,amplitude_ocog_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'waveform amplitude derived from the OCOG retracker');
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Waveform amplitude derived from the OCOG retracker.');

amplitude_tcog_name = 'amplitude_tcog';
id_aux = netcdf.defVar(ncid,amplitude_tcog_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'waveform amplitude derived from the TCOG retracker');
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Waveform amplitude derived from the TCOG retracker.');

amplitude_tfmra_name = 'amplitude_tfmra';
id_aux = netcdf.defVar(ncid,amplitude_tfmra_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'waveform amplitude derived from the TFMRA retracker');
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Waveform amplitude derived from the TFMRA retracker.');

amplitude_name = 'amplitude';
id_aux = netcdf.defVar(ncid,amplitude_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'waveform amplitude derived from the two-step retracker');
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Waveform amplitude derived from the physical retracker.');

%---------- Retracking retrievals variables --------------------------------------

range_ocog_name = 'range_ocog';
id_aux = netcdf.defVar(ncid,range_ocog_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'corrected altimeter range from the OCOG retracker');
netcdf.putAtt(ncid,id_aux,std_name_att,'altimeter_range');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Range computed from the OCOG retracker. It includes all instrumental corrections.');

range_tcog_name = 'range_tcog';
id_aux = netcdf.defVar(ncid,range_tcog_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'corrected altimeter range from the TCOG retracker');
netcdf.putAtt(ncid,id_aux,std_name_att,'altimeter_range');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Range computed from the TCOG retracker. It includes all instrumental corrections.');

range_tfmra_name = 'range_tfmra';
id_aux = netcdf.defVar(ncid,range_tfmra_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'corrected altimeter range from the TFMRA retracker');
netcdf.putAtt(ncid,id_aux,std_name_att,'altimeter_range');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Range computed from the TFMRA retracker. It includes all instrumental corrections.');

range_name = 'range';
id_aux = netcdf.defVar(ncid,range_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'corrected altimeter range from the two-step retracker');
netcdf.putAtt(ncid,id_aux,std_name_att,'altimeter_range');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Range computed from the physical retracker. It includes all instrumental corrections.');

sig0_ocog_name = 'sig0_ocog';
id_aux = netcdf.defVar(ncid,sig0_ocog_name,int16_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'corrected backscatter coefficient from the OCOG retracker');
netcdf.putAtt(ncid,id_aux,std_name_att,'surface_backwards_scattering_coefficient_of_radar_wave');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,comment_att,'Fully corrected backscatter from the OCOG retracker including instrument gain correction and bias.');

sig0_tcog_name = 'sig0_tcog';
id_aux = netcdf.defVar(ncid,sig0_tcog_name,int16_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'corrected backscatter coefficient from the TCOG retracker');
netcdf.putAtt(ncid,id_aux,std_name_att,'surface_backwards_scattering_coefficient_of_radar_wave');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,comment_att,'Fully corrected backscatter from the TCOG retracker including instrument gain correction and bias.');

sig0_tfmra_name = 'sig0_tfmra';
id_aux = netcdf.defVar(ncid,sig0_tfmra_name,int16_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'corrected backscatter coefficient from the TFMRA retracker');
netcdf.putAtt(ncid,id_aux,std_name_att,'surface_backwards_scattering_coefficient_of_radar_wave');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,comment_att,'Fully corrected backscatter from the TFMRA retracker including instrument gain correction and bias.');

sig0_name = 'sig0';
id_aux = netcdf.defVar(ncid,sig0_name,int16_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'corrected backscatter coefficient from the two-step retracker');
netcdf.putAtt(ncid,id_aux,std_name_att,'surface_backwards_scattering_coefficient_of_radar_wave');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,comment_att,'Fully corrected backscatter from the physical retracker including instrument gain correction and bias.');

%---------- Ocean geophysical retrievals variables --------------------------------------

ssha_name = 'ssha';
id_aux = netcdf.defVar(ncid,ssha_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'sea surface height from the ocean retracker');
netcdf.putAtt(ncid,id_aux,std_name_att,'sea_surface_height_above_mean_sea_level');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Sea surface height anomaly');

swh_name = 'swh';
id_aux = netcdf.defVar(ncid,swh_name,int16_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'significant wave height from the two-step retracker');
netcdf.putAtt(ncid,id_aux,std_name_att,'sea_surface_wave_significant_height');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'SWH computed from the physical retracker.');

wind_speed_alt_name = 'wind_speed_alt';
id_aux = netcdf.defVar(ncid,wind_speed_alt_name,int16_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'wind speed over ocean');
netcdf.putAtt(ncid,id_aux,std_name_att,'wind_speed');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Wind speed at 10 m height derived from altimeter measurements.');

%---------- Snow geophysical retrievals variables --------------------------------------

snow_depth_name = 'snow_depth';
id_aux = netcdf.defVar(ncid,snow_depth_name,int16_type,nr_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'depth of the snow from Ku and Ka band measurements');
netcdf.putAtt(ncid,id_aux,std_name_att,'surface_snow_thickness');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Depth of the snow retrieved with Ku and Ka band measurements, based on the different scattering horizon sensed by the altimeter at Ku and Ka bands within the snowpack. Set to zero if t if Ku and Ka are not co-located at L1B.');

snow_depth_correction_name = 'snow_depth_correction';
id_aux = netcdf.defVar(ncid,snow_depth_correction_name,int16_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'additional delay produced by the snow');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Additional delay due to the snow computed from the snow thickness. Applicable only to Ku-band.');

snow_interface_alignement_name = 'snow_interface_alignement';
id_aux = netcdf.defVar(ncid,snow_interface_alignement_name,int16_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'snow interface alignment');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Snow interface alignement retrieved from the pulse peakiness.Alignment to the snow-ice interface in Ku band, alignement to the snow surface in Ka band.');

%---------- Sea-ice geophysical retrievals variables --------------------------------------

sea_ice_concentration_name = 'sea_ice_concentration';
id_aux = netcdf.defVar(ncid,sea_ice_concentration_name,int16_type,nr_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'sea ice concentration from model');
netcdf.putAtt(ncid,id_aux,std_name_att,'sea_ice_area_fraction');
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,comment_att,'Sea ice concentration from NSIDC. Percentage between 0 and 100.');

sea_ice_type_name = 'sea_ice_type';
id_aux = netcdf.defVar(ncid,sea_ice_type_name,int8_type,nr_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'sea ice type');
netcdf.putAtt(ncid,id_aux,std_name_att,'status_flag');
netcdf.putAtt(ncid,id_aux, flag_values_att,'1, 2, 3, 4');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'open_water, first_year_ice, multi_year_ice, ambiguous');
netcdf.putAtt(ncid,id_aux,comment_att,'Sea ice type from OSISAF.');

waveform_classification_flag_name = 'waveform_classification_flag';
id_aux = netcdf.defVar(ncid,waveform_classification_flag_name,int8_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'waveform classification in sea ice domain');
netcdf.putAtt(ncid,id_aux,std_name_att,'status_flag');
netcdf.putAtt(ncid,id_aux, flag_values_att,'0, 1, 2, 3');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'lead, ocean, sea_ice, complex');
netcdf.putAtt(ncid,id_aux,comment_att,'Waveform classification determined from shape features of the waveforms and sea ice concentration.');

sea_ice_freeboard_name = 'sea_ice_freeboard';
id_aux = netcdf.defVar(ncid,sea_ice_freeboard_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'difference between the height of the surface of sea ice and the water in open leads');
netcdf.putAtt(ncid,id_aux,std_name_att,'sea_ice_freeboard');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Computed as sea_ice_freeboard = altitude of satellite (altitude) - range from retracker (range_tfmra or range, depending on selection through configuration parameter) - sea-ice geo corrections (seaice_geo_corrections) - mean sea surface (mean_sea_surface) - ssha interpolated (ssha_lead_interpolated).');

sea_ice_thickness_name = 'sea_ice_thickness';
id_aux = netcdf.defVar(ncid,sea_ice_thickness_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'derived thickness of the sea ice');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Computed from the sea-ice freeboard and the snow depth.');

ssha_lead_interpolated_name = 'ssha_lead_interpolated';
id_aux = netcdf.defVar(ncid,ssha_lead_interpolated_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'interpolated sea surface height anomaly');
netcdf.putAtt(ncid,id_aux,std_name_att,'sea_surface_height_above_mean_sea_level');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Interpolated at sea-ice freeboard location from ssha over ocean and leads.');

iceberg_detection_flag_name = 'iceberg_detection_flag';
id_aux = netcdf.defVar(ncid,iceberg_detection_flag_name,int8_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'iceberg_detection_flag');
netcdf.putAtt(ncid,id_aux, flag_values_att,'0, 1');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'no_iceberg, high_iceberg_probability');
netcdf.putAtt(ncid,id_aux,comment_att,'Iceberg detection flag determined from the thermal noise part of the HR waveform.');

iceberg_freeboard_name = 'iceberg_freeboard';
id_aux = netcdf.defVar(ncid,iceberg_freeboard_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'iceberg freeboard');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Iceberg freeboard determined from the range of the iceberg and from the angle of arrival in SARIN mode.');

iceberg_sig0_name = 'iceberg_sig0';
id_aux = netcdf.defVar(ncid,iceberg_sig0_name,int16_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'iceberg backscatter coefficient');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,comment_att,'Mean sigma0 of the pixels which contain the iceberg candidate.');

%---------- Land Ice geophysical retrievals variables --------------------------------------

land_ice_elevation_name = 'land_ice_elevation';
id_aux = netcdf.defVar(ncid,land_ice_elevation_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'land ice elevation');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Measured height of the surface above the reference ellipsoid (WGS84) at each measurement location. All instrumental and land ice geophysical corrections included.');

%---------- Inland waters retrievals variables --------------------------------------

land_ice_elevation_name = 'water_level_height';
id_aux = netcdf.defVar(ncid,land_ice_elevation_name,int32_type,nr_dimension_ku);
netcdf.putAtt(ncid,id_aux,long_name_att,'water level height');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,comment_att,'Measured height of the surface above the reference ellipsoid (WGS84) at each measurement location. All instrumental and inland water geophysical corrections included.');


netcdf.endDef(ncid);
netcdf.close(ncid);

end
