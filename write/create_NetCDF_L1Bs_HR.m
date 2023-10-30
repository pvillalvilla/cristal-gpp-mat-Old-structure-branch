%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% 
% This code implements the CODING & PACKING 
% algorithm for Level-1BS HR CRISTAL product 
%
% ---------------------------------------------------------
% Objective: Pack variables and write the NETCDF
% 
% INPUTs : Workspace
% OUTPUTs: NETCDF file with structure as defined in CRIS-DS-ISR-GS-0007-Annex-A_v5a_PFS_L1
%
% ----------------------------------------------------------
% Author:    Juan Pedro López-Zaragoza / isardSAT
% Reviewer:  
% Last rev.: 
% 
% Versions
% 1.0 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function [files] = create_NetCDF_L1Bs_HR(files, N_surfs_loc_estimated,N_bursts,max_beams_x_stack, cnf, chd, cst, resize)

date_creation = datestr(now, '_yyyymmddTHHMMSS_');

% (TO DO!!!) WRITE THE CORRECT NAMING OF THE OUTPUT PRODUCT
L1A_date_find = strfind(files.filename_L1A,'_201');
L1A_date_find = strfind(files.filename_L0,'_202'); %JPLZ: date for IRIS OB chronogram simulations is 202X
L1A_mode_find = strfind(files.filename_L1A, 'PICE_SIRS_');

% if (~resize)
%     files.filename_L1Bs = strcat(files.outputPath,'PICE_SIRS_',...
%         files.filename_L1A((L1A_mode_find(1)+9):(L1A_mode_find(end)+14)), '_1Bs_',...
%         files.filename_L1A((L1A_date_find(1)+1):(L1A_date_find(1)+30)),...
%         '_isd','_long','.nc');
% else
%     files.filename_L1Bs = strcat(files.outputPath,'PICE_SIRS_',...
%         files.filename_L1A((L1A_mode_find(1)+9):(L1A_mode_find(end)+14)), '_1Bs_',...
%         files.filename_L1A((L1A_date_find(1)+1):(L1A_date_find(1)+30)),...
%         '_isd','.nc');
% end

% (TO DO!) Review correct file naming prefixes/suffixes for all configurations
if (~resize)
    files.filename_L1Bs = strcat(files.outputPath,'CRA_IR_1S_HR_SIO_',... % add file prefixes
        files.filename_L1A((L1A_date_find(1)+1):(L1A_date_find(1)+31)),... % add data start and end Rx time
        date_creation, ...% add processing time
        '_____________ISRD_SIR____TST','_long','.nc'); %add final suffixes
else
    files.filename_L1Bs = strcat(files.outputPath,'CRA_IR_1S_HR_SIO_',... % add file prefixes
        files.filename_L1A((L1A_date_find(1)+1):(L1A_date_find(1)+31)),... % add data start and end Rx time
        date_creation, ...% add processing time
        '_____________ISRD_SIR____TST','.nc'); % add final suffixes
end


% files.filename_L1Bs = strcat(files.outputPath,'SR_1_SRA_BS_',...
%                                 files.filename_L1A(end-37:end-3),...
%                                 '.nc');

% Choose NETCDF writing mode: netcdf4 format and overwriting  
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLOBBER'));
ncid = netcdf.create(files.filename_L1Bs,cmode);

long_name_att = 'long_name';
std_name_att = 'standard_name';
calendar_name_att='calendar';
comment_att = 'comment';
units_att = 'units';
scale_factor_att = 'scale_factor';
add_offset_att = 'add_offset';
flag_values_att='flag_values';
flag_desc_att='flag_meanings';
dimensions_key = 'Dimensions';
format_key = 'Format';
data_type_key = 'DataType';
fill_value_key='FillValue';

netcdf_v4_format = 'netcdf4';

% nr_dimension = 'nr';
% nl_dimension = 'Nl';
% nl_dimension_size = N_max_beams_stack;
% ns_dimension = 'Ns';
% space_3D_dimension = 'space_3D';
% space_3D_dimension_size = 3;


n0_dimension = netcdf.defDim(ncid,'n0',1);
nr_dimension = netcdf.defDim(ncid,'nr',N_surfs_loc_estimated);
np_dimension = netcdf.defDim(ncid,'np',chd.N_pulses_burst); % Number of pulses
nl_dimension = netcdf.defDim(ncid,'nl',max_beams_x_stack);
ns_dimension = netcdf.defDim(ncid,'ns',chd.N_samples_sar*cnf.zp_fact_range);
space_3D_dimension = netcdf.defDim(ncid,'space_3D',3);

day_units = 'day';
seconds_units = 'seconds';
seconds_3dot125d64d1e9_units='3.125/64*1e-9 seconds';
seconds_3dot125d1024d1e9_units='3.125/1024*1e-9 seconds';
number_units = 'count';
degrees_units = 'degrees';
meters_units = 'meters';
meters_per_second_units = 'm/s';
rate_per_second_units='1/s';
dB_units = 'dB';
fft_pow_units='FFT power unit';
Hz_units = 'Hz';
T0d64_units = 'T0/64';
T0d16d64_units = 'T0/16/64';
W_per_count_units = 'Watt/#';
sqrtW_per_count_units = 'sqrt(Watt)/#';
rad_units = 'rad';
percent_units='percent';
Vd1_units= 'V/1';
celsius_units= 'celsius';

int8_type = 'NC_BYTE';
uint8_type = 'NC_BYTE';
int16_type = 'NC_SHORT';
uint16_type = 'NC_SHORT';
int32_type = 'NC_INT';
uint32_type = 'NC_INT';
uint64_type = 'uint64';
float_type = 'NC_FLOAT';
double_type= 'NC_DOUBLE';

%% CODING L1B

%--------- Global Attribtues of the netCDF ---------------------------
%                 altimeter_sensor_name='PICE RA'; 
        %gnss sensor name
%         if isempty(files.sph.gnss_info) 
%             gnss_sensor_name='Not available';
%         else 
%             gnss_sensor_name=files.sph.gnss_info;
%         end
%         %doris sensor name
%         if isempty(files.sph.gnss_info) 
%             doris_sensor_name='Not available';
%         else 
%             doris_sensor_name=files.sph.gnss_info;
%         end
        %Acquisition station
%         acq_station_name=files.sph.acq_station;
%         % UTC date of the first measurement
%         first_meas_time=files.sph.product_info.product_time_info.produc_start_time;
%         % UTC date of the last measurement
%         last_meas_time=files.sph.product_info.product_time_info.produc_stop_time;
%         % Name of the altimeter level 0 data file
%         xref_altimeter_level0=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'SIRAL_LEVEL_0',length('SIRAL_LEVEL_0'))).filename;
%         % Name of the file containing the DORIS-derived USO frequency
%         xref_altimeter_orbit=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'ORBIT',length('ORBIT'))).filename;
%         % Name of the file containing the DORIS-derived USO frequency
%         xref_doris_USO=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'DORIS',length('DORIS'))).filename;
%         % Name of the LTM  file containing the SAR mode CAL1 parameters
%         xref_altimeter_ltm_sar_cal1=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'CALIBRATION_TYPE_1',length('CALIBRATION_TYPE_1'))).filename;
%         % Name of the LTM  file containing the Ku-band CAL2 parameters
%         idx=find(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'CALIBRATION_TYPE_2',length('CALIBRATION_TYPE_2'))~=0);
%         if isempty(idx)
%            xref_altimeter_ltm_ku_cal2='Not available (not applied in the product)';
%         else
%            xref_altimeter_ltm_ku_cal2=files.sph.dsds(idx).filename;
%         end
%         % Name of the LTM file containing the C-band CAL2 parameters         
%         xref_altimeter_ltm_c_cal2='Not available for CR2';
%         % Name of the altimeter characterisation data file
%         xref_altimeter_characterisation=files.sph.dsds(strncmp(strsplit(strtrim([files.sph.dsds.ds_name])),'IPF_RA_DATABASE',length('IPF_RA_DATABASE'))).filename;
        % Semi-major axis of the reference ellipsoid
%         semi_major_ellipsoid_axis=num2str(cst.semi_major_axis,15);
        % Flattening coeffcient of the reference ellipsoid
%         ellipsoid_flattening=num2str(cst.flat_coeff,15);    
        % Attributes related to the OrbitalN_surfs_loc_estimatedormation required by Porto
%         orbit_phase_code=files.sph.orbit_info.phase_code;
%         orbit_cycle_num=((files.sph.orbit_info.cycle_num));
%         orbit_REL_Orbit=((files.sph.orbit_info.rel_orbit));
%         orbit_ABS_Orbit_Start=((files.sph.orbit_info.ABS_Orbit_Start));
%         orbit_Rel_Time_ASC_Node_Start=((files.sph.orbit_info.Rel_Time_ASC_Node_Start));
%         orbit_ABS_Orbit_Stop=((files.sph.orbit_info.ABS_Orbit_Stop));
%         orbit_Rel_Time_ASC_Node_Stop=((files.sph.orbit_info.Rel_Time_ASC_Node_Stop));
%         orbit_Equator_Cross_Time=(files.sph.orbit_info.Equator_Cross_Time);
%         orbit_Equator_Cross_Long=((files.sph.orbit_info.Equator_Cross_Long));
% 

%%----------0. Define groups ----------------------------------------------
global_ncid = netcdf.defGrp(ncid,'global');
global_ku_ncid = netcdf.defGrp(global_ncid,'ku');
global_ka_ncid = netcdf.defGrp(global_ncid,'ka');

data_ncid = netcdf.defGrp(ncid,'data');
data_ku_ncid = netcdf.defGrp(data_ncid,'ku');
data_ka_ncid = netcdf.defGrp(data_ncid,'ka');

if strcmp(chd.band,'Ku')
    global_currentband_ncid=global_ku_ncid;
    data_currentband_ncid=data_ku_ncid;
elseif strcmp(chd.band,'Ka')
    global_currentband_ncid=global_ka_ncid;
    data_currentband_ncid=data_ka_ncid;
end

%% PACKING L1Bs

%----------A. Time & counters variables ----------------------------------------------
looks_name = 'looks';
id_aux      = netcdf.defVar(ncid,looks_name,int16_type,nl_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Index of looks within a stack');
netcdf.putAtt(ncid,id_aux,comment_att,'Counter to index the looks within a stack.');

samples_ov_name = 'samples_ov';
id_aux      = netcdf.defVar(ncid,samples_ov_name,int16_type,ns_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Waveform sample number (with oversampling).');
netcdf.putAtt(ncid,id_aux,comment_att,'Counter to index the waveform samples after oversampling.');

space3d_name = 'space_3d';
id_aux      = netcdf.defVar(ncid,space3d_name,uint8_type,space_3D_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Index of carteesian coordinates.');
netcdf.putAtt(ncid,id_aux,comment_att,'Counter to index the cartesian coordinates.');

pulses_name = 'pulses';
id_aux      = netcdf.defVar(ncid,pulses_name,uint8_type,np_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Index of pulses within a burst.');
netcdf.putAtt(ncid,id_aux,comment_att,'Counter to index the pulses within a burst.');

l1b_record_counter_name = 'l1b_record_counter';
id_aux      = netcdf.defVar(data_currentband_ncid,l1b_record_counter_name,int32_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'l1b_record_counter');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'L1B record counter');
netcdf.putAtt(data_currentband_ncid,id_aux,calendar_name_att,'Gregorian');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'The L1B record counter, starting from 0 for each product.');

time_name = 'time';
id_aux      = netcdf.defVar(data_currentband_ncid,time_name,double_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'time');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'UTC seconds since 2000-01-01 00:00:00.0+00:00');
netcdf.putAtt(data_currentband_ncid,id_aux,calendar_name_att,'Gregorian');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Time refers to the instant the L1B waveform touches the surface. [tai_utc_difference] is the difference between TAI and UTC reference time (seconds) for the first measurement of the data set. [leap_second] is the UTC time at which a leap second occurs in the data set, if any. After this UTC time, the [tai_utc_difference] is increased by 1 second.');

time_tai_name = 'time_tai';
id_aux      = netcdf.defVar(data_currentband_ncid,time_tai_name,double_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'time');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'TAI seconds since 2000-01-01 00:00:00.0+00:00');
netcdf.putAtt(data_currentband_ncid,id_aux,calendar_name_att,'Gregorian');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Time refers to the instant the L1B waveform touches the surface. [tai_utc_difference] is the difference between TAI and UTC reference time (seconds) for the first measurement of the data set. [leap_second] is the UTC time at which a leap second occurs in the data set, if any. After this UTC time, the [tai_utc_difference] is increased by 1 second.');

tm_source_sequence_counter_name = 'tm_source_sequence_counter';
id_aux      = netcdf.defVar(data_currentband_ncid,tm_source_sequence_counter_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'tm_source_sequence_counter');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Instrument source sequence counter');
netcdf.putAtt(data_currentband_ncid,id_aux,calendar_name_att,'Gregorian');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'This is the instrument record counter, copied from ISP.');

%----------B. Orbit and attitude variables ------------------------------
altitude_name = 'altitude';
id_aux = netcdf.defVar(data_currentband_ncid,altitude_name,int32_type,nr_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'altitude');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Center of mass altitude above the reference ellipsoid');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_currentband_ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Altitude of the satellite center of mass above the reference ellipsoid (WGS-84).');


altitude_rate_name = 'altitude_rate';
id_aux = netcdf.defVar(data_currentband_ncid,altitude_rate_name,int32_type,nr_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'altitude rate');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Center of mass altitude rate with respect to the reference ellipsoid');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-4);
% netcdf.putAtt(data_currentband_ncid,id_aux,add_offset_att,0);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Instantaneous altitude rate at the center of mass with respect to the reference ellipsoid (WGS-84).');

latitude_name = 'latitude';
id_aux = netcdf.defVar(data_currentband_ncid,latitude_name,int32_type,nr_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
% netcdf.putAtt(data_currentband_ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'latitude');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'latitude (positive N, negative S)');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Latitude of measurement [-90, +90]. Positive latitude is North latitude, negative latitude is South latitude.');

longitude_name = 'longitude';
id_aux = netcdf.defVar(data_currentband_ncid,longitude_name,int32_type,nr_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'longitude');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'longitude (positive towards E)');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-6);
% netcdf.putAtt(data_currentband_ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Longitude of measurement [0, 360). East longitude relative to Greenwich meridian.');


if strcmp(chd.band,'Ku')
    off_nadir_roll_angle_ant1_name = 'off_nadir_roll_angle_ant1';
    id_aux = netcdf.defVar(data_ku_ncid,off_nadir_roll_angle_ant1_name,int32_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'off_nadir_roll_angle_ant1');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Off-nadir roll angle measured by the antenna 1');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,degrees_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-7);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'Roll angle of antenna 1 with respect to the nadir pointing, measured by the STRs and post-processed by Attitude and Orbit Control System (AOCS) or by ground facility. The attribute add_offset contains the estimated roll bias. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: roll = 0.');
    
    off_nadir_pitch_angle_ant1_name = 'off_nadir_pitch_angle_ant1';
    id_aux = netcdf.defVar(data_ku_ncid,off_nadir_pitch_angle_ant1_name,int32_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'off_nadir_pitch_angle_ant1');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Off-nadir pitch angle measured by the antenna 1');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,degrees_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-7);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'Pitch angle of antenna 1 with respect to the nadir pointing, measured by the STRs and post-processed by AOCS or by ground facility. The attribute add_offset contains the estimated pitch bias. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: pitch = 0.');
    
    off_nadir_yaw_angle_ant1_name = 'off_nadir_yaw_angle_ant1';
    id_aux = netcdf.defVar(data_ku_ncid,off_nadir_yaw_angle_ant1_name,int32_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'off_nadir_yaw_angle_ant1');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Off-nadir yaw angle measured by the antenna 1');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,degrees_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-7);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'Yaw angle of antenna 1 with respect to the nadir pointing, measured by the STRs and post-processed by AOCS or by ground facility. The attribute add_offset contains the estimated yaw bias. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: yaw = 0.');
end

if strcmp(cnf.processing_mode,'SIN')
    off_nadir_roll_angle_ant2_name = 'off_nadir_roll_angle_ant2';
    id_aux = netcdf.defVar(data_ku_ncid,off_nadir_roll_angle_ant2_name,int32_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'off_nadir_roll_angle_ant2');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Off-nadir roll angle measured by the antenna 2');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,degrees_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-7);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'Roll angle of antenna 2 with respect to the nadir pointing, measured by the STRs and post-processed by Attitude and Orbit Control System (AOCS) or by ground facility. The attribute add_offset contains the estimated roll bias. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: roll = 0.');
    
    off_nadir_pitch_angle_ant2_name = 'off_nadir_pitch_angle_ant2';
    id_aux = netcdf.defVar(data_ku_ncid,off_nadir_pitch_angle_ant2_name,int32_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'off_nadir_pitch_angle_ant2');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Off-nadir pitch angle measured by the antenna 2');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,degrees_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-7);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'Pitch angle of antenna 2 with respect to the nadir pointing, measured by the STRs and post-processed by AOCS or by ground facility. The attribute add_offset contains the estimated pitch bias. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: pitch = 0.');
    
    off_nadir_yaw_angle_ant2_name = 'off_nadir_yaw_angle_ant2';
    id_aux = netcdf.defVar(data_ku_ncid,off_nadir_yaw_angle_ant2_name,int32_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'off_nadir_yaw_angle_ant2');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Off-nadir yaw angle measured by the antenna 2');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,degrees_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-7);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'Yaw angle of antenna 2 with respect to the nadir pointing, measured by the STRs and post-processed by AOCS or by ground facility. The attribute add_offset contains the estimated yaw bias. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: yaw = 0.');
end

if strcmp(chd.band,'Ka')
    off_nadir_roll_angle_name = 'off_nadir_roll_angle';
    id_aux = netcdf.defVar(data_ka_ncid,off_nadir_roll_angle_name,int32_type,nr_dimension);
    netcdf.putAtt(data_ka_ncid,id_aux,std_name_att,'off_nadir_roll_angle');
    netcdf.putAtt(data_ka_ncid,id_aux,long_name_att,'Off-nadir roll angle measured by the Ka-band Tx-Rx antenna');
    netcdf.putAtt(data_ka_ncid,id_aux,units_att,degrees_units);
    netcdf.putAtt(data_ka_ncid,id_aux,scale_factor_att,1.e-7);
    netcdf.putAtt(data_ka_ncid,id_aux,comment_att,'Roll angle Ka-band Tx-Rx antenna with respect to the nadir pointing, measured by the STRs and post-processed by Attitude and Orbit Control System (AOCS) or by ground facility. The attribute add_offset contains the estimated roll bias. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: roll = 0.');
    
    off_nadir_pitch_angle_name = 'off_nadir_pitch_angle';
    id_aux = netcdf.defVar(data_ka_ncid,off_nadir_pitch_angle_name,int32_type,nr_dimension);
    netcdf.putAtt(data_ka_ncid,id_aux,std_name_att,'off_nadir_pitch_angle');
    netcdf.putAtt(data_ka_ncid,id_aux,long_name_att,'Off-nadir pitch angle measured by the Ka-band Tx-Rx antenna');
    netcdf.putAtt(data_ka_ncid,id_aux,units_att,degrees_units);
    netcdf.putAtt(data_ka_ncid,id_aux,scale_factor_att,1.e-7);
    netcdf.putAtt(data_ka_ncid,id_aux,comment_att,'Pitch angle of Ka-band Tx-Rx antenna with respect to the nadir pointing, measured by the STRs and post-processed by AOCS or by ground facility. The attribute add_offset contains the estimated pitch bias. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: pitch = 0.');
    
    off_nadir_yaw_angle_name = 'off_nadir_yaw_angle';
    id_aux = netcdf.defVar(data_ka_ncid,off_nadir_yaw_angle_name,int32_type,nr_dimension);
    netcdf.putAtt(data_ka_ncid,id_aux,std_name_att,'off_nadir_yaw_angle');
    netcdf.putAtt(data_ka_ncid,id_aux,long_name_att,'Off-nadir yaw angle measured by the Ka-band Tx-Rx antenna');
    netcdf.putAtt(data_ka_ncid,id_aux,units_att,degrees_units);
    netcdf.putAtt(data_ka_ncid,id_aux,scale_factor_att,1.e-7);
    netcdf.putAtt(data_ka_ncid,id_aux,comment_att,'Yaw angle of Ka-band Tx-Rx antenna with respect to the nadir pointing, measured by the STRs and post-processed by AOCS or by ground facility. The attribute add_offset contains the estimated yaw bias. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: yaw = 0.');
end

orbit_type_flag_name = 'orbit_type_flag';
id_aux = netcdf.defVar(global_currentband_ncid,orbit_type_flag_name,int8_type,n0_dimension);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'orbit_type');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Orbit type used');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'Orbit type used for the product generation.');

position_vector_name = 'position_vector';
id_aux = netcdf.defVar(data_currentband_ncid,position_vector_name,double_type,[space_3D_dimension ,nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'position_vector');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Center of mass position vector in ITRF, components: [1] x, [2] y, [3] z');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Position vector at the center of mass in ITRF.  The 3 components are given according to the space_3d dimension: [1] x, [2] y, [3] z.');

velocity_vector_name = 'velocity_vector';
id_aux = netcdf.defVar(data_currentband_ncid,velocity_vector_name,int32_type,[space_3D_dimension ,nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'velocity_vector');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Center of mass velocity vector in ITRF, components: [1] x, [2] y, [3] z');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Velocity vector at the center of mass in ITRF.  The 3 components are given according to the space_3d dimension: [1] x, [2] y, [3] z.');

%----------C. Configuration and quality variables --------------------------------------
mcd_flags_name = 'mcd_flags';
id_aux = netcdf.defVar(data_currentband_ncid,mcd_flags_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'mcd_flags');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Measurement confidence data (MCD): collection of error flags.');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Set of error flags linked to the information retrieval common for the different products. When equal to 0 the information has been retrieved correctly. If equal to 1 an error has occurred when retrieving the information.');

processing_configuration_flags_name = 'processing_configuration_flags';
id_aux = netcdf.defVar(global_currentband_ncid,processing_configuration_flags_name,int8_type,n0_dimension);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'processing_configuration_flags');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Collection of parameters from the processor configuration file.');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'uso_correction_activated (0 - deactivated, 1 - activated), cal1_correction_activated (0 - deactivated, 1 - activated), cal1_power_correction_max (0 - total power, 1 - maximum power), cal1_delay_correction_max (0 - cog, 1 - max power).');

hr_mcd_flags_name = 'hr_mcd_flags';
id_aux = netcdf.defVar(data_currentband_ncid,hr_mcd_flags_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'hr_mcd_flags');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Measurement confidence data (MCD): collection of error flags for the HR processor.');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Set of error flags linked to the information retrieval common for the different products. When equal to 0 the information has been retrieved correctly. If equal to 1 an error has occurred when retrieving the information.');

hr_processing_configuration_flags_name = 'hr_processing_configuration_flags';
id_aux = netcdf.defVar(global_currentband_ncid,hr_processing_configuration_flags_name,int8_type,n0_dimension);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'hr_processing_configuration_flags');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Collection of parameters for the HR chain from the processor configuration file.');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'emove_doppler_ambiguities (0 - not removed, 1 - removed), azimuth_weighing_applied (0 - not applied, 1 - applied), antenna_weighting_applied (0 - not applied, 1 - applied), surface_weighting_applied (0 - not applied, 1 - applied), avoid_zeros_in_multilooking (0 - not avoided, 1 - avoided), exact_azimuth_processing (0 - approximate method, 1 - exact method).');

iris_instrument_configuration_flags_name = 'iris_instrument_configuration_flags';
id_aux = netcdf.defVar(data_currentband_ncid,iris_instrument_configuration_flags_name,int8_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'iris_instrument_configuration_flags');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Collection of flags from the IRIS instrument source packets');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'It contains flags from ISPs related to IRIS configuration and tracking.');

iris_mode_flag_name = 'iris_mode_flag';
id_aux = netcdf.defVar(data_currentband_ncid,iris_mode_flag_name,int8_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'iris_mode_flag');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Mode identifier from the IRIS instrument source packets');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Flag containing the mode identifier from the IRIS instrument source packets.');

range_oversampling_factor_name = 'range_oversampling_factor';
id_aux = netcdf.defVar(global_currentband_ncid,range_oversampling_factor_name,int8_type,n0_dimension);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'range_oversampling_factor');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Oversampling factor used in the range compression (FFT)');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'The instrument samples the waveforms with a 395 MHz clock, providing a nominal_sampling = c / 395e6 / 2 = ~0.379m (with c=speed of light). In addition, the ground processor can apply an oversampling factor, providing a waveform_sampling = nominal_sampling / range_oversampling_factor. Note that the altimeter range resolution is fixed and given by the chirp bandwidth of 320 MHz: c / 320e6 / 2 = ~0.468m.');

telemetry_type_flag_name = 'telemetry_type_flag';
id_aux = netcdf.defVar(data_currentband_ncid,telemetry_type_flag_name,int8_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'telemetry_type_flag');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Telemetry type: RAW or RMC');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Telemetry packet type can be: RAW or RMC.');

%----------------D. Altimeter range variables ---------------------------
range_cor_com_ant1_name = 'range_cor_com_ant1';
id_aux = netcdf.defVar(data_currentband_ncid,range_cor_com_ant1_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'range_cor_com_ant1');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Range correction: 1-way distance CoM - antenna 1');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'It is the (1-way) distance in the z-component between the center of mass and the antenna 1 reference point. The attitude at the SRF is used for the computation of this correction. The z-component closest in time is used.');

if strcmp(cnf.processing_mode,'SIN')
    range_cor_com_ant2_name = 'range_cor_com_ant2';
    id_aux = netcdf.defVar(data_ku_ncid,range_cor_com_ant2_name,int16_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'range_cor_com_ant2');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Range correction: 1-way distance CoM - antenna 2');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'It is the (1-way) distance in the z-component between the center of mass and the antenna 2 reference point. The attitude at the SRF is used for the computation of this correction. The z-component closest in time is used.');
end

if strcmp(chd.band,'Ka')
    range_cor_com_name = 'range_cor_com';
    id_aux = netcdf.defVar(data_ka_ncid,range_cor_com_name,int16_type,nr_dimension);
    netcdf.putAtt(data_ka_ncid,id_aux,std_name_att,'range_cor_com');
    netcdf.putAtt(data_ka_ncid,id_aux,long_name_att,'range correction: 1-way distance CoM - Ka-band Tx-Rx antenna');
    netcdf.putAtt(data_ka_ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(data_ka_ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(data_ka_ncid,id_aux,comment_att,'It is the (1-way) distance in the z-component between the center of mass and the antenna 2 reference point. The attitude at the SRF is used for the computation of this correction. The z-component closest in time is used.');
end

range_cor_doppler_name = 'range_cor_doppler';
id_aux = netcdf.defVar(data_currentband_ncid,range_cor_doppler_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'range_cor_doppler');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Range correction: 1-way Doppler compensation');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'One-way Doppler range correction due to satellite altitude rate (only applied in LR and LR-OS).');

range_cor_internal_delay_cal_rx1_name = 'range_cor_internal_delay_cal_rx1';
id_aux = netcdf.defVar(data_currentband_ncid,range_cor_internal_delay_cal_rx1_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'range_cor_internal_delay_cal_rx1');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Range correction: 1-way internal calibration path delay (from CAL1), for antenna Rx1.');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'The internal path delay, for antenna Rx1, is periodically measured in flight with a dedicated CAL pulse per radar cycle, or by means CAL1 LR and HR sequences. The CAL1 closest in time is used.');

if strcmp(cnf.processing_mode,'SIN')
    range_cor_internal_delay_cal_rx2_name = 'range_cor_internal_delay_cal_rx2';
    id_aux = netcdf.defVar(data_ku_ncid,range_cor_internal_delay_cal_rx2_name,int16_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'range_cor_internal_delay_cal_rx2');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Range correction: 1-way internal calibration path delay (from CAL1), for antenna Rx2.');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'The internal path delay, for antenna Rx2, is periodically measured in flight with a dedicated CAL pulse per radar cycle, or by means CAL1 LR and HR sequences. The CAL1 closest in time is used.');
end

range_cor_internal_delay_att_rx1_name = 'range_cor_internal_delay_att_rx1';
id_aux = netcdf.defVar(data_currentband_ncid,range_cor_internal_delay_att_rx1_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'range_cor_internal_delay_att_rx1');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'range correction: 1-way range compensation for the delay introduced by the attenuator, for antenna Rx1');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'This is the 1-way range correction due to the delay introduced by the attenuator, for antenna Rx1.');

if strcmp(cnf.processing_mode,'SIN')
    range_cor_internal_delay_att_rx2_name = 'range_cor_internal_delay_att_rx2';
    id_aux = netcdf.defVar(data_ku_ncid,range_cor_internal_delay_att_rx2_name,int16_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'range_cor_internal_delay_att_rx2');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Range correction: 1-way range compensation for the delay introduced by the attenuator, for antenna Rx2.');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'This is the 1-way range correction due to the delay introduced by the attenuator, for antenna Rx2.');
end

range_cor_uso_name = 'range_cor_uso';
id_aux = netcdf.defVar(data_currentband_ncid,range_cor_uso_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'range_cor_uso');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Range correction: 1-way USO drift.');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-5);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'It is the USO Drift (1-way) range correction. The USO drift is applied directly to the clock. This variable has been included for monitoring purpose. The USO Drift closest in time is used.');


range_cor_reference_sample_name = 'range_cor_reference_sample';
id_aux = netcdf.defVar(data_currentband_ncid,range_cor_reference_sample_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'range_cor_reference_sample');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Range correction: 1-way compensation to shift the tracker range to the reference sample.');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'IRIS includes a 32768 samples IFFT on board. After the IFF there is a range selection block to allow for selecting a subset of these 32768 samples. This subset is then telemetered on ground. The first sample, and the width of the selected echo is configurable per measurement mode and available in telemetry. [N_FS-N_CS+samples]*T0 is the two-way correction, which multiplied by the speed of light and converted to one-way is equivalent to this variable.');

tracker_range_calibrated_name = 'tracker_range_calibrated';
id_aux = netcdf.defVar(data_currentband_ncid,tracker_range_calibrated_name,double_type,nr_dimension); %JPLZ: Official var type is int32, but for some reason MATLAB is not saving correctly the decimals. With double var type it works
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'tracker_range_calibrated');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Calibrated 1-way tracker range: CoM to middle range window (at sample ns/2 from 0).');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_currentband_ncid,id_aux,add_offset_att,7e5);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'This is the 1-way distance from the satellite center of mass to the middle of the range window (sample ns/2 from 0). It includes the following range calibrations: (a) range_corr_internal_delay_cal, (b) range_corr_internal_delay_diff, (c) range_corr_external_group_delay, (d) range_corr_com. Note: the actual altimeter clock (variable altimeter_clock) has been used to compute the altimeter range.');

tracker_range_calibrated_gnss_name = 'tracker_range_calibrated_gnss';
id_aux = netcdf.defVar(data_currentband_ncid,tracker_range_calibrated_gnss_name,int32_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'tracker_range_calibrated_gnss');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Calibrated 1-way tracker range given by GNSS+DEM: CoM to middle range window (at sample ns/2 from 0).');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(data_currentband_ncid,id_aux,add_offset_att,7e5);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'This is the 1-way distance from the satellite center of mass to the middle of the range window as given by the GNSS+DEM (sample ns/2 from 0). It includes the following range calibrations: (a) range_corr_internal_delay_cal, (b) range_corr_internal_delay_diff, (c) range_corr_external_group_delay, (d) range_corr_com. Note: the actual altimeter clock (variable altimeter_clock) has been used to compute the altimeter range.');

%----------E. Altimeter power variables ------------------------------
altimeter_power_drift_rx1_name = 'altimeter_power_drift_rx1';
id_aux = netcdf.defVar(data_currentband_ncid,altimeter_power_drift_rx1_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'altimeter_power_drift_rx1');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Power drift due to orbital and aging variations (from CAL1), for antenna Rx1.');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'This is a measure of the instrument power drift. It is mainly composed of a temperature-dependent harmonic variation along the orbit and a long-term decay. It is regularly measured from CAL1 impulse response with respect to a reference at the beginning of the mission. the instrument power drift closest in time is used. For antenna Rx1.');

if strcmp(cnf.processing_mode,'SIN')
    altimeter_power_drift_rx2_name = 'altimeter_power_drift_rx2';
    id_aux = netcdf.defVar(data_ku_ncid,altimeter_power_drift_rx2_name,int16_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'altimeter_power_drift_rx1');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Power drift due to orbital and aging variations (from CAL1), for antenna Rx2.');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'This is a measure of the instrument power drift. It is mainly composed of a temperature-dependent harmonic variation along the orbit and a long-term decay. It is regularly measured from CAL1 impulse response with respect to a reference at the beginning of the mission. the instrument power drift closest in time is used. For antenna Rx2.');
end

attenuation_calibrated_rx1_name = 'attenuation_calibrated_rx1';
id_aux = netcdf.defVar(data_currentband_ncid,attenuation_calibrated_rx1_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'attenuation_calibrated_rx1');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Power scaling: ATT calibrated (from ATT table), for antenna Rx1.');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'The ISPs contain the ATT command applied within a radar cycle. The actual ATT value, which slightly differs from the nominal one, is extracted from a look-up table. The Attenuation closest in time is used. For antenna Rx1.');


if strcmp(cnf.processing_mode,'SIN')
    attenuation_calibrated_rx2_name = 'attenuation_calibrated_rx2';
    id_aux = netcdf.defVar(data_ku_ncid,attenuation_calibrated_rx2_name,int16_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'attenuation_calibrated_rx2');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Power scaling: ATT calibrated (from ATT table), for antenna Rx2.');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'The ISPs contain the ATT command applied within a radar cycle. The actual ATT value, which slightly differs from the nominal one, is extracted from a look-up table. The Attenuation closest in time is used. For antenna Rx2.');
end

power_scaling_to_antenna_rx1_name = 'power_scaling_to_antenna_rx1';
id_aux = netcdf.defVar(data_currentband_ncid,power_scaling_to_antenna_rx1_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'power_scaling_to_antenna_rx1');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Power scaling: overall instrument gain, for antenna Rx1.');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'The overall instrument gain, for antenna Rx1, including: (a) attenuator_calibrated, (b) variable_digital_gain, (c) cal1_power_drift, (d) g_scaling, (e) antenna_gain. This variable is applied to the power waveform, and once it is applied the resulting L1B waveform is equivalent to the ratio of received vs transmitted power.');

if strcmp(cnf.processing_mode,'SIN')
    power_scaling_to_antenna_rx2_name = 'power_scaling_to_antenna_rx2';
    id_aux = netcdf.defVar(data_ku_ncid,power_scaling_to_antenna_rx2_name,int16_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'power_scaling_to_antenna_rx2');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Power scaling: overall instrument gain, for antenna Rx2.');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'The overall instrument gain, for antenna Rx2, including: (a) attenuator_calibrated, (b) variable_digital_gain, (c) cal1_power_drift, (d) g_scaling, (e) antenna_gain. This variable is applied to the power waveform, and once it is applied the resulting L1B waveform is equivalent to the ratio of received vs transmitted power.');
end

variable_digital_gain_name = 'variable_digital_gain';
id_aux = netcdf.defVar(data_currentband_ncid,variable_digital_gain_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'variable_digital_gain');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Power scaling: variable digital gain');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Variable component of the gain produced by the digital processing unit. The variable gain closest in time is used.');

cal1_power_rx1_name = 'cal1_power_rx1';
id_aux = netcdf.defVar(data_currentband_ncid,cal1_power_rx1_name,int32_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'cal1_power_rx1');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Instrument impulse response maximum or integrated power, for antenna Rx1');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Max or Total power of the instrument PTR from the CAL1 measurement, for antenna Rx1. The choice of Max vs Total power is specified in the mcd_flags. Either case, the cal1_power waveform has been corrected for (a) calibrated attenuation correction, (b) variable digital gain and (c) residual cal1 fixed digital gain.');

if strcmp(cnf.processing_mode,'SIN')
    cal1_power_rx2_name = 'cal1_power_rx2';
    id_aux = netcdf.defVar(data_ku_ncid,cal1_power_rx2_name,int32_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'cal1_power_rx2');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Instrument impulse response maximum or integrated power, for antenna Rx2');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'Max or Total power of the instrument PTR from the CAL1 measurement, for antenna Rx2. The choice of Max vs Total power is specified in the mcd_flags. Either case, the cal1_power waveform has been corrected for (a) calibrated attenuation correction, (b) variable digital gain and (c) residual cal1 fixed digital gain.'); 
end

%----------F. Altimeter engineering variables --------------------------
altimeter_clock_name = 'altimeter_clock';
id_aux = netcdf.defVar(data_currentband_ncid,altimeter_clock_name,int32_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'Altimeter_clock');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Altimeter clock');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,Hz_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_currentband_ncid,id_aux,add_offset_att,6e8);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'This is the actual altimeter clock. The altimeter clock is based upon the USO clock. The nominal USO clock is 10MHz, while the nominal altimeter clock is 600 MHz. The actual USO clock is provided regularly as a drift information: based on the USO drift, the actual altimeter clock is computed, assuming a linear dependency. The USO clock closest in time is used.');

hn_mean_name = 'hn_mean';
id_aux = netcdf.defVar(data_currentband_ncid,hn_mean_name,uint32_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'hn_mean');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Mean altitude instruction, for monitoring purposes');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,T0d16d64_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'This field is for monitoring the H0/COR2 decoding into CAI (Coarse Altitude instruction) and FAI (Fine Altitude Instruction). In LRM it is the altitude instruction (CAI+FAI) averaged and rounded over the radar cycle.  In SAR, it is the altitude instruction (CAI+FAI) averaged and rounded over the burst closest to the nadir of the surface location. NOTE: T0 = 1/altimeter_clock. The mean altitude closest in time is used.');

pulse_repetition_interval_name = 'pulse_repetition_interval';
id_aux = netcdf.defVar(data_currentband_ncid,pulse_repetition_interval_name,int32_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'pulse_repetition_interval');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Pulse Repetition Interval converted into seconds');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-12);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'The Pulse Repetition Interval. PRI is constant within all received pulses in a radar cycle, but it can change within consecutive radar cycles. It is provided in counters of (T0*8) by the altimeter, with T0 = 1/altimeter_clock. The PRI closest in time is used.');

tm_cor2_name = 'tm_cor2';
id_aux = netcdf.defVar(data_currentband_ncid,tm_cor2_name,uint32_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'tm_cor2');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'COR2 telemetry (same in 1 radar cycle)');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,T0d16d64_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'The COR2 (altitude rate estimation) value, copied from ISP. The altimeter estimates one COR2 value per radar cycle. It is used in combination with the H0 to compute the CAI and FAI instructions. NOTE: T0 = 1/altimeter_clock. The COR2 closest in time is used.');

tm_h0_name = 'tm_h0';
id_aux = netcdf.defVar(data_currentband_ncid,tm_h0_name,uint32_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'tm_h0');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'h0 telemetry (same in 1 rada cycle)');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,T0d64_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'The H0 (initial altitude instruction) value, copied from ISP. The altimeter provides one H0 value per radar cycle. It is used in combination with the COR2 to compute the CAI and FAI instructions. NOTE: T0 = 1/altimeter_clock. The closest H0 in time is used.');

%----------G. Altimeter characterization variables --------------------------
residual_fixed_digital_gain_raw_name = 'residual_fixed_digital_gain_raw';
id_aux = netcdf.defVar(global_currentband_ncid,residual_fixed_digital_gain_raw_name,int32_type,nr_dimension);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'residual_fixed_digital_gain_raw');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Power scaling: residual on-board fixed digital processing gain for RAW');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(global_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'This variable accounts for any potential residual gain to be applied to the HR RAW waveform.');

residual_fixed_digital_gain_rmc_name = 'residual_fixed_digital_gain_rmc';
id_aux = netcdf.defVar(global_currentband_ncid,residual_fixed_digital_gain_rmc_name,int32_type,nr_dimension);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'residual_fixed_digital_gain_rmc');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Power scaling: residual on-board fixed digital processing gain for RMC');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(global_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'This variable accounts for any potential residual gain to be applied to the HR RMC waveform.');

range_cor_external_group_delay_rx1_name = 'range_cor_external_group_delay_rx1';
id_aux = netcdf.defVar(global_currentband_ncid,range_cor_external_group_delay_rx1_name,int32_type,nr_dimension);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'range_cor_external_group_delay_rx1');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Range correction: 1-way external group delay from ground characterization, for antenna Rx1');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(global_currentband_ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'The absolute 1-way group delay contribution for parts outside the calibration loop, for antenna Rx1.');

range_cor_external_group_delay_rx2_name = 'range_cor_external_group_delay_rx2';
id_aux = netcdf.defVar(global_currentband_ncid,range_cor_external_group_delay_rx2_name,int32_type,nr_dimension);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'range_cor_external_group_delay_rx2_name');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Range correction: 1-way external group delay from ground characterization, for antenna Rx2');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(global_currentband_ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'The absolute 1-way group delay contribution for parts outside the calibration loop, for antenna Rx2.');

g_scaling_rx1_name = 'g_scaling_rx1';
id_aux = netcdf.defVar(global_currentband_ncid,g_scaling_rx1_name,int32_type,nr_dimension);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'g_scaling_rx1');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Ratio of losses between the RA transmission/reception path and the calibration path, for antenna Rx1');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(global_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'Ratio of losses between the RA transmission/reception and calibration paths from on-ground characterization, for antenna Rx1.');

if strcmp(cnf.processing_mode,'SIN')
    g_scaling_rx2_name = 'g_scaling_rx2';
    id_aux = netcdf.defVar(global_ku_ncid,g_scaling_rx2_name,int32_type,nr_dimension);
    netcdf.putAtt(global_ku_ncid,id_aux,std_name_att,'g_scaling_rx2');
    netcdf.putAtt(global_ku_ncid,id_aux,long_name_att,'Ratio of losses between the RA transmission/reception path and the calibration path, for antenna Rx2');
    netcdf.putAtt(global_ku_ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(global_ku_ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(global_ku_ncid,id_aux,comment_att,'Ratio of losses between the RA transmission/reception and calibration paths from on-ground characterization, for antenna Rx2.');
end

antenna_gain_ant1_name = 'antenna_gain_ant1';
id_aux = netcdf.defVar(global_currentband_ncid,antenna_gain_ant1_name,int32_type,nr_dimension);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'antenna_gain_ant1');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Maximum antenna 1 gain');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(global_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'Maximum antenna 1 gain from on-ground characterization.');

if strcmp(cnf.processing_mode,'SIN')
    antenna_gain_ant2_name = 'antenna_gain_ant2';
    id_aux = netcdf.defVar(global_ku_ncid,antenna_gain_ant2_name,int32_type,nr_dimension);
    netcdf.putAtt(global_ku_ncid,id_aux,std_name_att,'antenna_gain_ant2');
    netcdf.putAtt(global_ku_ncid,id_aux,long_name_att,'Maximum antenna 2 gain');
    netcdf.putAtt(global_ku_ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(global_ku_ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(global_ku_ncid,id_aux,comment_att,'Maximum antenna 2 gain from on-ground characterization.');
end

if strcmp(chd.band,'Ka')
    antenna_gain_ant_name = 'antenna_gain_ant';
    id_aux = netcdf.defVar(global_ka_ncid,antenna_gain_ant_name,int32_type,nr_dimension);
    netcdf.putAtt(global_ka_ncid,id_aux,std_name_att,'antenna_gain_ant');
    netcdf.putAtt(global_ka_ncid,id_aux,long_name_att,'Maximum Ka-band Tx-Rx antenna gain');
    netcdf.putAtt(global_ka_ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(global_ka_ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(global_ka_ncid,id_aux,comment_att,'Maximum Ka-band Tx-Rx antenna gain from on-ground characterization.');
end

%------------------H. Flag variables --------------------------
manoeuvre_flag_name = 'manoeuvre_flag';
id_aux = netcdf.defVar(data_currentband_ncid,manoeuvre_flag_name,uint8_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'manoeuvre_flag');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Flag indicating whether the satellite is in a manoeuvre');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Flags contained within the NAVATT file.');

surface_classification_flag_name = 'surface_classification_flag';
id_aux = netcdf.defVar(data_currentband_ncid,surface_classification_flag_name,uint8_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'surface_classification_flag');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Flag indicating the surface type');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'7-state surface type classification from MODIS and GlobCover.');

%------------------J. Waveform related variables --------------------------
burst_phase_array_cor_rx1_name = 'burst_phase_array_cor_rx1';
id_aux = netcdf.defVar(global_currentband_ncid,burst_phase_array_cor_rx1_name,int32_type,[np_dimension ,n0_dimension]);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'burst_phase_array_cor_rx1');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Burst phase correction, for antenna Rx1');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,rad_units);
netcdf.putAtt(global_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'Burst phase correction applied to the SAR Ku waveforms. It is obtained by averaging the burst phase along a configurable time window to reduce noise, for antenna Rx1. The closest burst phase correction in time is used.');

if strcmp(cnf.processing_mode,'SIN')
    burst_phase_array_cor_rx2_name = 'burst_phase_array_cor_rx2';
    id_aux = netcdf.defVar(global_ku_ncid,burst_phase_array_cor_rx2_name,int32_type,[np_dimension ,n0_dimension]);
    netcdf.putAtt(global_ku_ncid,id_aux,std_name_att,'burst_phase_array_cor_rx2');
    netcdf.putAtt(global_ku_ncid,id_aux,long_name_att,'Burst phase correction, for antenna Rx2');
    netcdf.putAtt(global_ku_ncid,id_aux,units_att,rad_units);
    netcdf.putAtt(global_ku_ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(global_ku_ncid,id_aux,comment_att,'Burst phase correction applied to the SAR Ku waveforms. It is obtained by averaging the burst phase along a configurable time window to reduce noise, for antenna Rx2. The closest burst phase correction in time is used.');
end

burst_power_array_cor_rx1_name = 'burst_power_array_cor_rx1';
id_aux = netcdf.defVar(global_currentband_ncid,burst_power_array_cor_rx1_name,int32_type,[np_dimension ,n0_dimension]);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'burst_power_array_cor_rx1');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Burst power correction, for antenna Rx1');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(global_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'Burst power correction applied to the SAR Ku waveforms. It is obtained by averaging the burst power along a configurable time window to reduce noise, for antenna Rx1. The closest burst power in time is used.');

if strcmp(cnf.processing_mode,'SIN')
    burst_power_array_cor_rx2_name = 'burst_power_array_cor_rx2';
    id_aux = netcdf.defVar(global_ku_ncid,burst_power_array_cor_rx2_name,int32_type,[np_dimension ,n0_dimension]);
    netcdf.putAtt(global_ku_ncid,id_aux,std_name_att,'burst_power_array_cor_rx2');
    netcdf.putAtt(global_ku_ncid,id_aux,long_name_att,'Burst power correction, for antenna Rx1');
    netcdf.putAtt(global_ku_ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(global_ku_ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(global_ku_ncid,id_aux,comment_att,'Burst power correction applied to the SAR Ku waveforms. It is obtained by averaging the burst power along a configurable time window to reduce noise, for antenna Rx2. The closest burst power in time is used.');
end

if strcmp(cnf.processing_mode,'SIN')
instr_ext_phase_cor_name = 'instr_ext_phase_cor';
    id_aux = netcdf.defVar(global_ku_ncid,instr_ext_phase_cor_name,int32_type,n0_dimension);
    netcdf.putAtt(global_ku_ncid,id_aux,std_name_att,'instr_ext_phase_cor');
    netcdf.putAtt(global_ku_ncid,id_aux,long_name_att,'External phase correction');
    netcdf.putAtt(global_ku_ncid,id_aux,units_att,rad_units);
    netcdf.putAtt(global_ku_ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(global_ku_ncid,id_aux,comment_att,'External phase correction (SARIN only) to be added to the internal phase correction term. The external phase correction is the temperature-averaged component of external inter-channel phase difference derived from phase difference sensitive antenna subsystem, waveguides and instrument waveguide switches. The external phase correction doesn’t contain internal instrument effects of calibration coupler and duplexer which are included in the internal phase difference correction.');
end

if strcmp(cnf.processing_mode,'SIN')
instr_int_phase_cor_name = 'instr_int_phase_cor';
    id_aux = netcdf.defVar(data_ku_ncid,instr_int_phase_cor_name,int32_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'instr_int_phase_cor');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Internal phase correction computed from the CAL4');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,rad_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'Internal phase difference correction computed from the CAL-4 packets. It is set from the closest available CAL-4 packet. Applicable to SARIN only.');
end

if strcmp(cnf.processing_mode,'SIN')
phase_slope_cor_name = 'phase_slope_cor';
    id_aux = netcdf.defVar(data_ku_ncid,phase_slope_cor_name,int32_type,nr_dimension);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'phase_slope_cor');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Phase slope correction');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,rad_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'Differential group delay phase difference slope correction across the whole bandwidth (SARIN only). It is composed by fixed contributions from CCDB and by variable contributions covering differences between the CAL1 and CAL4 paths.');
end

if strcmp(cnf.processing_mode,'SIN')
interf_base_vector_name = 'interf_base_vector';
    id_aux = netcdf.defVar(data_ku_ncid,interf_base_vector_name,int32_type,[space_3D_dimension ,nr_dimension]);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'interf_base_vector');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Interferometric baseline direction vector in the satellite reference frame');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'Interferometer baseline direction vector. This is the direction vector from Tx-Rx antenna (antenna 1) reference point to Rx only antenna (antenna 2) reference point described in the satellite reference frame. The 3 components are given according to the space_3d dimension: [1] x, [2] y, [3] z.');
end

pnr_estimation_name = 'pnr_estimation_name';
id_aux = netcdf.defVar(data_currentband_ncid,pnr_estimation_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'pnr_estimation_name');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Peak to noise ratio estimation on the L1A burst');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Peak to Noise Ratio estimation done on L1A burst-averaged waveform as: peak / noise_level, where noise_level is estimated from the first n (configurable) samples.');

ptr_main_lobe_width_name = 'ptr_main_lobe_width_name';
id_aux = netcdf.defVar(data_currentband_ncid,ptr_main_lobe_width_name,int32_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'ptr_main_lobe_width_name');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'PTR main lobe width');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Width of the CAL1 waveform main lobe at -3dB.');

sig0_scaling_factor_name = 'sig0_scaling_factor';
id_aux = netcdf.defVar(data_currentband_ncid,sig0_scaling_factor_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'sig0_scaling_factor');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Power scaling: scaling factor for sigma-0 evaluation (from antenna flange)');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'The scaling factor in order to retrieve sigma-0 from L1B waveform. It includes the antenna gains and all the geometry factors satellite-surface, according to the radar equation. It is not applied to the l1b waveforms. It has to be applied to the re-tracked value of the (power_waveform*waveform_scale_factor).');

snr_instr_estimation_name = 'snr_instr_estimation';
id_aux = netcdf.defVar(data_currentband_ncid,snr_instr_estimation_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'snr_instr_estimation');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Instrument estimated SNR from telemetry');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,dB_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Instrument estimated SNR from Telemetry, converted from engineering to physical units. The SNR is estimated on tracking waveform (N-2). In the L1B HR it is the average of the snr_instr_estimation at 140Hz of all contributing Doppler beams within the stack.');

%----------------N. Look related variables --------------------------
look_counter_name = 'look_counter';
id_aux = netcdf.defVar(data_currentband_ncid,look_counter_name,int32_type,[nl_dimension ,nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'look_counter');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Burst number from L1A');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Look identification. Copied from the L1A variable tm_source_sequence_counter.');

look_time_name = 'look_time';
id_aux = netcdf.defVar(data_currentband_ncid,look_time_name,double_type,[nl_dimension ,nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'look_time');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Time in UTC: seconds since 1 Jan 2000');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Look time stamp. Copied from the L1A variable time. [tai_utc_difference] is the difference between TAI and UTC reference time (seconds) for the first measurement of the data set.');

look_i_samples_rx1_name = 'look_i_samples_rx1';
id_aux = netcdf.defVar(data_currentband_ncid,look_i_samples_rx1_name,int8_type,[ns_dimension nl_dimension nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'look_i_samples_rx1');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'I-samples for SAR L1B-S looks for antenna Rx1, arranged in stacks of [looks x samples] elements. I-samples are scaled to range [- 127, +127]');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'The in-phase component of each L1B-S look for antenna Rx1. Each look is a fully calibrated, high resolution complex waveform. Each look within the stack is: (a) given in the time domain, (b) aligned within the stack (slant range, Doppler range, window delay misalignments corrections applied), (c) fully calibrated. A final scaling, given in the variable iq_scale_factor, is applied in order to best fit the i- component into 1 byte.');

if strcmp(cnf.processing_mode,'SIN')
    look_i_samples_rx2_name = 'look_i_samples_rx2';
    id_aux = netcdf.defVar(data_ku_ncid,look_i_samples_rx2_name,int8_type,[ns_dimension nl_dimension nr_dimension]);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'look_i_samples_rx2');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'I-samples for SAR L1B-S looks for antenna Rx2, arranged in stacks of [looks x samples] elements. I-samples are scaled to range [- 127, +127]');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,number_units);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'The in-phase component of each L1B-S look for antenna Rx2. Each look is a fully calibrated, high resolution complex waveform. Each look within the stack is: (a) given in the time domain, (b) aligned within the stack (slant range, Doppler range, window delay misalignments corrections applied), (c) fully calibrated. A final scaling, given in the variable iq_scale_factor, is applied in order to best fit the i- component into 1 byte.');
end

look_q_samples_rx1_name = 'look_q_samples_rx1';
id_aux = netcdf.defVar(data_currentband_ncid,look_q_samples_rx1_name,int8_type,[ns_dimension nl_dimension nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'look_q_samples_rx1');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Q-samples for SAR L1B-S looks for antenna Rx1, arranged in stacks of looksxsamples elements. Q-samples are scaled to range [-127, +127]');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'The quadrature component of each L1B-S look for antenna Rx1. Each look is a fully calibrated, high resolution complex waveform. Each look within the stack is: (a) given in the time domain, (b) aligned within the stack (slant range, Doppler range, window delay misalignments corrections applied), (c) fully calibrated. A final scaling, given in the variable iq_scale_factor, is applied in order to best fit the q- component into 1 byte.');

if strcmp(cnf.processing_mode,'SIN')
    look_q_samples_rx2_name = 'look_q_samples_rx2';
    id_aux = netcdf.defVar(data_ku_ncid,look_q_samples_rx2_name,int8_type,[ns_dimension nl_dimension nr_dimension]);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'look_q_samples_rx2');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Q-samples for SAR L1B-S looks for antenna Rx2, arranged in stacks of looksxsamples elements. Q-samples are scaled to range [-127, +127]');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,number_units);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'The quadrature component of each L1B-S look for antenna Rx2. Each look is a fully calibrated, high resolution complex waveform. Each look within the stack is: (a) given in the time domain, (b) aligned within the stack (slant range, Doppler range, window delay misalignments corrections applied), (c) fully calibrated. A final scaling, given in the variable iq_scale_factor, is applied in order to best fit the q- component into 1 byte.');
end

iq_scale_factor_rx1_name = 'iq_scale_factor_rx1';
id_aux = netcdf.defVar(data_currentband_ncid,iq_scale_factor_rx1_name,float_type,[nl_dimension ,nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'iq_scale_factor_rx1');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'I and Q scale factor for antenna Rx1, to convert I and Q samples from [-127, +127] [or [-32766, +32766]] to amplitude at antenna flange');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,Vd1_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'The scaling factor in order to convert the I and Q samples from count to V for antenna Rx1. Note that the same scaling factor applies to both I and Q samples. The scaling is applied as follows: look_[iq]_samples_V(time, looks, samples) = look_[iq]_samples(time,looks, samples) * iq_scale_factor(time,looks).');

if strcmp(cnf.processing_mode,'SIN')
    iq_scale_factor_rx2_name = 'iq_scale_factor_rx2';
    id_aux = netcdf.defVar(data_ku_ncid,iq_scale_factor_rx2_name,float_type,[nl_dimension nr_dimension]);
    netcdf.putAtt(data_ku_ncid,id_aux,std_name_att,'iq_scale_factor_rx2');
    netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'I and Q scale factor for antenna Rx2, to convert I and Q samples from [-127, +127] [or [-32766, +32766]] to amplitude at antenna flange');
    netcdf.putAtt(data_ku_ncid,id_aux,units_att,Vd1_units);
    netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'The scaling factor in order to convert the I and Q samples from count to V for antenna Rx2. Note that the same scaling factor applies to both I and Q samples. The scaling is applied as follows: look_[iq]_samples_V(time, looks, samples) = look_[iq]_samples(time,looks, samples) * iq_scale_factor(time,looks).');
end

%----------------O. Look characterization variables --------------------------
look_index_name = 'look_index';
id_aux = netcdf.defVar(data_currentband_ncid,look_index_name,int8_type,[nl_dimension nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'look_index');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'The look index (-32, +31)');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Number that indicates, for each contributing look in the stack, the position of the beam within the burst it belongs.');

look_angle_name = 'look_angle';
id_aux = netcdf.defVar(data_currentband_ncid,look_angle_name,int16_type,[nl_dimension nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'look_angle');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Look angle associated to the look');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,rad_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'It is the angle between: (a) perpendicular from the satellite CoM to the surface, (b) direction satellite - surface location. The look angle depends purely on geometry.');

doppler_angle_name = 'doppler_angle';
id_aux = netcdf.defVar(data_currentband_ncid,doppler_angle_name,int16_type,[nl_dimension nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'doppler_angle');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Doppler angle associated to the look');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,rad_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'It is the angle between: (a) perpendicular to the velocity vector, (b) direction satellite - surface location. The Doppler angle depends on velocity vector and on geometry.');

pointing_angle_name = 'pointing_angle';
id_aux = netcdf.defVar(data_currentband_ncid,pointing_angle_name,int16_type,[nl_dimension nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'pointing_angle');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Pointing angle associated to the look');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,rad_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'It is the angle between: (a) antenna boresight direction, (b) direction satellite - surface location. The pointing angle depends on geometry and attitude (roll and pitch).');

slant_range_correction_applied_name = 'slant_range_correction_applied';
id_aux = netcdf.defVar(data_currentband_ncid,slant_range_correction_applied_name,int32_type,[nl_dimension nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'slant_range_correction_applied');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Slant range correction applied to the look');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-3);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Slant range correction applied to the look.');

doppler_correction_applied_name = 'doppler_correction_applied';
id_aux = netcdf.defVar(data_currentband_ncid,doppler_correction_applied_name,int32_type,[nl_dimension nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'doppler_correction_applied');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Doppler correction applied to the look');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,meters_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-3);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Doppler range correction applied to the look.');

stack_mask_name = 'stack_mask';
id_aux = netcdf.defVar(data_currentband_ncid,stack_mask_name,uint8_type,[nl_dimension nr_dimension]);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'stack_mask');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Zero-mask for all looks in the stack');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,2.e0);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'The zero-mask applied to the stack before multilooking. It is NOT applied at L1B-S. Each element of the mask refers to a look in the stack and indicates the index of the first sample set to zero.');

%----------R. Rate & collocation variables --------------------------
posting_rate_name = 'posting_rate';
id_aux = netcdf.defVar(global_currentband_ncid,posting_rate_name,float_type,n0_dimension);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'posting_rate_');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Waveform spacing in frequency');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,Hz_units);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'The approximate waveform spacing in frequency (Hz). The actual spacing will vary slightly along the orbit depending on orbital and processing parameters.');

ku_ka_collocation_flag_name = 'ku_ka_collocation_flag';
id_aux = netcdf.defVar(global_currentband_ncid,ku_ka_collocation_flag_name,int8_type,n0_dimension);
netcdf.putAtt(global_currentband_ncid,id_aux,std_name_att,'ku_ka_collocation_flag');
netcdf.putAtt(global_currentband_ncid,id_aux,long_name_att,'Flag indicating whether the Ku and Ka measurements are collocated');
netcdf.putAtt(global_currentband_ncid,id_aux,units_att,number_units);
netcdf.putAtt(global_currentband_ncid,id_aux,comment_att,'Flags indicating whether the Ku and Ka measurements are collocated. Only used in HR; the default value for the rest of modes is collocated. The same value is used for the Ku and Ka subgroups.');

%----------T. Thermistor temperature variables --------------------------
thr0_txrf_ku_name = 'thr0_txrf_ku';
id_aux = netcdf.defVar(data_currentband_ncid,thr0_txrf_ku_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr0_txrf_ku');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR0: TxRF_Ku');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR0: TxRF_Ku. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr1_txrf_ka_name = 'thr1_txrf_ka';
id_aux = netcdf.defVar(data_currentband_ncid,thr1_txrf_ka_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr1_txrf_ka');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR1: TxRF_Ka');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR1: TxRF_Ka. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr2_rxrf_ku1_name = 'thr2_rxrf_ku1';
id_aux = netcdf.defVar(data_currentband_ncid,thr2_rxrf_ku1_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr2_rxrf_ku1');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR2: RxRF_Ku');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR2: RxRF_Ku. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr3_rxrf_ku2_name = 'thr3_rxrf_ku2';
id_aux = netcdf.defVar(data_currentband_ncid,thr3_rxrf_ku2_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr3_rxrf_ku2');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR3: RxRF_Ku');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR3: RxRF_Ku. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr4_rxrf_ka_name = 'thr4_rxrf_ka';
id_aux = netcdf.defVar(data_currentband_ncid,thr4_rxrf_ka_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr4_rxrf_ka');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR4: RxRF_Ka');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR4: RxRF_Ka. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr5_lo_name = 'thr5_lo';
id_aux = netcdf.defVar(data_currentband_ncid,thr5_lo_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr5_lo');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR5: LO');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR5: LO. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr6_sspa_ku_name = 'thr6_sspa_ku';
id_aux = netcdf.defVar(data_currentband_ncid,thr6_sspa_ku_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr6_sspa_ku');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR6: SSPA_Ku');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR6: SSPA_Ku. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr7_sspa_ka_name = 'thr7_sspa_ka';
id_aux = netcdf.defVar(data_currentband_ncid,thr7_sspa_ka_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr7_sspa_ka');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR7: SSPA_Ka');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR7: SSPA_Ka. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr8_txnum_name = 'thr8_txnum';
id_aux = netcdf.defVar(data_currentband_ncid,thr8_txnum_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr8_txnum');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR8: TxNUM');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR8: TxNUM. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr9_rxnum_ku_name = 'thr9_rxnum_ku';
id_aux = netcdf.defVar(data_currentband_ncid,thr9_rxnum_ku_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr9_rxnum_ku');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR9: RxNUM_Ku');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR9: RxNUM_Ku. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr10_rxnum_ka_name = 'thr10_rxnum_ka';
id_aux = netcdf.defVar(data_currentband_ncid,thr10_rxnum_ka_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr10_rxnum_ka');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR10: RxNUM_Ka');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR10: RxNUM_Ka. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr11_dcdc_isps_1_name = 'thr11_dcdc_isps_1';
id_aux = netcdf.defVar(data_currentband_ncid,thr11_dcdc_isps_1_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr11_dcdc_isps_1');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR11: DCDC (ISPS 1)');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR11: DCDC (ISPS 1). Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr12_dcdc_isps_2_name = 'thr12_dcdc_isps_2';
id_aux = netcdf.defVar(data_currentband_ncid,thr12_dcdc_isps_2_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr12_dcdc_isps_2');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR12: DCDC (ISPS 2)');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR12: DCDC (ISPS 2). Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr13_formatter_name = 'thr13_formatter';
id_aux = netcdf.defVar(data_currentband_ncid,thr13_formatter_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr13_formatter');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR13: FORMATTER');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR13: FORMATTER. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

thr14_sequencer_name = 'thr14_sequencer';
id_aux = netcdf.defVar(data_currentband_ncid,thr14_sequencer_name,int16_type,nr_dimension);
netcdf.putAtt(data_currentband_ncid,id_aux,std_name_att,'thr14_sequencer');
netcdf.putAtt(data_currentband_ncid,id_aux,long_name_att,'Temperature reading from THR14: SEQUENCER');
netcdf.putAtt(data_currentband_ncid,id_aux,units_att,celsius_units);
netcdf.putAtt(data_currentband_ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(data_currentband_ncid,id_aux,comment_att,'Temperature reading from THR14: SEQUENCER. Converted from digital values to degree celsius using a dedicated table in the Characterisation array file.');

%--- End of variables definition ---

netcdf.endDef(ncid);

% var_id = netcdf.inqVarID(ncid,'zero_padding_l1b_echo');
% netcdf.putVar(ncid,var_id,0,cnf.zp_fact_range);

netcdf.close(ncid);
% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);
