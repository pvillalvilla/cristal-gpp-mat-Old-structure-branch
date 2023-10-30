%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the CODING & PACKING 
% algorithm for Level-1BS product using Sentinel-3 like format
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
function [files] = create_NetCDF_FFL1B(files, N_surfs_loc_estimated, type, cnf, chd, cst, resize)

%global N_samples   semi_major_axis_cst flat_coeff_cst;
global bw_ku_chd 
% zp_fact_range_cnf  sec_in_day_cst pi_cst
%global mission mode N_max_beams_stack_chd
%global compute_L1_stadistics_cnf include_wfms_aligned
%global optional_ext_file_flag file_ext_string
%added by EM: 04.10.2016 
%global ACDC_application_cnf cnf_p_ACDC
%global netcdf_type
%global processing_mode_cnf
% t6 = tic;

if(strcmp(type,'FF-SL')||strcmp(type,'FF-ML'))
    disp(['writting ' type ' product']);
else
    disp('type should indicate FF-SL or FF-ML Fully focused products');
    return;
end

date_creation = datestr(now, '_yyyymmddTHHMMSS_');

% JPLZ: the following part is for creating the L1B-FF filename, but it is
% duplicated? It was already done in GPPICE_processor.m. I comment it for
% the moment. I see no use for it here, we already bring the filenames as
% inputs

% L1A_date_find = strfind(files.filename_L1A,'_201');
% L1A_date_find = strfind(files.filename_L1A,'_202'); %JPLZ: date for IRIS OB chronogram simulations is 202X
% L1A_mode_find = strfind(files.filename_L1A, 'PICE_SIRS_');

% if (~resize)
%     files.filename_L1BFF = strcat(files.outputPath,'PICE_SIRS_FF_',...
%         files.filename_L1A((L1A_mode_find(1)+9):(L1A_mode_find(end)+13)), '_1B_',...
%         files.filename_L1A((L1A_date_find(1)+1):(L1A_date_find(1)+30)),...
%         '_isd','_long','.nc');
% else
%     files.filename_L1BFF = strcat(files.outputPath,'PICE_SIRS_FF_',...
%         files.filename_L1A((L1A_mode_find(1)+9):(L1A_mode_find(end)+13)), '_1B_',...
%         files.filename_L1A((L1A_date_find(1)+1):(L1A_date_find(1)+30)),...
%         '_isd','.nc');
% end

%{
files.filename_L1B = strcat(files.outputPath,'SR_1_SRA____',...
                                files.filename_L1A(end-37:end-3),...
                                '_isd','.nc');



%}
if(strcmp(type,'FF-SL'))
    ncid = netcdf.create(files.filename_L1BFF_SL,'CLOBBER');
elseif (strcmp(type,'FF-ML'))
    ncid = netcdf.create(files.filename_L1BFF_ML,'CLOBBER');
end

% switch netcdf_type
%     case 'netcdf4'
%         ncid = netcdf.create(files.filename_netCDF,'NETCDF4');
%     case 'netcdf3'
%         ncid = netcdf.create(files.filename_netCDF,'CLASSIC_MODEL');        
% end

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

% nb_dimension = 'time_l1b_echo';
% nl_dimension = 'Nl';
% nl_dimension_size = N_max_beams_stack_chd;
% ns_dimension = 'Ns';
% space_3D_dimension = 'space_3D';
% space_3D_dimension_size = 3;

nl_dimension = netcdf.defDim(ncid,'nl',chd.N_max_beams_stack);
nb_dimension = netcdf.defDim(ncid,'nb',N_surfs_loc_estimated); % Number of records
np_dimension = netcdf.defDim(ncid,'np',chd.N_pulses_burst); % Number of pulses
ns_dimension = netcdf.defDim(ncid,'echo_sample_ind',chd.N_samples_sar*cnf.zp_fact_range);
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
rad_units = 'rad';
percent_units='percent';


int8_type = 'NC_BYTE';
uint8_type = netcdf.getConstant('ubyte');
int16_type = 'NC_SHORT';
uint16_type = netcdf.getConstant('ushort');
int32_type = 'NC_INT';
uint32_type = netcdf.getConstant('uint');
%int64_type = netcdf.getConstant('int64');
int64_type = 'NC_INT64';
float_type = 'NC_FLOAT';
double_type= 'NC_DOUBLE';



%% PACKING L1B


%----------A. Time variables ----------------------------------------------
l1_mode_id_name = 'l1_mode_id';
id_aux      = netcdf.defVar(ncid,l1_mode_id_name,int8_type,nb_dimension);
            netcdf.putAtt(ncid,id_aux,long_name_att,'L1 mode ID');
            netcdf.putAtt(ncid,id_aux,comment_att,'L1 Mode Identifier. Each L1A record type has a unique ID, as follows: 0 SAR, 1 SARIn');



time_l1b_echo_name = 'time_l1b_echo';
id_aux      = netcdf.defVar(ncid,time_l1b_echo_name,double_type,nb_dimension);
            netcdf.putAtt(ncid,id_aux,long_name_att,'UTC ');
            netcdf.putAtt(ncid,id_aux,units_att,seconds_units);

seq_count_l1b_echo_name = 'seq_count_l1b_echo';
id_aux 		= netcdf.defVar(ncid,seq_count_l1b_echo_name,int16_type, nb_dimension);
            netcdf.putAtt(ncid,id_aux,'FillValue',65535);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Source Sequence Count: Same as Burst number ');
            netcdf.putAtt(ncid,id_aux,units_att,number_units);			
			

UTC_day_l1b_echo_name = 'UTC_day_l1b_echo';
id_aux      = netcdf.defVar(ncid,UTC_day_l1b_echo_name,int16_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Days since 2000-01-01 00:00:00.0+00:00 (Ku-band)');
            netcdf.putAtt(ncid,id_aux,units_att,day_units);
			netcdf.putAtt(ncid,id_aux,comment_att,'days elapsed since 2000-01-01. To be used to link with L1 and L2 records (time_l1b provides the number of seconds since 2000-01-01).');



UTC_sec_l1b_echo_name = 'UTC_sec_l1b_echo';
id_aux = netcdf.defVar(ncid,UTC_sec_l1b_echo_name,double_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,1.844674407370960e+19);
netcdf.putAtt(ncid,id_aux,long_name_att,'Seconds in the day UTC, with microsecond resolution (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,comment_att,'seconds in the day. To be used to link L1 and L2 records (time_l1b provides the number of seconds since 2000-01-01).');

if(strcmp(type,'FF-ML'))
    
    time_surf_ML_name = 'time_surf_ML';
    id_aux      = netcdf.defVar(ncid,time_surf_ML_name,double_type,nb_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'UTC ');
    netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
    
end
%----------B. Orbit and attitude variables ------------------------------
lat_l1b_echo_name = 'lat_l1b_echo';
id_aux = netcdf.defVar(ncid,lat_l1b_echo_name,int32_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,std_name_att,'latitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'latitude (positive N, negative S) (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Latitude of measurement [-90, +90]: Positive at Nord, Negative at South');


lon_l1b_echo_name = 'lon_l1b_echo';
id_aux = netcdf.defVar(ncid,lon_l1b_echo_name,int32_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,std_name_att,'longitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'longitude (positive E, negative W) (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'longitude of measurement [-180, +180]: Positive at East, Negative at West');


alt_l1b_echo_name = 'alt_l1b_echo';
id_aux = netcdf.defVar(ncid,alt_l1b_echo_name,int32_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'altitude of satellite');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Altitude of the satellite Centre of Mass');


orb_alt_rate_l1b_echo_name = 'orb_alt_rate_l1b_echo';
id_aux = netcdf.defVar(ncid,orb_alt_rate_l1b_echo_name,int16_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
netcdf.putAtt(ncid,id_aux,long_name_att,'orbital altitude rate');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Instantaneous altitude rate at the Centre of Mass');


satellite_mispointing_l1b_echo_name = 'satellite_mispointing_l1b_echo';
id_aux = netcdf.defVar(ncid,satellite_mispointing_l1b_echo_name,int32_type,[space_3D_dimension ,nb_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'Mispointing angle, measures by STRs: [1] Roll, [2] Pitch, [3] Yaw (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-7);
netcdf.putAtt(ncid,id_aux,comment_att,'Attitude mispointing, measured by STRs and post-processed by AOCS or by ground facility. The 3 components are given according to the ''space_3D'' dimension: [1] Roll, [2] Pitch, [3] Yaw. This variable includes the "mispointing bias" given by the variable mispointing_bias_ku. Note: nominal pointing is at satellite nadir (antenna perpendicular to ellipsoid) and corresponds to: roll = pitch = yaw = 0');

if(strcmp(type,'FF-ML'))
    lat_surf_ML_name = 'lat_surf_ML';
    id_aux = netcdf.defVar(ncid,lat_surf_ML_name,int32_type,nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,214748364);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,std_name_att,'latitude');
    netcdf.putAtt(ncid,id_aux,long_name_att,'latitude (positive N, negative S) (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(ncid,id_aux,comment_att,'ML Latitude of measurement [-90, +90]: Positive at Nord, Negative at South');
    
    
    lon_surf_ML_name = 'lon_surf_ML';
    id_aux = netcdf.defVar(ncid,lon_surf_ML_name,int32_type,nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,214748364);
    netcdf.putAtt(ncid,id_aux,std_name_att,'longitude ML');
    netcdf.putAtt(ncid,id_aux,long_name_att,'longitude (positive E, negative W) (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,comment_att,'ML longitude of measurement [-180, +180]: Positive at East, Negative at West');
    
    
    alt_surf_ML_name = 'alt_surf_ML';
    id_aux = netcdf.defVar(ncid,alt_surf_ML_name,int32_type,nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,214748364);
    netcdf.putAtt(ncid,id_aux,long_name_att,'altitude of satellite ML');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
    netcdf.putAtt(ncid,id_aux,comment_att,'ML Altitude of the ML surface');
    
    
    alt_rate_surf_ML_name = 'alt_rate_surf_ML';
    id_aux = netcdf.defVar(ncid,alt_rate_surf_ML_name,int16_type,nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,32767);
    netcdf.putAtt(ncid,id_aux,long_name_att,'orbital altitude rate ML');
    netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0);
    netcdf.putAtt(ncid,id_aux,comment_att,'ML Instantaneous altitude rate of the ML surface');
end

%----------C. Flag time variables --------------------------------------

%----------D. Position/Velocity variables ------------------------------
x_pos_l1b_echo_name = 'x_pos_l1b_echo';
id_aux = netcdf.defVar(ncid,x_pos_l1b_echo_name,double_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite altitude-x component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

y_pos_l1b_echo_name = 'y_pos_l1b_echo';
id_aux = netcdf.defVar(ncid,y_pos_l1b_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite altitude-y component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

z_pos_l1b_echo_name = 'z_pos_l1b_echo';
id_aux = netcdf.defVar(ncid,z_pos_l1b_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite altitude-z component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

x_vel_l1b_echo_name = 'x_vel_l1b_echo';
id_aux = netcdf.defVar(ncid,x_vel_l1b_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite velocity-x component');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);

y_vel_l1b_echo_name = 'y_vel_l1b_echo';
id_aux = netcdf.defVar(ncid,y_vel_l1b_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite velocity-y component');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);

z_vel_l1b_echo_name = 'z_vel_l1b_echo';
id_aux = netcdf.defVar(ncid,z_vel_l1b_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite velocity-z component');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);

if(strcmp(type,'FF-ML'))    
    x_vel_sat_surf_ML_name = 'x_vel_sat_surf_ML';
    id_aux = netcdf.defVar(ncid,x_vel_sat_surf_ML_name,double_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ML surface satellite velocity-x component');
    netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);
    
    y_vel_sat_surf_ML_name = 'y_vel_sat_surf_ML';
    id_aux = netcdf.defVar(ncid,y_vel_sat_surf_ML_name,double_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ML surface satellite velocity-y component');
    netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);
    
    z_vel_sat_surf_ML_name = 'z_vel_sat_surf_ML';
    id_aux = netcdf.defVar(ncid,z_vel_sat_surf_ML_name,double_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ML surface satellite velocity-z component');
    netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units); 
end

%----------E. Navigation Bulletin ------------------------------
% seq_count_l1b_echo_name = 'seq_count_l1b_echo';
% id_aux = netcdf.defVar(ncid,seq_count_l1b_echo_name,int32_type,nb_dimension);        
% switch netcdf_type
%     case 'netcdf4'
%           id_aux = netcdf.defVar(ncid,seq_count_l1b_echo_name,uint16_type,nb_dimension);        
%     case 'netcdf3'
%            %increase the size
%           id_aux = netcdf.defVar(ncid,seq_count_l1b_echo_name,int32_type,nb_dimension);
% end    
%             netcdf.defVarFill(ncid,id_aux,false,65535);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Sequence count');
% netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

%----------F. Operating instrument & tracking --------------------------
oper_instr_l1b_echo_name = 'oper_instr_l1b_echo';
id_aux = netcdf.defVar(ncid,oper_instr_l1b_echo_name,int8_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,127);
netcdf.putAtt(ncid,id_aux,long_name_att,'Operating instrument');
netcdf.putAtt(ncid,id_aux,flag_values_att,'0b,1b');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'A,B (Sentinel-3) / Nominal, Redundant (CryoSat-2)');
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement. For Sentinel-3: Instrument A stands for SRAL Nominal and instrument B stands for SRAL Redundant');    

SAR_mode_l1b_echo_name='SAR_mode_l1b_echo';
id_aux = netcdf.defVar(ncid,SAR_mode_l1b_echo_name,int8_type, nb_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'SAR mode identifier');
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');
netcdf.putAtt(ncid,id_aux,flag_values_att,'0b,1b,2b (Sentinel-3) / 0b, 1b (Cryosat-2)');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'closed_loop, open_loop, open_loop_fixed_gain (Sentinel-3) / closed, open (CryoSat-2)');

%---------- G. H0, COR2 and AGC ----------------------------------------
h0_applied_l1b_echo_name = 'h0_applied_l1b_echo';
id_aux = netcdf.defVar(ncid,h0_applied_l1b_echo_name,int32_type, nb_dimension);
% switch netcdf_type
%     case 'netcdf4'
%         id_aux = netcdf.defVar(ncid,h0_applied_l1b_echo_name,int64_type, nb_dimension);
%     case 'netcdf3'
%         id_aux = netcdf.defVar(ncid,h0_applied_l1b_echo_name,int64_type, nb_dimension);
% end
%             netcdf.defVarFill(ncid,id_aux,false,4294967295);
netcdf.putAtt(ncid,id_aux,long_name_att,'Applied altitude command H0');
netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d64d1e9_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

cor2_applied_l1b_echo_name = 'cor2_applied_l1b_echo';
id_aux = netcdf.defVar(ncid,cor2_applied_l1b_echo_name,int16_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
netcdf.putAtt(ncid,id_aux,long_name_att,'Applied altitude command COR2');
netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d1024d1e9_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');


agccode_ku_l1b_echo_name = 'agccode_ku_l1b_echo';
id_aux = netcdf.defVar(ncid,agccode_ku_l1b_echo_name,int8_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'AGCCODE for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

%---------- H. Surface type -----------------------------------------------
surf_type_l1b_echo_name = 'surf_type_l1b_echo';
id_aux = netcdf.defVar(ncid,surf_type_l1b_echo_name,int8_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,127);
netcdf.putAtt(ncid,id_aux,long_name_att,'Altimeter surface type');
netcdf.putAtt(ncid,id_aux,flag_values_att,'0,1,2,3');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'open_ocean or semi-enclosed_seas, enclosed_seas or lakes, continental_ice, land,Transponder');
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

%---------- I. Altimeter range and Corrections ----------------------------
range_ku_l1b_echo_name = 'range_ku_l1b_echo';
id_aux = netcdf.defVar(ncid,range_ku_l1b_echo_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Corrected range for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Reference range corrected for USO frequency drift and internal path correction');


uso_cor_l1b_echo_name = 'uso_cor_l1b_echo';
id_aux = netcdf.defVar(ncid,uso_cor_l1b_echo_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'USO frequency drift correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

int_path_cor_ku_l1b_echo_name = 'int_path_cor_ku_l1b_echo';
id_aux = netcdf.defVar(ncid,int_path_cor_ku_l1b_echo_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Internal path correction for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

range_rate_l1b_echo_name = 'range_rate_l1b_echo';
id_aux = netcdf.defVar(ncid,range_rate_l1b_echo_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Range rate');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

if(strcmp(type,'FF-ML'))
    
    range_ML_name = 'range_ML';
    id_aux = netcdf.defVar(ncid,range_ML_name,int32_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ML corrected range for Ku-Ka band');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
    netcdf.putAtt(ncid,id_aux,comment_att,'ML reference range corrected for USO frequency drift and internal path correction');
    
end

%---------- J. AGC and Sigma0 scalings --------------------------------
scale_factor_ku_l1b_echo_name = 'scale_factor_ku_l1b_echo';
id_aux = netcdf.defVar(ncid,scale_factor_ku_l1b_echo_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Scaling factor for sigma0 evaluation');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'This is a scaling factor in order to retrieve sigma-0 from the L1B waveform. It includes antenna gains and geometry satellite - surface. It is not applied to the L1B waveforms');

%---------- K. Stack characterization--------------------------------------
nb_stack_l1b_echo_name = 'nb_stack_l1b_echo';
id_aux = netcdf.defVar(ncid,nb_stack_l1b_echo_name,int32_type, nb_dimension);
% switch netcdf_type
%     case 'netcdf4'
%         id_aux = netcdf.defVar(ncid,nb_stack_l1b_echo_name,uint16_type, nb_dimension);
%     case 'netcdf3'
%         id_aux = netcdf.defVar(ncid,nb_stack_l1b_echo_name,int32_type, nb_dimension);
% end
%             netcdf.defVarFill(ncid,id_aux,false,65535);
netcdf.putAtt(ncid,id_aux,long_name_att,'Number of waveforms summed in stack (contributing beams: effective number of looks or beams from stack used in multilooking; if beams with all samples set to zero not to be included in multilooking they are accordingly not accounted in this number; if there are some gaps in-between mask set to zero all the samples related to those beams and they will not be contributing to the multilooking )');
netcdf.putAtt(ncid,id_aux,units_att,number_units);

nb_stack_start_stop_l1b_echo_name = 'nb_stack_start_stop_l1b_echo';
id_aux = netcdf.defVar(ncid,nb_stack_start_stop_l1b_echo_name,int32_type, nb_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Number of waveforms in stack (considering the number of beams/looks from start and stop beams: they correspond to the first and last beams at edges of stack. If option of discarding beams with all-zeros samples these correspond to the first and last beams with non-all zeros samples (if there exist gaps in between: mask setting to zero all samples, these in-between beams are considered anyway in the "nb_stack_start_stop_l1b_echo"), otherwise they correspond to the very first and very last beam in the stack ). This number of beams will be useful to construct accordingly the modelled stack for retracking.');
netcdf.putAtt(ncid,id_aux,units_att,number_units);

look_angle_start_l1b_echo_name = 'look_angle_start_l1b_echo';
id_aux = netcdf.defVar(ncid,look_angle_start_l1b_echo_name,int16_type,nb_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of first look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Look angle of the first contributing look (non-0 weight) to the L1B waveform');

look_angle_stop_l1b_echo_name = 'look_angle_stop_l1b_echo';
id_aux = netcdf.defVar(ncid,look_angle_stop_l1b_echo_name,int16_type,nb_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of last look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Look angle of the last contributing look (non-0 weight) to the L1B waveform');

doppler_angle_start_l1b_echo_name = 'doppler_angle_start_l1b_echo';
id_aux = netcdf.defVar(ncid,doppler_angle_start_l1b_echo_name,int16_type,nb_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of first look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Doppler angle of the first contributing look (non-0 weight) to the L1B waveform');

doppler_angle_stop_l1b_echo_name = 'doppler_angle_stop_l1b_echo';
id_aux = netcdf.defVar(ncid,doppler_angle_stop_l1b_echo_name,int16_type,nb_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of last contributing look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Doppler angle of the last contributing look to the L1B waveform');

pointing_angle_start_l1b_echo_name = 'pointing_angle_start_l1b_echo';
id_aux = netcdf.defVar(ncid,pointing_angle_start_l1b_echo_name,int16_type,nb_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of first contributing look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Pointing angle of the first contributing look to the L1B waveform');

pointing_angle_stop_l1b_echo_name = 'pointing_angle_stop_l1b_echo';
id_aux = netcdf.defVar(ncid,pointing_angle_stop_l1b_echo_name,int16_type,nb_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of last contributing look  (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,rad_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Pointing angle of the last contributing look to the L1B waveform');

if cnf.compute_L1_stadistics
    skew_stack_l1b_echo_name = 'skew_stack_l1b_echo';
    id_aux = netcdf.defVar(ncid,skew_stack_l1b_echo_name,int32_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Skewness of stack');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,comment_att,'Skewness of the Gaussian that fits the integrated power of the looks within a stack. The skewness indicates how symmetric or asymmetric the power within the stack is');
    
    kurt_stack_l1b_echo_name = 'kurt_stack_l1b_echo';
    id_aux = netcdf.defVar(ncid,kurt_stack_l1b_echo_name,int32_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Kurtosis of stack');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,comment_att,'Kurtosis of the Gaussian that fits the integrated power of the looks within a stack. Kurtosis is a measure of peakiness');
    
    stdev_stack_l1b_echo_name = 'stdev_stack_l1b_echo';
    id_aux = netcdf.defVar(ncid,stdev_stack_l1b_echo_name,int32_type, nb_dimension);
    % switch netcdf_type
    %     case 'netcdf4'
    %         id_aux = netcdf.defVar(ncid,stdev_stack_l1b_echo_name,uint32_type, nb_dimension);
    %     case 'netcdf3'
    %         id_aux = netcdf.defVar(ncid,stdev_stack_l1b_echo_name,int64_type, nb_dimension);
    % end
    %             netcdf.defVarFill(ncid,id_aux,false,4294967295);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Gaussian Power fitting: STD wrt look angle (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,rad_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,comment_att,'Standard deviation of the Gaussian that fits the integrated power of the looks within a stack. It is given with respect to the look angle. The width at -3dB of this Gaussian can be retrieved the following way: width_3db = 2*sqrt(2*ln2)*gaussian_fitting_std');
    
    gaussian_fitting_centre_look_l1b_echo_name = 'gaussian_fitting_centre_look_l1b_echo';
    id_aux = netcdf.defVar(ncid,gaussian_fitting_centre_look_l1b_echo_name,int16_type,nb_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Gaussian Power fitting: centre wrt look angle (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,rad_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,comment_att,'Position of the center of the Gaussian that fits the integrated power of the looks within a stack, with respect to look angle');
    
    gaussian_fitting_centre_pointing_l1b_echo_name = 'gaussian_fitting_centre_pointing_l1b_echo';
    id_aux = netcdf.defVar(ncid,gaussian_fitting_centre_pointing_l1b_echo_name,int16_type,nb_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Gaussian Power fitting: centre wrt pointing angle (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,rad_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
    netcdf.putAtt(ncid,id_aux,comment_att,'Position of the centre of the Gaussian that fits the integrated power of the looks within a stack, with respect to pointing angle');
    
    beam_form_l1b_echo_name = 'beam_form_l1b_echo';
    id_aux = netcdf.defVar(ncid,beam_form_l1b_echo_name,int32_type,nb_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Flag on beam formation quality in stack');
    netcdf.putAtt(ncid,id_aux,units_att,percent_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-2);
    netcdf.putAtt(ncid,id_aux,comment_att,'Beam formation quality in percentage');
end


%------------- L. Altimeter engineering variables -------------------------
altimeter_clock_l1b_echo_name = 'altimeter_clock_l1b_echo';
id_aux = netcdf.defVar(ncid,altimeter_clock_l1b_echo_name,int32_type,nb_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Altimeter clock (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,Hz_units);
netcdf.putAtt(ncid,id_aux,add_offset_att,bw_ku_chd);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-9);
netcdf.putAtt(ncid,id_aux,comment_att,'This is the actual altimeter clock.');

pulse_repetition_interval_l1b_echo_name = 'pri_l1b_echo';
id_aux = netcdf.defVar(ncid,pulse_repetition_interval_l1b_echo_name,int32_type,nb_dimension);
% switch netcdf_type
%     case 'netcdf4'
%         id_aux = netcdf.defVar(ncid,pulse_repetition_interval_l1b_echo_name,uint32_type,nb_dimension);
%     case 'netcdf3'
%         id_aux = netcdf.defVar(ncid,pulse_repetition_interval_l1b_echo_name,int64_type,nb_dimension);
% end
netcdf.putAtt(ncid,id_aux,long_name_att,'PRI converted into seconds (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-12);
netcdf.putAtt(ncid,id_aux,comment_att,'The ''Pulse Repetition Interval''. PRI is constant within all received pulses in a radar cycle, but it can change within consecutive radar cycles.');

%------------- M. Waveform related variables -----------------------------   
if(strcmp(type,'FF-ML'))  
    % Variable M04
    waveform_scale_factor_name = 'waveform_scale_factor';
    id_aux = netcdf.defVar(ncid,waveform_scale_factor_name,float_type,nb_dimension);
    netcdf.putAtt(ncid,id_aux,std_name_att,'Echo scale factor, to convert from [0-65534] to power at antenna flange (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
    netcdf.putAtt(ncid,id_aux,'coordinates','longitude latitude'); 
    netcdf.putAtt(ncid,id_aux,comment_att,'The L1B waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(time, samples) = power_waveform(time, samples) * waveform_scale_factor(time)');
    
    % Variable M10
    power_waveform_name = 'power_waveform';
    id_aux = netcdf.defVar(ncid,power_waveform_name,double_type,[ns_dimension nb_dimension]);
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,long_name_att,'L1B Fully Focused Multi-Looked power waveform: scaled 0-65534 (Ku-band)');
    netcdf.putAtt(ncid,id_aux,'coordinates','longitude latitude');
    netcdf.putAtt(ncid,id_aux,comment_att,'The L1B FF-ML Power waveform is a fully calibrated, high resolution, multi-looked waveform in the time domain. It includes: (a) all calibrations, which have been applied at L1A, (b) configuration according to processing_configuration_flags and ff_processing_configuration_flag, (c) final scaling, given in the variable ''waveform_scale_factor'', in order to best fit the waveform into 2 bytes.');
    
    if strcmp(cnf.processing_mode,'SIN')       
        waveform_scale_factor_2_name = 'waveform_scale_factor_2';
        id_aux = netcdf.defVar(ncid,waveform_scale_factor_2_name,float_type,nb_dimension);
        netcdf.putAtt(ncid,id_aux,std_name_att,'Echo scale factor, to convert from [0-65534] to power at antenna flange (Ku-band)');
        netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
        netcdf.putAtt(ncid,id_aux,'coordinates','longitude latitude');
        netcdf.putAtt(ncid,id_aux,comment_att,'The L1B waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(time, samples) = power_waveform(time, samples) * waveform_scale_factor(time)');
        
        % Variable M10
        power_waveform_2_name = 'power_waveform_2';
        id_aux = netcdf.defVar(ncid,power_waveform_2_name,double_type,[ns_dimension nb_dimension]);
        netcdf.putAtt(ncid,id_aux,units_att,number_units);
        netcdf.putAtt(ncid,id_aux,long_name_att,'L1B Fully Focused Multi-Looked power waveform: scaled 0-65534 (Ku-band)');
        netcdf.putAtt(ncid,id_aux,'coordinates','longitude latitude');
        netcdf.putAtt(ncid,id_aux,comment_att,'The L1B FF-ML Power waveform is a fully calibrated, high resolution, multi-looked waveform in the time domain. It includes: (a) all calibrations, which have been applied at L1A, (b) configuration according to processing_configuration_flags and ff_processing_configuration_flag, (c) final scaling, given in the variable ''waveform_scale_factor'', in order to best fit the waveform into 2 bytes.');
         
        waveform_scale_factor_comb_name = 'waveform_scale_factor_comb';
        id_aux = netcdf.defVar(ncid,waveform_scale_factor_comb_name,float_type,nb_dimension);
        netcdf.putAtt(ncid,id_aux,std_name_att,'Echo scale factor, to convert from [0-65534] to power at antenna flange (Ku-band)');
        netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
        netcdf.putAtt(ncid,id_aux,'coordinates','longitude latitude');
        netcdf.putAtt(ncid,id_aux,comment_att,'The L1B waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(time, samples) = power_waveform(time, samples) * waveform_scale_factor(time)');
        
        % Variable M10
        power_waveform_comb_name = 'power_waveform_comb';
        id_aux = netcdf.defVar(ncid,power_waveform_comb_name,double_type,[ns_dimension nb_dimension]);
        netcdf.putAtt(ncid,id_aux,units_att,number_units);
        netcdf.putAtt(ncid,id_aux,long_name_att,'L1B Fully Focused Multi-Looked power waveform combined: scaled 0-65534 (Ku-band)');
        netcdf.putAtt(ncid,id_aux,'coordinates','longitude latitude');
        netcdf.putAtt(ncid,id_aux,comment_att,'The L1B FF-ML Power waveform is a fully calibrated, high resolution, multi-looked waveform in the time domain. It includes: (a) all calibrations, which have been applied at L1A, (b) configuration according to processing_configuration_flags and ff_processing_configuration_flag, (c) final scaling, given in the variable ''waveform_scale_factor'', in order to best fit the waveform into 2 bytes.');
          
        
        phase_diff_ML_name = 'phase_diff_ML';
        id_aux = netcdf.defVar(ncid,phase_diff_ML_name,double_type,[ns_dimension nb_dimension]);
        netcdf.putAtt(ncid,id_aux,long_name_att,'SARIn Phase Difference');
        netcdf.putAtt(ncid,id_aux,units_att,rad_units);
        netcdf.putAtt(ncid,id_aux,comment_att,'Phase difference array, computed from the complex echoes on the 2 Rx channels. Values are between -pi and +pi radians, units are micro radians');
        
        wvfm_coh_ML_name = 'wvfm_coh_ML';
        id_aux = netcdf.defVar(ncid,wvfm_coh_ML_name,double_type,[ns_dimension nb_dimension]);
        netcdf.putAtt(ncid,id_aux,long_name_att,'SARIn Coherence');
        netcdf.putAtt(ncid,id_aux,units_att,number_units);
        netcdf.putAtt(ncid,id_aux,comment_att,'Coherence array , computed from the complex echoes on the 2 Rx channels');        
    end
    

    % Variable M05
%     pnr_estimation_name = 'pnr_estimation';
%     id_aux = netcdf.defVar(data_ku_ncid,pnr_estimation_name,int16_type,ku_rec_dimension);
%     netcdf.defVarFill(data_ku_ncid, id_aux, false, fillVal_short);
%     netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Peak to Noise Ratio estimation on the L1B waveform (Ku-band)');
%     netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1e-2);
%     netcdf.putAtt(data_ku_ncid,id_aux,units_att,dB_units);
%     netcdf.putAtt(data_ku_ncid,id_aux,'coordinates','longitude latitude'); 
%     netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'Peak to Noise Ratio estimation done on the L1b waveform as: peak / noise_level. Where noise_level is estimated from the first n (configurable) samples.');

    % Variable M06A (MULTI-LOOK VERSION)
%     sig0_scaling_factor_name = 'sig0_scaling_factor';
%     id_aux = netcdf.defVar(data_ku_ncid,sig0_scaling_factor_name,int16_type,ku_rec_dimension);
%     netcdf.defVarFill(data_ku_ncid, id_aux, false, fillVal_short);
%     netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Power scaling: scaling factor for sigma-0 evaluation (from antenna flange) (Ku-band)');
%     netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1e-2);
%     netcdf.putAtt(data_ku_ncid,id_aux,units_att,dB_units);
%     netcdf.putAtt(data_ku_ncid,id_aux,'coordinates','longitude latitude'); 
%     netcdf.putAtt(data_ku_ncid,id_aux,comment_att,' The scaling factor in order to retrieve sigma-0 from the L1B waveform. It includes the antenna gains and all the geometry factors satellite–surface, according to the radar equation. It is not applied to the L1b waveforms, it has to be applied to the re-tracked value of the (power_waveform*waveform_scale_factor)');

%if(strcmp(type,'FF-SL'))
    % Variable M06B (SINGLE-LOOK VERSION)
%     sig0_scaling_factor_name = 'sig0_scaling_factor';
%     id_aux = netcdf.defVar(data_ku_ncid,sig0_scaling_factor_name,int16_type,ku_rec_dimension);
%     netcdf.defVarFill(data_ku_ncid, id_aux, false, fillVal_short);
%     netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'Power scaling: scaling factor for sigma-0 evaluation (from antenna flange) (Ku-band)');
%     netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1e-2);
%     netcdf.putAtt(data_ku_ncid,id_aux,units_att,dB_units);
%     netcdf.putAtt(data_ku_ncid,id_aux,'coordinates','longitude latitude'); 
%     netcdf.putAtt(data_ku_ncid,id_aux,comment_att,' The scaling factor in order to retrieve sigma-0 from the L1B waveform. It includes the antenna gains and all the geometry factors satellite–surface, according to the radar equation. It is not applied to the L1b waveforms. It can be applied to the L1B SL waveform after power extraction');
%end

% % Variable M09
% snr_instr_estimation_name = 'snr_instr_estimation';
% id_aux = netcdf.defVar(data_ku_ncid,snr_instr_estimation_name,int16_type,ku_rec_dimension);
% netcdf.defVarFill(data_ku_ncid, id_aux, false, fillVal_short);
% netcdf.putAtt(data_ku_ncid,id_aux,long_name_att,'instrument estimated SNR from Telemetry');
% netcdf.putAtt(data_ku_ncid,id_aux,scale_factor_att,1e-2);
% netcdf.putAtt(data_ku_ncid,id_aux,units_att,dB_units);
% netcdf.putAtt(data_ku_ncid,id_aux,'coordinates','longitude latitude'); 
% netcdf.putAtt(data_ku_ncid,id_aux,comment_att,'instrument estimated SNR from Telemetry, converted from engineering to physical unit. The SNR is estimated on tracking waveform (N-2). In the L1B HR and L1B-S it is the average of the snr_instr_estimation of all the contributing looks within the stack.');

end
 
if(strcmp(type,'FF-SL'))
    
    waveform_scale_factor_l1b_echo_name = 'waveform_scale_factor_l1b_echo';
    id_aux = netcdf.defVar(ncid,waveform_scale_factor_l1b_echo_name,float_type,nb_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Echo Scale Factor, to convert from [0-65535]	to Power at	antenna	flange');
    netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The L1B waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = i2q2_meas_ku_l1b_echo(ku_rec, Ns) * waveform_scale_factor_l1b_echo(ku_rec)');
       
    i2q2_meas_ku_l1b_echo_name = 'i2q2_meas_ku_l1b_echo';
    id_aux = netcdf.defVar(ncid,i2q2_meas_ku_l1b_echo_name,double_type,[ns_dimension nb_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'SAR Power Echo waveform: scaled	0-65535 (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The SAR L1B Power waveforms is a fully calibrated, FF, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_l1b_echo'', in order to best fit the waveform into 2 bytes');
    
    if strcmp(cnf.processing_mode,'SIN')
        waveform_scale_factor_l1b_echo_name_2 = 'waveform_scale_factor_l1b_echo_2';
        id_aux = netcdf.defVar(ncid,waveform_scale_factor_l1b_echo_name_2,float_type,nb_dimension);
        netcdf.putAtt(ncid,id_aux,long_name_att,'Echo Scale Factor, to convert from [0-65535]	to Power at	antenna	flange');
        netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
        netcdf.putAtt(ncid,id_aux,comment_att,'The L1B waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = i2q2_meas_ku_l1b_echo(ku_rec, Ns) * waveform_scale_factor_l1b_echo(ku_rec)');
        
        i2q2_meas_ku_l1b_echo_name_2 = 'i2q2_meas_ku_l1b_echo_2';
        id_aux = netcdf.defVar(ncid,i2q2_meas_ku_l1b_echo_name_2,double_type,[ns_dimension nb_dimension]);
        netcdf.putAtt(ncid,id_aux,long_name_att,'SARIn Power Echo waveform: scaled	0-65535 (Ku-band)');
        netcdf.putAtt(ncid,id_aux,units_att,number_units);
        netcdf.putAtt(ncid,id_aux,comment_att,'The SARIn Rx2 L1B Power waveforms is a fully calibrated, FF, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SARIn processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_l1b_echo'', in order to best fit the waveform into 2 bytes');
                
        phase_diff_meas_ku_l1b_echo_name = 'phase_diff_meas_ku_l1b_echo';
        id_aux = netcdf.defVar(ncid,phase_diff_meas_ku_l1b_echo_name,double_type,[ns_dimension nb_dimension]);
        netcdf.putAtt(ncid,id_aux,long_name_att,'SARIn Phase Difference');
        netcdf.putAtt(ncid,id_aux,units_att,rad_units);
        netcdf.putAtt(ncid,id_aux,comment_att,'Phase difference array, computed from the complex echoes on the 2 Rx channels. Values are between -pi and +pi radians, units are micro radians');
        
        coherence_meas_ku_l1b_echo_name = 'coherence_meas_ku_l1b_echo';
        id_aux = netcdf.defVar(ncid,coherence_meas_ku_l1b_echo_name,double_type,[ns_dimension nb_dimension]);
        netcdf.putAtt(ncid,id_aux,long_name_att,'SARIn Coherence');
        netcdf.putAtt(ncid,id_aux,units_att,number_units);
        netcdf.putAtt(ncid,id_aux,comment_att,'Coherence array , computed from the complex echoes on the 2 Rx channels');     
    end
end

if cnf.include_wfms_aligned
%EM: 14.04.2016
i2q2_meas_ku_wdcorr_l1b_echo_name = 'i2q2_meas_ku_wdcorr_l1b_echo';
id_aux = netcdf.defVar(ncid,i2q2_meas_ku_wdcorr_l1b_echo_name,int32_type,[ns_dimension nb_dimension]);
% switch netcdf_type
%     case 'netcdf4'
%         id_aux = netcdf.defVar(ncid,i2q2_meas_ku_wdcorr_l1b_echo_name,uint16_type,[ns_dimension nb_dimension]);
%     case 'netcdf3'
%         id_aux = netcdf.defVar(ncid,i2q2_meas_ku_wdcorr_l1b_echo_name,int32_type,[ns_dimension nb_dimension]);
% end
netcdf.putAtt(ncid,id_aux,long_name_att,'SAR	Power Echo waveform aligned surfaces: scaled	0-65535 (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'The SAR L1B Power waveforms (aligned w.r.t reference surface) is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_l1b_echo'', in order to best fit the waveform into 2 bytes. The waveforms have been algined taking into account window delay correction w.r.t first surface.');
end

stack_mask_range_bin_l1b_echo_name = 'stack_mask_range_bin_l1b_echo';
id_aux = netcdf.defVar(ncid,stack_mask_range_bin_l1b_echo_name,int16_type,[nl_dimension nb_dimension]);
%netcdf.defVarFill(ncid,id_aux,false,0);
netcdf.putAtt(ncid,id_aux,long_name_att,'Range bin stack mask ');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,cnf.zp_fact_range);
netcdf.putAtt(ncid,id_aux,comment_att,'The zero-mask applied to the stack before multilooking. Each element of the mask refers to a look in the stack and indicates the index of the first sample set to zero. The first num_looks_start_stop elements of the mask are valid, while the remaining ones are filled with 0.');

zero_padding_l1b_echo_name = 'zero_padding_l1b_echo';
id_aux = netcdf.defVar(ncid,zero_padding_l1b_echo_name,int16_type,[]); %scalar
netcdf.putAtt(ncid,id_aux,long_name_att,'Oversampling factor used in the range compression (FFT)');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'the ground processor can apply an oversampling factor, providing a waveform_sampling = nominal_sampling / range_oversampling_factor. Note that the altimeter range resolution is fixed and given by the chirp bandwidth');

%EM: 22.09.2016
%{
if strcmp(mode,'SIN') && strcmp(mission,'CR2') && strcmp(processing_mode_cnf,'SIN') 
    i2q2_meas_ku_l1b_echo_2_name = 'i2q2_meas_ku_l1b_echo_2';
    id_aux = netcdf.defVar(ncid,i2q2_meas_ku_l1b_echo_2_name,int32_type,[ns_dimension nb_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'SAR Power Echo waveform 2: scaled	0-65535 (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The SAR L1B Power waveform for Channel 2 (SARin). It is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_l1b_echo'', in order to best fit the waveform into 2 bytes');
    
    waveform_scale_factor_l1b_echo_2_name = 'waveform_scale_factor_l1b_echo_2';
    id_aux = netcdf.defVar(ncid,waveform_scale_factor_l1b_echo_2_name,float_type,nb_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Echo Scale Factor Channel 2, to convert from [0-65535]	to Power at	antenna	flange');
    netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The L1B waveform for Channel 2 (SARin) scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = i2q2_meas_ku_l1b_echo(ku_rec, Ns) * waveform_scale_factor_l1b_echo(ku_rec)');
    
    coherence_meas_ku_l1b_echo_name = 'coherence_meas_ku_l1b_echo';
    id_aux = netcdf.defVar(ncid,coherence_meas_ku_l1b_echo_name,int16_type,[ns_dimension nb_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Interferometric Coherence');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Coherence between the two interferometric channels operating in SARIn, both channels are fully calibrated.');
    
    phase_difference_meas_ku_l1b_echo_name = 'phase_difference_meas_ku_l1b_echo';
    id_aux = netcdf.defVar(ncid,phase_difference_meas_ku_l1b_echo_name,int32_type,[ns_dimension nb_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Phase difference SARIn  (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,rad_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
    netcdf.putAtt(ncid,id_aux,comment_att,'Interferometric Phase multilooked obtained from the pair of channels operating in SARIn');
   
    
    if cnf.include_wfms_aligned
        
        i2q2_meas_ku_wdcorr_l1b_echo_2_name = 'i2q2_meas_ku_wdcorr_l1b_echo_2';
        id_aux = netcdf.defVar(ncid,i2q2_meas_ku_wdcorr_l1b_echo_2_name,int32_type,[ns_dimension nb_dimension]);
        netcdf.putAtt(ncid,id_aux,long_name_att,'SAR Power Echo waveform channel 2 aligned surfaces: scaled	0-65535 (Ku-band)');
        netcdf.putAtt(ncid,id_aux,units_att,number_units);
        netcdf.putAtt(ncid,id_aux,comment_att,'The SAR L1B Power waveform channel 2 (aligned w.r.t given reference surface). It is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_l1b_echo'', in order to best fit the waveform into 2 bytes. The waveforms have been algined taking into account window delay correction w.r.t first surface.');
        
        coherence_meas_ku_wdcorr_l1b_echo_name = 'coherence_meas_ku_wdcorr_l1b_echo';
        id_aux = netcdf.defVar(ncid,coherence_meas_ku_wdcorr_l1b_echo_name,int16_type,[ns_dimension nb_dimension]);
        netcdf.putAtt(ncid,id_aux,long_name_att,'Interferometric Coherence (aligned surfaces)');
        netcdf.putAtt(ncid,id_aux,units_att,number_units);
        netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
        netcdf.putAtt(ncid,id_aux,comment_att,'Coherence between the two interferometric channels operating in SARIn (aligned w.r.t given reference surface), both channels are fully calibrated.');
        
        phase_difference_meas_ku_wdcorr_l1b_echo_name = 'phase_difference_meas_ku_wdcorr_l1b_echo';
        id_aux = netcdf.defVar(ncid,phase_difference_meas_ku_wdcorr_l1b_echo_name,int32_type,nb_dimension);
        netcdf.putAtt(ncid,id_aux,long_name_att,'Phase difference SARIn aligned surfaces  (Ku-band)');
        netcdf.putAtt(ncid,id_aux,units_att,rad_units);
        netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
        netcdf.putAtt(ncid,id_aux,comment_att,'Interferometric Phase multilooked obtained from the pair of channels operating in SARIn (aligned w.r.t given reference surface)');

    end
           
end
%}
%------------- N. Geophysical Corrections variables -----------------------

dry_tropo_correction_name = 'dry_tropo_correction_l1b_echo';
id_aux = netcdf.defVar(ncid,dry_tropo_correction_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Dry Tropospheric Correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

wet_tropo_correction_name = 'wet_tropo_correction_l1b_echo';
id_aux = netcdf.defVar(ncid,wet_tropo_correction_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Wet Tropospheric correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

inverse_baro_correction_name = 'inverse_baro_correction_l1b_echo';
id_aux = netcdf.defVar(ncid,inverse_baro_correction_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Inverse Barometric Correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

Dynamic_atmospheric_correction_name = 'Dynamic_atmospheric_correction_l1b_echo';
id_aux = netcdf.defVar(ncid,Dynamic_atmospheric_correction_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Dynamic Atmospheric Correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

GIM_iono_correction_name = 'GIM_iono_correction_l1b_echo';
id_aux = netcdf.defVar(ncid,GIM_iono_correction_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'GIM Ionospheric Correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

model_iono_correction_name = 'model_iono_correction_l1b_echo';
id_aux = netcdf.defVar(ncid,model_iono_correction_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Model Ionospheric Correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

ocean_equilibrium_tide_name = 'ocean_equilibrium_tide_l1b_echo';
id_aux = netcdf.defVar(ncid,ocean_equilibrium_tide_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Ocean Equilibrium Tide');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

long_period_tide_height_name = 'long_period_tide_l1b_echo';
id_aux = netcdf.defVar(ncid,long_period_tide_height_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Long Period Ocean Tide');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

ocean_loading_tide_name = 'ocean_loading_tide_l1b_echo';
id_aux = netcdf.defVar(ncid,ocean_loading_tide_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Ocean Loading Tide');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

solid_earth_tide_name = 'solid_earth_tide_l1b_echo';
id_aux = netcdf.defVar(ncid,solid_earth_tide_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Solid Earth Tide');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

geocentric_polar_tide_name = 'geocentric_polar_tide_l1b_echo';
id_aux = netcdf.defVar(ncid,geocentric_polar_tide_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Geocentric Polar Tide');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

%------------- O. Processing parameters used ------------------------------


%added EM 04.10.2016
% %------------- P. ACDC variables ----------------------------------------
%{
if ACDC_application_cnf
    %-------------- Waveform -----------------------------------------------
    %multilooked ACDC waveform
    i2q2_meas_ku_ACDC_echo_name = 'i2q2_meas_ku_ACDC_echo';
    id_aux = netcdf.defVar(ncid,i2q2_meas_ku_ACDC_echo_name,int32_type,[ns_dimension nb_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC SAR Power Echo waveform: scaled 0-65535 (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The ACDC SAR L1B Power waveform. It is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_ACDC_echo'', in order to best fit the waveform into 2 bytes');
    %fitted ACDC waveform
    i2q2_fit_ku_ACDC_echo_name = 'i2q2_fit_ku_ACDC_echo';
    id_aux = netcdf.defVar(ncid,i2q2_fit_ku_ACDC_echo_name,int32_type,[ns_dimension nb_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC SAR Fitted Power Echo waveform: scaled 0-65535 (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The fitted ACDC SAR L1B Power waveform. It is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_ACDC_echo'', in order to best fit the waveform into 2 bytes');
    %fitted L1B SAR waveform
    i2q2_fit_ku_echo_name = 'i2q2_fit_ku_echo';
    id_aux = netcdf.defVar(ncid,i2q2_fit_ku_echo_name,int32_type,[ns_dimension nb_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'SAR Fitted Power Echo waveform: scaled 0-65535 (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The fitted SAR L1B Power waveform. It is a fully calibrated, high resolution, multilooked waveform. It includes: (a) all calibrations, which have been applied at L1A, (b) SAR processor configuration according to the L1B processing flags, (c) final scaling, given in the variable ''waveform_scale_factor_ACDC_echo'', in order to best fit the waveform into 2 bytes');
        
    
    %-------------- Scaling -----------------------------------------------
    %scaling ACDC waveform
    waveform_scale_factor_ACDC_echo_name = 'waveform_scale_factor_ACDC_echo';
    id_aux = netcdf.defVar(ncid,waveform_scale_factor_ACDC_echo_name,float_type,nb_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Echo Scale Factor, to convert from [0-65535]	to Power at	antenna	flange');
    netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The ACDC waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = i2q2_meas_ku_ACDC_echo(ku_rec, Ns) * waveform_scale_factor_ACDC_echo(ku_rec)');
    %scaling fitted ACDC waveform
    waveform_scale_factor_fit_ACDC_echo_name = 'waveform_scale_factor_fit_ACDC_echo';
    id_aux = netcdf.defVar(ncid,waveform_scale_factor_fit_ACDC_echo_name,float_type,nb_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC fit Scale Factor, to convert from [0-65535]	to Power at	antenna	flange');
    netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The ACDC fitted waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = i2q2_meas_ku_ACDC_echo(ku_rec, Ns) * waveform_scale_factor_ACDC_echo(ku_rec)');    
    %scaling fitted ACDC waveform
    waveform_scale_factor_fit_echo_name = 'waveform_scale_factor_fit_echo';
    id_aux = netcdf.defVar(ncid,waveform_scale_factor_fit_echo_name,float_type,nb_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'SAR fit Scale Factor, to convert from [0-65535]	to Power at	antenna	flange');
    netcdf.putAtt(ncid,id_aux,units_att,W_per_count_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'The SAR fitted waveform scaling factor, computed in order to best fit each waveform within 2 bytes. The scaling, needed to convert the L1B waveform into Watt, is applied as follows: power_waveform_watt(ku_rec, Ns) = i2q2_meas_ku_ACDC_echo(ku_rec, Ns) * waveform_scale_factor_ACDC_echo(ku_rec)');
    %-------------- Range index -------------------------------------------
    range_index_ACDC_echo_name = 'range_index_ACDC_echo';
    id_aux = netcdf.defVar(ncid,range_index_ACDC_echo_name,int32_type,[ns_dimension nb_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC range index (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'The ACDC range index for ACDC SAR Power waveform.');
    
    %-------------- Geophysical retrievals --------------------------------
    %******************* retracked range **********************************
    retracked_range_ACDC_20_ku_name = 'retracked_range_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,retracked_range_ACDC_20_ku_name,int32_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Retracked range for Ku band (2-step analytical retracker)');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC corrected range by the retracker offset, the reference range includes instrumental corrections (already the USO frequency drift and the internal/instrument corrections).');
    %******************* epoch ********************************************
    epoch_ACDC_20_ku_name = 'epoch_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,epoch_ACDC_20_ku_name,int32_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Epoch (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'Estimated ACDC epoch in seconds w.r.t center of the window (window delay is given to the center of the window) using the 2-step analytical retracker. This corresponds to zero-padded sample value.');
    %******************** SWH  ********************************************
    swh_ACDC_20_ku_name = 'swh_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,swh_ACDC_20_ku_name,int16_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Significant waveheight (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC Fitted significant waveheight. This corresponds to 4 times the fitted standard deviation of the surface height.');
    swh_ACDC_PRE_20_ku_name = 'swh_ACDC_PRE_20_ku';
    id_aux = netcdf.defVar(ncid,swh_ACDC_PRE_20_ku_name,int16_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Pre-retracker ACDC Significant waveheight (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'Pre-retracker ACDC Fitted significant waveheight. This corresponds to 4 times the fitted standard deviation of the surface height.');
    %******************** SSH  ********************************************
    ssh_ACDC_20_ku_name = 'ssh_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,ssh_ACDC_20_ku_name,int32_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Surface height (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC Surface heigth above the elliposid of reference and extracted using the orbital height and the corrected range (retracked range).');
    ssh_ACDC_PRE_20_ku_name = 'ssh_ACDC_PRE_20_ku';
    id_aux = netcdf.defVar(ncid,ssh_ACDC_PRE_20_ku_name,int32_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Pre-retracker ACDC Surface height (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'Pre-retracker ACDC Surface heigth above the elliposid of reference and extracted using the orbital height and the corrected range (retracked range).');
    %******************** SIGMA0 ******************************************    
    sig0_ACDC_20_ku_name = 'sig0_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,sig0_ACDC_20_ku_name,int16_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,32767);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Backscattering coefficient (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC Backscattering coefficient extracted as the fitted peak power once corrected by the sigma0 scaling factor.');
    
    sig0_ACDC_PRE_20_ku_name = 'sig0_ACDC_PRE_20_ku';
    id_aux = netcdf.defVar(ncid,sig0_ACDC_PRE_20_ku_name,int16_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,32767);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Pre-retracker ACDC Backscattering coefficient (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'Pre-retracker ACDC Backscattering coefficient extracted as the fitted peak power once corrected by the sigma0 scaling factor.');
    
    %********************* Pu: fitted peak power **************************    
    Pu_analytical_20_ku_name = 'Pu_analytical_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,Pu_analytical_20_ku_name,int16_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Fitted peak power (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC Peak power of the fitted waveform.');
    %********************** Pearson Correlation Coefficient ***************
    Pearson_corr_ACDC_20_ku_name = 'Pearson_corr_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,Pearson_corr_ACDC_20_ku_name,int16_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Pearson coefficient (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,percent_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC Pearson correlation coefficient as percentage indicating the goodness of fitting between the real ACDC waveform and the fitted one.');

    Pearson_corr_ACDC_PRE_20_ku_name = 'Pearson_corr_ACDC_PRE_20_ku';
    id_aux = netcdf.defVar(ncid,Pearson_corr_ACDC_PRE_20_ku_name,int16_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Pre-retracker ACDC Pearson coefficient (Ku-band)');
    netcdf.putAtt(ncid,id_aux,units_att,percent_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
    netcdf.putAtt(ncid,id_aux,comment_att,'Pre-retracker ACDC Pearson correlation coefficient as percentage indicating the goodness of fitting between the real ACDC waveform and the fitted one.');    
    
    %************************ Fitting flag ********************************    
    Flag_fitting_ACDC_20_ku_name = 'Flag_fitting_ACDC_20_ku';
    id_aux = netcdf.defVar(ncid,Flag_fitting_ACDC_20_ku_name,int8_type, nb_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,127);
    netcdf.putAtt(ncid,id_aux,long_name_att,'ACDC Exit flag fitting');
    netcdf.putAtt(ncid,id_aux,flag_values_att,'1,2,3,4,0,-1,-2');
    netcdf.putAtt(ncid,id_aux,flag_desc_att,'1: Function converged to a solution x; 2: Change in x was less than the specified tolerance; 3: Change in the residual was less than the specified tolerance; 4: Magnitude of search direction was smaller than the specified tolerance; 0: Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations; -1: Output function terminated the algorithm; -2: Problem is infeasible: the bounds lb and ub are inconsistent');
    netcdf.putAtt(ncid,id_aux,comment_att,'ACDC Flag on the reason the LSE solver stopped.');
      
end
%}

%----------  Global Attributes definition -----------------------------------
%---- attributes inherited from Sentinel-3 product description-------------
id_aux = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,id_aux,'creation_time',date_creation);
netcdf.putAtt(ncid,id_aux,'Conventions',netcdf_v4_format);
%{
netcdf.putAtt(ncid,id_aux,'mission_name',mission);
netcdf.putAtt(ncid,id_aux,'altimeter_sensor_name',altimeter_sensor_name);
netcdf.putAtt(ncid,id_aux,'gnss_sensor_name',gnss_sensor_name);
netcdf.putAtt(ncid,id_aux,'doris_sensor_name',doris_sensor_name);
netcdf.putAtt(ncid,id_aux,'acq_station_name',acq_station_name);
netcdf.putAtt(ncid,id_aux,'doris_sensor_name',acq_station_name);
netcdf.putAtt(ncid,id_aux,'first_meas_time',first_meas_time);
netcdf.putAtt(ncid,id_aux,'last_meas_time',last_meas_time);
netcdf.putAtt(ncid,id_aux,'xref_altimeter_level0',xref_altimeter_level0);
netcdf.putAtt(ncid,id_aux,'xref_altimeter_orbit',xref_altimeter_orbit);
netcdf.putAtt(ncid,id_aux,'xref_doris_USO',xref_doris_USO);
netcdf.putAtt(ncid,id_aux,'xref_altimeter_ltm_sar_cal1',xref_altimeter_ltm_sar_cal1);
netcdf.putAtt(ncid,id_aux,'xref_altimeter_ltm_ku_cal2',xref_altimeter_ltm_ku_cal2);
netcdf.putAtt(ncid,id_aux,'xref_altimeter_ltm_c_cal2',xref_altimeter_ltm_c_cal2);
netcdf.putAtt(ncid,id_aux,'xref_altimeter_characterisation',xref_altimeter_characterisation);
netcdf.putAtt(ncid,id_aux,'semi_major_ellipsoid_axis',semi_major_ellipsoid_axis);
netcdf.putAtt(ncid,id_aux,'ellipsoid_flattening',ellipsoid_flattening);

%--------------- add the attributes related to intermediate product--------
netcdf.putAtt(ncid,id_aux,'orbit_phase_code',orbit_phase_code);
netcdf.putAtt(ncid,id_aux,'orbit_cycle_num',orbit_cycle_num);
netcdf.putAtt(ncid,id_aux,'orbit_REL_Orbit',orbit_REL_Orbit);
netcdf.putAtt(ncid,id_aux,'orbit_ABS_Orbit_Start',orbit_ABS_Orbit_Start);
netcdf.putAtt(ncid,id_aux,'orbit_Rel_Time_ASC_Node_Start',orbit_Rel_Time_ASC_Node_Start);
netcdf.putAtt(ncid,id_aux,'orbit_ABS_Orbit_Stop',orbit_ABS_Orbit_Stop);
netcdf.putAtt(ncid,id_aux,'orbit_Rel_Time_ASC_Node_Stop',orbit_Rel_Time_ASC_Node_Stop);
netcdf.putAtt(ncid,id_aux,'orbit_Equator_Cross_Time',orbit_Equator_Cross_Time);
netcdf.putAtt(ncid,id_aux,'orbit_Equator_Cross_Long',orbit_Equator_Cross_Long);
netcdf.putAtt(ncid,id_aux,'orbit_Ascending_Flag',orbit_Ascending_Flag);
%}
netcdf.endDef(ncid);

var_id = netcdf.inqVarID(ncid,'zero_padding_l1b_echo');
netcdf.putVar(ncid,var_id,cnf.zp_fact_range);

netcdf.close(ncid);
% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);
