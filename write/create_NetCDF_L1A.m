%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: Create L1A
%
% Calling: 
% INPUTs:
%
%
% OUTPUTs:
%
%
% ----------------------------------------------------------
% Author:    Albert Garcia  / isardSAT
%            Eduard Makhoul / isardSAT
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function [files] = create_NetCDF_L1A(files, N_bursts, cnf, chd, cst)


date_creation = datestr(now, '_yyyymmddTHHMMSS_');

L0_date_find = strfind(files.filename_L0,'_201');
L0_date_find = strfind(files.filename_L0,'_202'); %JPLZ: date for IRIS OB chronogram simulations is 202X

L0_mode_find = strfind(files.filename_L0, 'S3NGT_SIRS___');

files.filename_L1A = strcat(files.outputPath,'S3NGT_SIRS_', ...
                                files.filename_L0((L0_mode_find(1)+11):(L0_mode_find(end)+15)), '_1A_',...
                                files.filename_L0((L0_date_find(1)+1):(L0_date_find(1)+30)),...
                                '_isd','.nc');
%{
files.filename_L1A = strcat(files.outputPath,'SR_1_SRA__A_',...
                                files.filename_L0(end-37:end-3),...
                                '_isd','.nc');
%}
ncid = netcdf.create(files.filename_L1A,'CLOBBER');

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


nb_dimension = netcdf.defDim(ncid,'nb',N_bursts); % Number of burst/records
np_dimension = netcdf.defDim(ncid,'np',chd.N_pulses_burst); % Number of pulses


if(cnf.onboard_reversion_flag)
    ns_dimension = netcdf.defDim(ncid,'ns',chd.N_samples_sar*2); % Number of samples
else
    ns_dimension = netcdf.defDim(ncid,'ns',chd.N_samples_sar); % Number of samples
end
space_3D_dimension = netcdf.defDim(ncid,'space_3D',3);
n1_dimension= netcdf.defDim(ncid,'n1',1); % a dimension 1 variable

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
rad_units = 'radian';
percent_units='percent';


int8_type = 'NC_BYTE';
uint8_type = 'NC_BYTE';
int16_type = 'NC_SHORT';
uint16_type = 'NC_SHORT';
int32_type = 'NC_INT';
uint32_type = 'NC_INT';
uint64_type = 'uint64';
float_type = 'NC_FLOAT';
double_type= 'NC_DOUBLE';



%% PACKING L1A


%----------A. Time variables ----------------------------------------------
l1_mode_id_name = 'l1_mode_id';
id_aux      = netcdf.defVar(ncid,l1_mode_id_name,uint8_type,nb_dimension);
            netcdf.putAtt(ncid,id_aux,long_name_att,'L1 mode ID');
            netcdf.putAtt(ncid,id_aux,comment_att,'L1 Mode Identifier. Each L1A record type has a unique ID, as follows: 0 SAR, 1 SARIn');



time_l1a_echo_name = 'time_l1a_echo';
id_aux      = netcdf.defVar(ncid,time_l1a_echo_name,double_type,nb_dimension);
            netcdf.putAtt(ncid,id_aux,long_name_att,'UTC ');
            netcdf.putAtt(ncid,id_aux,units_att,seconds_units);


seq_count_l1a_echo_name = 'seq_count_l1a_echo';
id_aux      = netcdf.defVar(ncid,seq_count_l1a_echo_name,int16_type, nb_dimension);
            netcdf.putAtt(ncid,id_aux,'FillValue',65535);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Source Sequence Count: Same as Burst number ');
            netcdf.putAtt(ncid,id_aux,units_att,number_units);

UTC_day_l1a_echo_name = 'UTC_day_l1a_echo';
id_aux 		= netcdf.defVar(ncid,UTC_day_l1a_echo_name,double_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,1.844674407370960e+19);
			netcdf.putAtt(ncid,id_aux,long_name_att,'day UTC ');
			netcdf.putAtt(ncid,id_aux,units_att,day_units);


UTC_sec_l1a_echo_name = 'UTC_sec_l1a_echo';
id_aux 		= netcdf.defVar(ncid,UTC_sec_l1a_echo_name,double_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,1.844674407370960e+19);
			netcdf.putAtt(ncid,id_aux,long_name_att,'seconds in the day UTC ');
			netcdf.putAtt(ncid,id_aux,units_att,seconds_units);

uso_cor_l1a_echo_name = 'uso_cor_l1a_echo';
id_aux = netcdf.defVar(ncid,uso_cor_l1a_echo_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'USO frequency drift correction ');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
% netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');

            
            
burst_count_prod_l1a_echo_name = 'burst_count_prod_l1a_echo';
id_aux = netcdf.defVar(ncid,burst_count_prod_l1a_echo_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'bursts counter within the product ');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
% netcdf.putAtt(ncid,id_aux,comment_att,'range 1 to number of records in the product');            
            
cog_cor_l1a_echo_name = 'cog_cor_l1a_echo';
id_aux = netcdf.defVar(ncid,cog_cor_l1a_echo_name,int16_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'Z-Distance antenna-CoG correction ');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Distance  in the z-component between the centre of mass of the satellite and the altimeter antenna reference point');



%----------D. Position/Velocity variables ------------------------------
x_pos_l1a_echo_name = 'x_pos_l1a_echo';
id_aux = netcdf.defVar(ncid,x_pos_l1a_echo_name,double_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite-x component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

y_pos_l1a_echo_name = 'y_pos_l1a_echo';
id_aux = netcdf.defVar(ncid,y_pos_l1a_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite altitude-y component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

z_pos_l1a_echo_name = 'z_pos_l1a_echo';
id_aux = netcdf.defVar(ncid,z_pos_l1a_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite altitude-z component');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);

x_vel_l1a_echo_name = 'x_vel_l1a_echo';
id_aux = netcdf.defVar(ncid,x_vel_l1a_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite velocity-x component');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);

y_vel_l1a_echo_name = 'y_vel_l1a_echo';
id_aux = netcdf.defVar(ncid,y_vel_l1a_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite velocity-y component');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);

z_vel_l1a_echo_name = 'z_vel_l1a_echo';
id_aux = netcdf.defVar(ncid,z_vel_l1a_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
netcdf.putAtt(ncid,id_aux,long_name_att,'Satellite velocity-z component');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);


alt_l1a_echo_name = 'alt_l1a_echo';
id_aux = netcdf.defVar(ncid,alt_l1a_echo_name,int32_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'altitude of satellite');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Altitude of the satellite Centre of Mass');


orb_alt_rate_l1a_echo_name = 'alt_rate_l1a_echo';
id_aux = netcdf.defVar(ncid,orb_alt_rate_l1a_echo_name,int32_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
netcdf.putAtt(ncid,id_aux,long_name_att,'orbital altitude rate');
netcdf.putAtt(ncid,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Instantaneous altitude rate at the Centre of Mass');

roll_sral_mispointing_l1a_echo_name = 'roll_mispointing_l1a_echo';
id_aux = netcdf.defVar(ncid,roll_sral_mispointing_l1a_echo_name,int16_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'Mispointing angle - roll ');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
% netcdf.putAtt(ncid,id_aux,comment_att,'value for the closest in time from the burst time tag, given in the nadir pointing reference frame');

pitch_sral_mispointing_l1a_echo_name = 'pitch_mispointing_l1a_echo';
id_aux = netcdf.defVar(ncid,pitch_sral_mispointing_l1a_echo_name,int16_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'Mispointing angle - pitch ');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
% netcdf.putAtt(ncid,id_aux,comment_att,'value for the closest in time from the burst time tag, given in the nadir pointing reference frame');

yaw_sral_mispointing_l1a_echo_name = 'yaw_mispointing_l1a_echo';
id_aux = netcdf.defVar(ncid,yaw_sral_mispointing_l1a_echo_name,int16_type,nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'Mispointing angle - yaw ');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
% netcdf.putAtt(ncid,id_aux,comment_att,'value for the closest in time from the burst time tag, given in the nadir pointing reference frame');



%---------- G. H0, COR2 and AGC ----------------------------------------
h0_applied_l1a_echo_name = 'h0_applied_l1a_echo';
id_aux = netcdf.defVar(ncid,h0_applied_l1a_echo_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,4294967295);
netcdf.putAtt(ncid,id_aux,long_name_att,'Applied altitude command H0');
netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d64d1e9_units);
netcdf.putAtt(ncid,id_aux,comment_att,'From ISP. Applied altitude command H0');

cor2_applied_l1a_echo_name = 'cor2_applied_l1a_echo';
id_aux = netcdf.defVar(ncid,cor2_applied_l1a_echo_name,int16_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
netcdf.putAtt(ncid,id_aux,long_name_att,'Applied altitude command COR2');
netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d1024d1e9_units);
netcdf.putAtt(ncid,id_aux,comment_att,'From ISP. Applied altitude variation');

% agc_l1a_echo_name = 'agc_l1a_echo';
% id_aux = netcdf.defVar(ncid,agc_l1a_echo_name,int32_type, nb_dimension);
% %             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
% netcdf.putAtt(ncid,id_aux,long_name_att,'corrected AGC');
% netcdf.putAtt(ncid,id_aux,units_att,dB_units);
% netcdf.putAtt(ncid,id_aux,scale_factor_att,0.01);
% netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
% netcdf.putAtt(ncid,id_aux,comment_att,'AGC corrected for instrumental errors (calibration)');


burst_count_cycle_l1a_echo_name = 'burst_count_cycle_l1a_echo';
id_aux = netcdf.defVar(ncid,burst_count_cycle_l1a_echo_name,int8_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'bursts counter within the cycle ');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Cyclic counter from 1 to 4');            
           

% sig0_cal_l1a_echo_name = 'sig0_cal_l1a_echo';
% id_aux = netcdf.defVar(ncid,sig0_cal_l1a_echo_name,int32_type, nb_dimension);
% %             netcdf.defVarFill(ncid,id_aux,false,18446744073709551616);
% netcdf.putAtt(ncid,id_aux,long_name_att,'internal calibration correction on Sigma0');
% netcdf.putAtt(ncid,id_aux,units_att,dB_units);
% netcdf.putAtt(ncid,id_aux,scale_factor_att,0.01);
% netcdf.putAtt(ncid,id_aux,add_offset_att,0.);



%---------- I. Altimeter range and Corrections ----------------------------
range_l1a_echo_name = 'range_l1a_echo';
id_aux = netcdf.defVar(ncid,range_l1a_echo_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Corrected range for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Reference range corrected for USO frequency drift and internal path correction');

int_path_cor_l1a_echo_name = 'int_path_cor_l1a_echo';
id_aux = netcdf.defVar(ncid,int_path_cor_l1a_echo_name,int32_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Internal path correction for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
% netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement');


%------------- M. BURST Timming -----------------------------
pri_l1a_name = 'pri_l1a_echo';
id_aux = netcdf.defVar(ncid,pri_l1a_name,double_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Pulse repetition interval');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-12);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
% netcdf.putAtt(ncid,id_aux,comment_att,'Pulse Repetition Interval');
bri_l1a_name = 'bri_l1a_echo';
id_aux = netcdf.defVar(ncid,bri_l1a_name,double_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Burst repetition interval');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);

ambiguity_order_l1a_name = 'ambiguity_order_l1a_echo';
id_aux = netcdf.defVar(ncid,ambiguity_order_l1a_name,int16_type, nb_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'l');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'The Pulse Ambiguity Rank is the number of pulses that are transmitted between the transmission and the reception of each pulse ');


%------------- M. Looks related variables -----------------------------
i_meas_l1a_echo_name = 'rx1_complex_waveforms_i_samples';
id_aux = netcdf.defVar(ncid,i_meas_l1a_echo_name,double_type,[ns_dimension np_dimension nb_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'Chain 1 echoes, i measurements ');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Chain 1 echoes of a burst, I values (64*256 samples) in the time domain. The pulses are not corrected for AGC');

q_meas_l1a_echo_name = 'rx1_complex_waveforms_q_samples';
id_aux = netcdf.defVar(ncid,q_meas_l1a_echo_name,double_type,[ns_dimension np_dimension nb_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'Chain 1 echoes, q measurements ');
netcdf.putAtt(ncid,id_aux,units_att,number_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Chain 1 echoes of a burst, Q values (64*256 samples) in the time domain. The pulses are not corrected for AGC');


% i_scale_factor_name = 'rx1_i_scale_factor';
% id_aux = netcdf.defVar(ncid,i_scale_factor_name,double_type, nb_dimension);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Chain 1 I scale factor, to convert I samples from [-127,+127] to amplitude at antenna flange');
% netcdf.putAtt(ncid,id_aux,units_att,sqrtW_per_count_units);
% netcdf.putAtt(ncid,id_aux,comment_att,'Chain 1 The i_samples scaling factor, computed in order to best fit the i_samples within 1 byte. The scaling, needed to convert the i_samples into sqrt(watt), is applied as follows: i_samples_sqr_watt(Nb,Np, Ns) = i_samples(Nb,Np, Ns) * i_scale_factor(Nb,Np) ');
% 
% q_scale_factor_name = 'rx1_q_scale_factor';
% id_aux = netcdf.defVar(ncid,q_scale_factor_name,double_type, nb_dimension);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Chain 1 Q scale factor, to convert I samples from [-127,+127] to amplitude at antenna flange');
% netcdf.putAtt(ncid,id_aux,units_att,sqrtW_per_count_units);
% netcdf.putAtt(ncid,id_aux,comment_att,'Chain 1 The q_samples scaling factor, computed in order to best fit the q_samples within 1 byte. The scaling, needed to convert the q_samples into sqrt(watt), is applied as follows: q_samples_sqr_watt(Nb,Np, Ns) = q_samples(Nb,Np, Ns) * q_scale_factor(Nb,Np) ');
% 


if strcmp(cnf.processing_mode,'SIN')
    
    i_meas_l1a_echo_name = 'rx2_complex_waveforms_i_samples';
    id_aux = netcdf.defVar(ncid,i_meas_l1a_echo_name,double_type,[ns_dimension np_dimension nb_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Chain 2 echoes, i measurements ');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'Chain 2 echoes of a burst, I values (64*256 samples) in the time domain. The pulses are not corrected for AGC');
    
    q_meas_l1a_echo_name = 'rx2_complex_waveforms_q_samples';
    id_aux = netcdf.defVar(ncid,q_meas_l1a_echo_name,double_type,[ns_dimension np_dimension nb_dimension]);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Chain 2 echoes, q measurements ');
    netcdf.putAtt(ncid,id_aux,units_att,number_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'Chain 2 echoes of a burst, Q values (64*256 samples) in the time domain. The pulses are not corrected for AGC');
    
%     i_scale_factor_name = 'rx2_i_scale_factor';
%     id_aux = netcdf.defVar(ncid,i_scale_factor_name,double_type, nb_dimension);
%     netcdf.putAtt(ncid,id_aux,long_name_att,'Chain 2 I scale factor, to convert I samples from [-127,+127] to amplitude at antenna flange');
%     netcdf.putAtt(ncid,id_aux,units_att,sqrtW_per_count_units);
%     netcdf.putAtt(ncid,id_aux,comment_att,'Chain 2 The i_samples scaling factor, computed in order to best fit the i_samples within 1 byte. The scaling, needed to convert the i_samples into sqrt(watt), is applied as follows: i_samples_sqr_watt(Nb, Np, Ns) = i_samples(Nb, Np, Ns) * i_scale_factor(Nb,Np) ');
%     
%     q_scale_factor_name = 'rx2_q_scale_factor';
%     id_aux = netcdf.defVar(ncid,q_scale_factor_name,double_type, nb_dimension);
%     netcdf.putAtt(ncid,id_aux,long_name_att,' Chain 2 Q scale factor, to convert I samples from [-127,+127] to amplitude at antenna flange');
%     netcdf.putAtt(ncid,id_aux,units_att,sqrtW_per_count_units);
%     netcdf.putAtt(ncid,id_aux,comment_att,'Chain 2 The q_samples scaling factor, computed in order to best fit the q_samples within 1 byte. The scaling, needed to convert the q_samples into sqrt(watt), is applied as follows: q_samples_sqr_watt(Nb, Np, Ns) = q_samples(Nb, Np, Ns) * q_scale_factor(Nb,Np) ');

    
end

%----OB CAL MASK VARIABLES---

OB_cal_flag_name='OB_cal_flag';
id_aux = netcdf.defVar(ncid,OB_cal_flag_name,double_type,nb_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'OB_cal_flag');
netcdf.putAtt(ncid,id_aux,units_att,'');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'Flag indicating the presence of CAL pulses in OB processing');

OB_cal_mask_name='OB_cal_mask';
id_aux = netcdf.defVar(ncid,OB_cal_mask_name,double_type,[n1_dimension np_dimension nb_dimension]);
netcdf.putAtt(ncid,id_aux,long_name_att,'OB_cal_flag');
netcdf.putAtt(ncid,id_aux,units_att,'');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'Mask indicating the position of the CAL pulses inside the burst (0=CAL, 1=no CAL)');


% var_id = netcdf.inqVarID(files.ncid_L1A,'OB_cal_flag');
% netcdf.putVar(files.ncid_L1A,var_id,i_burst,L1A.OB_cal_flag);
% 
% var_id = netcdf.inqVarID(files.ncid_L1A,'OB_cal_mask');
% netcdf.putVar(files.ncid_L1A,var_id,i_burst,L1A.OB_cal_mask);
% netcdf.endDef(ncid);
% 
% 
% alt_l1a_echo_name = 'alt_l1a_echo';
% id_aux = netcdf.defVar(ncid,alt_l1a_echo_name,int32_type,nb_dimension);
% %             netcdf.defVarFill(ncid,id_aux,false,214748364);
% netcdf.putAtt(ncid,id_aux,long_name_att,'altitude of satellite');
% netcdf.putAtt(ncid,id_aux,units_att,meters_units);
% netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
% netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
% netcdf.putAtt(ncid,id_aux,comment_att,'Altitude of the satellite Centre of Mass');
% 
% h0_applied_l1a_echo_name = 'h0_applied_l1a_echo';
% id_aux = netcdf.defVar(ncid,h0_applied_l1a_echo_name,int32_type, nb_dimension);
% %             netcdf.defVarFill(ncid,id_aux,false,4294967295);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Applied altitude command H0');
% netcdf.putAtt(ncid,id_aux,units_att,seconds_3dot125d64d1e9_units);
% netcdf.putAtt(ncid,id_aux,comment_att,'From ISP. Applied altitude command H0');

netcdf.close(ncid);
% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);
