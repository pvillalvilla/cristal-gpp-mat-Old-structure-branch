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
function [files] = create_NetCDF_L1ASS(files, N_bursts,N_surf, cnf, chd, cst)


date_creation = datestr(now, '_yyyymmddTHHMMSS_');

L0_date_find = strfind(files.filename_L0,'_201');
L0_mode_find = strfind(files.filename_L0, 'PIC_SIRS___');

files.filename_L1AS = strcat(files.outputPath,'PICE_SIRS_', ...
                                files.filename_L0((L0_mode_find(1)+11):(L0_mode_find(end)+15)), '_1AS',...
                                files.filename_L0((L0_date_find(1)+1):(L0_date_find(1)+30)),...
                                '_isd','.nc');
%{
files.filename_L1AS = strcat(files.outputPath,'SR_1_SRA__A_',...
                                files.filename_L0(end-37:end-3),...
                                '_isd','.nc');
%}
files.ncid_L1AS = netcdf.create(files.filename_L1AS,'CLOBBER');

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


nb_dimension = netcdf.defDim(files.ncid_L1AS,'nb',N_bursts); % Number of burst/records
nsurf_dimension = netcdf.defDim(files.ncid_L1AS,'nss',N_surf); % Number of surfaces
nchd_dim = netcdf.defDim(files.ncid_L1AS,'nchd',1); % Number of surfaces
np_dimension = netcdf.defDim(files.ncid_L1AS,'np',chd.N_pulses_burst); % Number of pulses
ns_dimension = netcdf.defDim(files.ncid_L1AS,'ns',chd.N_samples_sar); % Number of samples
space_3D_dimension = netcdf.defDim(files.ncid_L1AS,'space_3D',3);

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



%% PACKING L1ASS


%----------A. Time variables ----------------------------------------------

time_l1a_echo_name = 'time_l1a_echo';
id_aux      = netcdf.defVar(files.ncid_L1AS,time_l1a_echo_name,double_type,nb_dimension);
            netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Burst Timestamp at bouncing');
            netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,seconds_units);

%------------- B. Looks related variables -----------------------------
i_meas_l1a_echo_name = 'rx1_complex_waveforms_i_samples';
id_aux = netcdf.defVar(files.ncid_L1AS,i_meas_l1a_echo_name,double_type,[ns_dimension np_dimension nb_dimension]);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Chain 1 echoes, i measurements ');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,number_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'Chain 1 echoes of a burst, I values (64*256 samples) in the time domain.');

q_meas_l1a_echo_name = 'rx1_complex_waveforms_q_samples';
id_aux = netcdf.defVar(files.ncid_L1AS,q_meas_l1a_echo_name,double_type,[ns_dimension np_dimension nb_dimension]);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Chain 1 echoes, q measurements ');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,number_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'Chain 1 echoes of a burst, Q values (64*256 samples) in the time domain.');


% i_scale_factor_name = 'rx1_i_scale_factor';
% id_aux = netcdf.defVar(files.ncid_L1AS,i_scale_factor_name,double_type, nb_dimension);
% netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Chain 1 I scale factor, to convert I samples from [-127,+127] to amplitude at antenna flange');
% netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,sqrtW_per_count_units);
% netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'Chain 1 The i_samples scaling factor, computed in order to best fit the i_samples within 1 byte. The scaling, needed to convert the i_samples into sqrt(watt), is applied as follows: i_samples_sqr_watt(Nb,Np, Ns) = i_samples(Nb,Np, Ns) * i_scale_factor(Nb,Np) ');
% 
% q_scale_factor_name = 'rx1_q_scale_factor';
% id_aux = netcdf.defVar(files.ncid_L1AS,q_scale_factor_name,double_type, nb_dimension);
% netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Chain 1 Q scale factor, to convert I samples from [-127,+127] to amplitude at antenna flange');
% netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,sqrtW_per_count_units);
% netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'Chain 1 The q_samples scaling factor, computed in order to best fit the q_samples within 1 byte. The scaling, needed to convert the q_samples into sqrt(watt), is applied as follows: q_samples_sqr_watt(Nb,Np, Ns) = q_samples(Nb,Np, Ns) * q_scale_factor(Nb,Np) ');
% 


if strcmp(cnf.processing_mode,'SIN')
    
    i_meas_l1a_echo_name = 'rx2_complex_waveforms_i_samples';
    id_aux = netcdf.defVar(files.ncid_L1AS,i_meas_l1a_echo_name,double_type,[ns_dimension np_dimension nb_dimension]);
    netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Chain 2 echoes, i measurements ');
    netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,number_units);
    netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'Chain 2 echoes of a burst, I values (64*256 samples) in the time domain.');
    
    q_meas_l1a_echo_name = 'rx2_complex_waveforms_q_samples';
    id_aux = netcdf.defVar(files.ncid_L1AS,q_meas_l1a_echo_name,double_type,[ns_dimension np_dimension nb_dimension]);
    netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Chain 2 echoes, q measurements ');
    netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,number_units);
    netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'Chain 2 echoes of a burst, Q values (64*256 samples) in the time domain.');
    
%     i_scale_factor_name = 'rx2_i_scale_factor';
%     id_aux = netcdf.defVar(files.ncid_L1AS,i_scale_factor_name,double_type, nb_dimension);
%     netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Chain 2 I scale factor, to convert I samples from [-127,+127] to amplitude at antenna flange');
%     netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,sqrtW_per_count_units);
%     netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'Chain 2 The i_samples scaling factor, computed in order to best fit the i_samples within 1 byte. The scaling, needed to convert the i_samples into sqrt(watt), is applied as follows: i_samples_sqr_watt(Nb, Np, Ns) = i_samples(Nb, Np, Ns) * i_scale_factor(Nb,Np) ');
%     
%     q_scale_factor_name = 'rx2_q_scale_factor';
%     id_aux = netcdf.defVar(files.ncid_L1AS,q_scale_factor_name,double_type, nb_dimension);
%     netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,' Chain 2 Q scale factor, to convert I samples from [-127,+127] to amplitude at antenna flange');
%     netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,sqrtW_per_count_units);
%     netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'Chain 2 The q_samples scaling factor, computed in order to best fit the q_samples within 1 byte. The scaling, needed to convert the q_samples into sqrt(watt), is applied as follows: q_samples_sqr_watt(Nb, Np, Ns) = q_samples(Nb, Np, Ns) * q_scale_factor(Nb,Np) ');

end

%----------C. Position/Velocity variables ------------------------------
x_pos_l1a_echo_name = 'x_pos_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,x_pos_l1a_echo_name,double_type,nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,18446744073709551616);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Satellite-x component');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,meters_units);

y_pos_l1a_echo_name = 'y_pos_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,y_pos_l1a_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,18446744073709551616);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Satellite altitude-y component');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,meters_units);

z_pos_l1a_echo_name = 'z_pos_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,z_pos_l1a_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,18446744073709551616);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Satellite altitude-z component');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,meters_units);

x_vel_l1a_echo_name = 'x_vel_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,x_vel_l1a_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,18446744073709551616);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Satellite velocity-x component');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,scale_factor_att,1.e-4);

y_vel_l1a_echo_name = 'y_vel_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,y_vel_l1a_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,18446744073709551616);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Satellite velocity-y component');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,scale_factor_att,1.e-4);

z_vel_l1a_echo_name = 'z_vel_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,z_vel_l1a_echo_name,double_type, nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,18446744073709551616);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Satellite velocity-z component');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,scale_factor_att,1.e-4);


alt_l1a_echo_name = 'alt_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,alt_l1a_echo_name,int32_type,nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,214748364);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'altitude of satellite');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,meters_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(files.ncid_L1AS,id_aux,add_offset_att,700000);
netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'Altitude of the satellite Centre of Mass');


orb_alt_rate_l1a_echo_name = 'alt_rate_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,orb_alt_rate_l1a_echo_name,int16_type,nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,32767);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'orbital altitude rate');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,meters_per_second_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(files.ncid_L1AS,id_aux,add_offset_att,0);
netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'Instantaneous altitude rate at the Centre of Mass');

roll_sral_mispointing_l1a_echo_name = 'roll_mispointing_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,roll_sral_mispointing_l1a_echo_name,int16_type,nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,214748364);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Mispointing angle - roll ');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,degrees_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(files.ncid_L1AS,id_aux,add_offset_att,0.);
% netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'value for the closest in time from the burst time tag, given in the nadir pointing reference frame');

pitch_sral_mispointing_l1a_echo_name = 'pitch_mispointing_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,pitch_sral_mispointing_l1a_echo_name,int16_type,nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,214748364);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Mispointing angle - pitch ');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,degrees_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(files.ncid_L1AS,id_aux,add_offset_att,0.);
% netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'value for the closest in time from the burst time tag, given in the nadir pointing reference frame');

yaw_sral_mispointing_l1a_echo_name = 'yaw_mispointing_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,yaw_sral_mispointing_l1a_echo_name,int16_type,nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,214748364);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Mispointing angle - yaw ');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,degrees_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(files.ncid_L1AS,id_aux,add_offset_att,0.);
% netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'value for the closest in time from the burst time tag, given in the nadir pointing reference frame');



%---------- I. Altimeter range and Corrections ----------------------------
range_l1a_echo_name = 'range_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,range_l1a_echo_name,int32_type, nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,2147483647);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Corrected range for Ku band');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,meters_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(files.ncid_L1AS,id_aux,add_offset_att,700000);
netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'Reference range corrected for USO frequency drift and internal path correction');

%------------- M. BURST Timming -----------------------------
pri_l1a_name = 'pri_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,pri_l1a_name,double_type, nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,2147483647);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Pulse repetition interval');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,seconds_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,scale_factor_att,1e-12);
netcdf.putAtt(files.ncid_L1AS,id_aux,add_offset_att,0.);
% netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'Pulse Repetition Interval');
bri_l1a_name = 'bri_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,bri_l1a_name,double_type, nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,2147483647);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Burst repetition interval');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,seconds_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,scale_factor_att,1e-6);
netcdf.putAtt(files.ncid_L1AS,id_aux,add_offset_att,0.);

ambiguity_order_l1a_name = 'ambiguity_order_l1a_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,ambiguity_order_l1a_name,int16_type, nb_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,2147483647);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'l');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,number_units);
netcdf.putAtt(files.ncid_L1AS,id_aux,scale_factor_att,1.);
netcdf.putAtt(files.ncid_L1AS,id_aux,add_offset_att,0.);
netcdf.putAtt(files.ncid_L1AS,id_aux,comment_att,'The Pulse Ambiguity Rank is the number of pulses that are transmitted between the transmission and the reception of each pulse ');

time_l1ass_echo_name = 'time_l1ass_echo';
id_aux      = netcdf.defVar(files.ncid_L1AS,time_l1ass_echo_name,double_type,nsurf_dimension);
            netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Surface Timestamp');
            netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,seconds_units);


x_pos_l1ass_echo_name = 'x_pos_l1ass_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,x_pos_l1ass_echo_name,double_type,nsurf_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,18446744073709551616);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Surface-x component');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,meters_units);

y_pos_l1ass_echo_name = 'y_pos_l1ass_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,y_pos_l1ass_echo_name,double_type, nsurf_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,18446744073709551616);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Surface altitude-y component');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,meters_units);

z_pos_l1ass_echo_name = 'z_pos_l1ass_echo';
id_aux = netcdf.defVar(files.ncid_L1AS,z_pos_l1ass_echo_name,double_type, nsurf_dimension);
%             netcdf.defVarFill(files.ncid_L1AS,id_aux,false,18446744073709551616);
netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Surface altitude-z component');
netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,meters_units);

%------------- N. CHD Parameters -----------------------------

freq_name = 'freq';
id_aux      = netcdf.defVar(files.ncid_L1AS,freq_name,double_type,nchd_dim);
            netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Carrier Frequency');
            netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,Hz_units);
bw_name = 'bw';
id_aux      = netcdf.defVar(files.ncid_L1AS,bw_name,double_type,nchd_dim);
            netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Tx chirp bandwidth');
            netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,Hz_units);
pulse_length_name = 'pulse_length';
id_aux      = netcdf.defVar(files.ncid_L1AS,pulse_length_name,double_type,nchd_dim);
            netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Tx chirp length');
            netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,seconds_units);
chirp_slope_name = 'chirp_slope';
id_aux      = netcdf.defVar(files.ncid_L1AS,chirp_slope_name,double_type,nchd_dim);
            netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'Tx chirp slope');
            netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,'Hz/s');
alt_freq_name = 'alt_freq';
id_aux      = netcdf.defVar(files.ncid_L1AS,alt_freq_name,double_type,nchd_dim);
            netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'raw echoes sampling frequency');
            netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,Hz_units);
N_bursts_cycle_name = 'N_bursts_cycle';
id_aux      = netcdf.defVar(files.ncid_L1AS,N_bursts_cycle_name,double_type,nchd_dim);
            netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'tracking cycle length');
            netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,number_units);
N_samples_sar_name = 'N_samples_sar';
id_aux      = netcdf.defVar(files.ncid_L1AS,N_samples_sar_name,double_type,nchd_dim);
            netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'number of range samples');
            netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,number_units);
N_pulses_burst_name = 'N_pulses_burst';
id_aux      = netcdf.defVar(files.ncid_L1AS,N_pulses_burst_name,double_type,nchd_dim);
            netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'number of pulses in a burst');
            netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,number_units);
antenna_beamwidth_along_track_name = 'antenna_beamwidth_along_track';
id_aux      = netcdf.defVar(files.ncid_L1AS,antenna_beamwidth_along_track_name,double_type,nchd_dim);
            netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'along track -3dB aperture');
            netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,number_units);
antenna_beamwidth_across_track_name = 'antenna_beamwidth_across_track';
id_aux      = netcdf.defVar(files.ncid_L1AS,antenna_beamwidth_across_track_name,double_type,nchd_dim);
            netcdf.putAtt(files.ncid_L1AS,id_aux,long_name_att,'across track -3dB aperture');
            netcdf.putAtt(files.ncid_L1AS,id_aux,units_att,number_units);


    



netcdf.endDef(files.ncid_L1AS);


% netcdf.close(files.ncid_L1AS);
% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);
