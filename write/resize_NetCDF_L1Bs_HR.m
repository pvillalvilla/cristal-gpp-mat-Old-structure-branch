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

function resize_NetCDF_L1Bs_HR(files,N_total_surfaces,max_beams_x_stack,cnf,chd,cst)

ncid = netcdf.open(files.filename_L1Bs,'WRITE');
ncid_2remove=netcdf.open(files.filename_netCDF_2remove_Bs,'NOWRITE');

ns_dimension = chd.N_samples_sar*cnf.zp_fact_range;
nl_dimension = max_beams_x_stack;
np_dimension = chd.N_pulses_burst;
space_3D_dimension = 3;

%%----------0. Define groups ----------------------------------------------
% Define the subgroups of the short netCDF that we will keep
global_ncid = netcdf.inqNcid(ncid,'global');
global_ku_ncid = netcdf.inqNcid(global_ncid,'ku');
global_ka_ncid = netcdf.inqNcid(global_ncid,'ka');

data_ncid = netcdf.inqNcid(ncid,'data');
data_ku_ncid = netcdf.inqNcid(data_ncid,'ku');
data_ka_ncid = netcdf.inqNcid(data_ncid,'ka');

if strcmp(chd.band,'Ku')
    global_currentband_ncid=global_ku_ncid;
    data_currentband_ncid=data_ku_ncid;
elseif strcmp(chd.band,'Ka')
    global_currentband_ncid=global_ka_ncid;
    data_currentband_ncid=data_ka_ncid;
end

% Define the subgroups of the long netCDF that we will remove
global_ncid_2remove = netcdf.inqNcid(ncid_2remove,'global');
global_ku_ncid_2remove = netcdf.inqNcid(global_ncid_2remove,'ku');
global_ka_ncid_2remove = netcdf.inqNcid(global_ncid_2remove,'ka');

data_ncid_2remove = netcdf.inqNcid(ncid_2remove,'data');
data_ku_ncid_2remove = netcdf.inqNcid(data_ncid_2remove,'ku');
data_ka_ncid_2remove = netcdf.inqNcid(data_ncid_2remove,'ka');

if strcmp(chd.band,'Ku')
    global_currentband_ncid_2remove=global_ku_ncid_2remove;
    data_currentband_ncid_2remove=data_ku_ncid_2remove;
elseif strcmp(chd.band,'Ka')
    global_currentband_ncid_2remove=global_ka_ncid_2remove;
    data_currentband_ncid_2remove=data_ka_ncid_2remove;
end

%% PACKING L1Bs
%----------A. Time & counters variables ----------------------------------------------
var_id = netcdf.inqVarID(ncid,'looks');
aux=netcdf.getVar(ncid_2remove,var_id,0,nl_dimension);
netcdf.putVar(ncid,var_id,0, nl_dimension,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'samples_ov');
aux=netcdf.getVar(ncid_2remove,var_id,0,ns_dimension);
netcdf.putVar(ncid,var_id,0, ns_dimension,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'space_3d');
aux=netcdf.getVar(ncid_2remove,var_id,0,space_3D_dimension);
netcdf.putVar(ncid,var_id,0, space_3D_dimension,aux); clear aux;

var_id = netcdf.inqVarID(ncid,'pulses');
aux=netcdf.getVar(ncid_2remove,var_id,0,np_dimension);
netcdf.putVar(ncid,var_id,0, np_dimension,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'l1b_record_counter');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0, N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'time');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'time_tai');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0, N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'tm_source_sequence_counter');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

%----------B. Orbit and attitude variables ------------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'altitude');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'altitude_rate');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'latitude');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'longitude');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

if strcmp(chd.band,'Ku')
    var_id = netcdf.inqVarID(data_ku_ncid,'off_nadir_roll_angle_ant1');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    var_id = netcdf.inqVarID(data_ku_ncid,'off_nadir_pitch_angle_ant1');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    var_id = netcdf.inqVarID(data_ku_ncid,'off_nadir_yaw_angle_ant1');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'off_nadir_roll_angle_ant2');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    var_id = netcdf.inqVarID(data_ku_ncid,'off_nadir_pitch_angle_ant2');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    var_id = netcdf.inqVarID(data_ku_ncid,'off_nadir_yaw_angle_ant2');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

if strcmp(chd.band,'Ka')
    var_id = netcdf.inqVarID(data_ka_ncid,'off_nadir_roll_angle');
    aux=netcdf.getVar(data_ka_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ka_ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    var_id = netcdf.inqVarID(data_ka_ncid,'off_nadir_pitch_angle');
    aux=netcdf.getVar(data_ka_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ka_ncid,var_id,0,N_total_surfaces,aux); clear aux;
    
    var_id = netcdf.inqVarID(data_ka_ncid,'off_nadir_yaw_angle');
    aux=netcdf.getVar(data_ka_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ka_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

var_id = netcdf.inqVarID(global_currentband_ncid,'orbit_type_flag');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,0,1);
netcdf.putVar(global_currentband_ncid,var_id,0,1,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'position_vector');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0],[3 N_total_surfaces]);
netcdf.putVar(data_currentband_ncid,var_id,[0 0],[3 N_total_surfaces],aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'velocity_vector');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0],[3 N_total_surfaces]);
netcdf.putVar(data_currentband_ncid,var_id,[0 0],[3 N_total_surfaces],aux); clear aux;

%----------C. Configuration and quality variables --------------------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'mcd_flags');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(global_currentband_ncid,'processing_configuration_flags');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,0,1);
netcdf.putVar(global_currentband_ncid,var_id,0,1,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'hr_mcd_flags');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(global_currentband_ncid,'hr_processing_configuration_flags');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,0,1);
netcdf.putVar(global_currentband_ncid,var_id,0,1,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'iris_instrument_configuration_flags');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'iris_mode_flag');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(global_currentband_ncid,'range_oversampling_factor');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,0,1);
netcdf.putVar(global_currentband_ncid,var_id,0,1,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'telemetry_type_flag');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

%----------------D. Altimeter range variables ---------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'range_cor_com_ant1');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'range_cor_com_ant2');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);    
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

if strcmp(chd.band,'Ka')
    var_id = netcdf.inqVarID(data_ka_ncid,'range_cor_com');
    aux=netcdf.getVar(data_ka_ncid_2remove,var_id,0,N_total_surfaces);    
    netcdf.putVar(data_ka_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

var_id = netcdf.inqVarID(data_currentband_ncid,'range_cor_doppler');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'range_cor_internal_delay_cal_rx1');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'range_cor_internal_delay_cal_rx2');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);    
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

var_id = netcdf.inqVarID(data_currentband_ncid,'range_cor_internal_delay_att_rx1');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'range_cor_internal_delay_cal_rx2');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);    
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

var_id = netcdf.inqVarID(data_currentband_ncid,'range_cor_uso');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'range_cor_reference_sample');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'tracker_range_calibrated');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'tracker_range_calibrated_gnss');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

%----------E. Altimeter power variables ------------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'altimeter_power_drift_rx1');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'altimeter_power_drift_rx2');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

var_id = netcdf.inqVarID(data_currentband_ncid,'attenuation_calibrated_rx1');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'attenuation_calibrated_rx2');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

var_id = netcdf.inqVarID(data_currentband_ncid,'power_scaling_to_antenna_rx1');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'power_scaling_to_antenna_rx2');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

var_id = netcdf.inqVarID(data_currentband_ncid,'variable_digital_gain');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'cal1_power_rx1');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'cal1_power_rx2');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

%----------F. Altimeter engineering variables --------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'altimeter_clock');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'hn_mean');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'pulse_repetition_interval');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'tm_cor2');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'tm_h0');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

%----------G. Altimeter characterization variables --------------------------
var_id = netcdf.inqVarID(global_currentband_ncid,'residual_fixed_digital_gain_raw');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,0,1);
netcdf.putVar(global_currentband_ncid,var_id,0,1,aux); clear aux;

var_id = netcdf.inqVarID(global_currentband_ncid,'residual_fixed_digital_gain_rmc');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,0,1);
netcdf.putVar(global_currentband_ncid,var_id,0,1,aux); clear aux;

var_id = netcdf.inqVarID(global_currentband_ncid,'range_cor_external_group_delay_rx1');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,0,1);
netcdf.putVar(global_currentband_ncid,var_id,0,1,aux); clear aux;

var_id = netcdf.inqVarID(global_currentband_ncid,'range_cor_external_group_delay_rx2');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,0,1);
netcdf.putVar(global_currentband_ncid,var_id,0,1,aux); clear aux;

var_id = netcdf.inqVarID(global_currentband_ncid,'g_scaling_rx1');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,0,1);
netcdf.putVar(global_currentband_ncid,var_id,0,1,aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(global_ku_ncid,'g_scaling_rx2');
    aux=netcdf.getVar(global_ku_ncid_2remove,var_id,0,1);
    netcdf.putVar(global_ku_ncid,var_id,0,1,aux); clear aux;
end

var_id = netcdf.inqVarID(global_currentband_ncid,'antenna_gain_ant1');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,0,1);
netcdf.putVar(global_currentband_ncid,var_id,0,1,aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(global_ku_ncid,'antenna_gain_ant2');
    aux=netcdf.getVar(global_ku_ncid_2remove,var_id,0,1);
    netcdf.putVar(global_ku_ncid,var_id,0,1,aux); clear aux;
end

if strcmp(chd.band,'Ka')
    var_id = netcdf.inqVarID(global_ka_ncid,'antenna_gain_ant');
    aux=netcdf.getVar(global_ka_ncid_2remove,var_id,0,1);
    netcdf.putVar(global_ka_ncid,var_id,0,1,aux); clear aux;
end

%------------------H. Flag variables --------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'manoeuvre_flag');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'surface_classification_flag');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

%------------------J. Waveform related variables --------------------------
var_id = netcdf.inqVarID(global_currentband_ncid,'burst_phase_array_cor_rx1');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,[0 0], [chd.N_pulses_burst 1]);
netcdf.putVar(global_currentband_ncid,var_id,[0 0],[chd.N_pulses_burst 1],aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(global_ku_ncid,'burst_phase_array_cor_rx2');
    aux=netcdf.getVar(global_ku_ncid_2remove,var_id,[0 0], [chd.N_pulses_burst 1]);
    netcdf.putVar(global_ku_ncid,var_id,[0 0],[chd.N_pulses_burst 1],aux); clear aux;
end

var_id = netcdf.inqVarID(global_currentband_ncid,'burst_power_array_cor_rx1');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,[0 0], [chd.N_pulses_burst 1]);
netcdf.putVar(global_currentband_ncid,var_id,[0 0],[chd.N_pulses_burst 1],aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(global_ku_ncid,'burst_power_array_cor_rx2');
    aux=netcdf.getVar(global_ku_ncid_2remove,var_id,[0 0], [chd.N_pulses_burst 1]);    
    netcdf.putVar(global_ku_ncid,var_id,[0 0],[chd.N_pulses_burst 1],aux); clear aux;
end

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(global_ku_ncid,'instr_ext_phase_cor');
    aux=netcdf.getVar(global_ku_ncid_2remove,var_id,0,1);    
    netcdf.putVar(global_ku_ncid,var_id,0,1,aux); clear aux;
end

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'instr_int_phase_cor');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);   
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'phase_slope_cor');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,0,N_total_surfaces);    
    netcdf.putVar(data_ku_ncid,var_id,0,N_total_surfaces,aux); clear aux;
end

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'interf_base_vector');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,[0 0],[3 N_total_surfaces]);    
    netcdf.putVar(data_ku_ncid,var_id,[0 0],[3 N_total_surfaces],aux); clear aux;
end

var_id = netcdf.inqVarID(data_currentband_ncid,'pnr_estimation_name');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'ptr_main_lobe_width_name');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'sig0_scaling_factor');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'snr_instr_estimation');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid,var_id,0,N_total_surfaces,aux); clear aux;

%----------------N. Look related variables --------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'look_counter');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0],[nl_dimension N_total_surfaces]); 
netcdf.putVar(data_currentband_ncid,var_id,[0 0],[nl_dimension N_total_surfaces],aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'look_time');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0],[nl_dimension N_total_surfaces]); 
netcdf.putVar(data_currentband_ncid,var_id,[0 0],[nl_dimension N_total_surfaces],aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'look_i_samples_rx1');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0 0],[ns_dimension nl_dimension N_total_surfaces]); 
netcdf.putVar(data_currentband_ncid,var_id,[0 0 0],[ns_dimension nl_dimension N_total_surfaces],aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'look_i_samples_rx2');
    aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0 0],[ns_dimension nl_dimension N_total_surfaces]); 
    netcdf.putVar(data_ku_ncid,var_id,[0 0 0],[ns_dimension nl_dimension N_total_surfaces],aux); clear aux;
end

var_id = netcdf.inqVarID(data_currentband_ncid,'look_q_samples_rx1');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0 0],[ns_dimension nl_dimension N_total_surfaces]); 
netcdf.putVar(data_currentband_ncid,var_id,[0 0 0],[ns_dimension nl_dimension N_total_surfaces],aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'look_q_samples_rx2');
    aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0 0],[ns_dimension nl_dimension N_total_surfaces]); 
    netcdf.putVar(data_ku_ncid,var_id,[0 0 0],[ns_dimension nl_dimension N_total_surfaces],aux); clear aux;
end

var_id = netcdf.inqVarID(data_currentband_ncid,'iq_scale_factor_rx1');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0],[nl_dimension N_total_surfaces]); 
netcdf.putVar(data_currentband_ncid,var_id,[0 0],[nl_dimension N_total_surfaces],aux); clear aux;

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'iq_scale_factor_rx2');
    aux=netcdf.getVar(data_ku_ncid_2remove,var_id,[0 0],[nl_dimension N_total_surfaces]); 
    netcdf.putVar(data_ku_ncid ,var_id,[0 0],[nl_dimension N_total_surfaces],aux); clear aux;
end

%----------------O. Look characterization variables --------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'look_index');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0],[nl_dimension N_total_surfaces]); 
netcdf.putVar(data_currentband_ncid ,var_id,[0 0],[nl_dimension N_total_surfaces],aux);

var_id = netcdf.inqVarID(data_currentband_ncid,'look_angle');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0],[nl_dimension N_total_surfaces]); 
netcdf.putVar(data_currentband_ncid ,var_id,[0 0],[nl_dimension N_total_surfaces],aux);

var_id = netcdf.inqVarID(data_currentband_ncid,'doppler_angle');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0],[nl_dimension N_total_surfaces]); 
netcdf.putVar(data_currentband_ncid ,var_id,[0 0],[nl_dimension N_total_surfaces],aux);

var_id = netcdf.inqVarID(data_currentband_ncid,'pointing_angle');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0],[nl_dimension N_total_surfaces]); 
netcdf.putVar(data_currentband_ncid ,var_id,[0 0],[nl_dimension N_total_surfaces],aux);

var_id = netcdf.inqVarID(data_currentband_ncid,'slant_range_correction_applied');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0],[nl_dimension N_total_surfaces]); 
netcdf.putVar(data_currentband_ncid ,var_id,[0 0],[nl_dimension N_total_surfaces],aux);

var_id = netcdf.inqVarID(data_currentband_ncid,'doppler_correction_applied');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0],[nl_dimension N_total_surfaces]); 
netcdf.putVar(data_currentband_ncid ,var_id,[0 0],[nl_dimension N_total_surfaces],aux);

var_id = netcdf.inqVarID(data_currentband_ncid,'stack_mask');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,[0 0],[nl_dimension N_total_surfaces]); 
netcdf.putVar(data_currentband_ncid ,var_id,[0 0],[nl_dimension N_total_surfaces],aux);

%----------R. Rate & collocation variables --------------------------
var_id = netcdf.inqVarID(global_currentband_ncid,'posting_rate');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,0,1);
netcdf.putVar(global_currentband_ncid,var_id,0,1,aux); clear aux;

var_id = netcdf.inqVarID(global_currentband_ncid,'ku_ka_collocation_flag');
aux=netcdf.getVar(global_currentband_ncid_2remove,var_id,0,1);
netcdf.putVar(global_currentband_ncid ,var_id,0,1,aux); clear aux;

%----------T. Thermistor temperature variables --------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'thr0_txrf_ku');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr1_txrf_ka');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr2_rxrf_ku1');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr3_rxrf_ku2');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr4_rxrf_ka');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr5_lo');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr6_sspa_ku');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr7_sspa_ka');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr8_txnum');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr9_rxnum_ku');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr10_rxnum_ka');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr11_dcdc_isps_1');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr12_dcdc_isps_2');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr13_formatter');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

var_id = netcdf.inqVarID(data_currentband_ncid,'thr14_sequencer');
aux=netcdf.getVar(data_currentband_ncid_2remove,var_id,0,N_total_surfaces);
netcdf.putVar(data_currentband_ncid ,var_id,0,N_total_surfaces,aux); clear aux;

netcdf.close(ncid_2remove);
netcdf.close(ncid);
% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);
% netcdf.close(files.ncid);