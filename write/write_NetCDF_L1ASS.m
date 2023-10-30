%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the CODING & PACKING 
% algorithm for Level-1 ASS product for input in Fully Focused Chain
%
% ---------------------------------------------------------
% Objective: Pack variables and write the NETCDF
% 
% INPUTs : Workspace
% OUTPUTs: TM Structure as defined on isardSAT_JasonCS_DPM
%
% ----------------------------------------------------------
% Author:    Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT ()
% 
% Versions
% 1.0 
% 1.1 Updated time conversion for data between 2010 and 2016 (CR2)
% 2.0 Transformed to a function. Writting one record
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function write_NetCDF_L1ASS(files, L1A, i_burst, cnf, chd, cst)
i_burst=i_burst-1;
%% CODING L1A
%----------A. Time variables ----------------------------------------------


switch L1A.ProcessID
    case 'SAR'
%         switch L1A.ins_id
%           case
        l1_mode_id = 0;
%           case
    case 'SIN'
        %         switch L1A.ins_id
%           case
        l1_mode_id = 1;
%           case        
end

time_day = floor(double(L1A.time / cst.sec_in_day));
time_seconds = int64((L1A.time - double(time_day) * cst.sec_in_day) * 1e6);
seq_count_l1a_echo = int16(L1A.burst);
uso_corr=0; % TBD

burst_counter = int16(L1A.burst);
burst_counter_cycle=mod(burst_counter,4);

%----------B. Orbit and attitude variables ------------------------------
latitude = int32(L1A.lat_sar_sat * 1e7);
for i_lon = 1:size(L1A.lon_sar_sat,1)
    if (L1A.lon_sar_sat(i_lon) > 180)
        longitude = uint32((L1A.lon_sar_sat +360)* 1e7); 
    else
        longitude = uint32((L1A.lon_sar_sat )* 1e7); 
    end
end

com_altitude = int32((L1A.alt_sar_sat - 0.7e6) * 1e4);
com_altitude_rate = int32(L1A.alt_rate_sar_sat * 1e4);
com_velocity_vector = int32([L1A.x_vel_sat_sar * 1e4; L1A.y_vel_sat_sar * 1e4; L1A.z_vel_sat_sar * 1e4]);
satellite_mispointing = int32([L1A.pitch_sar * 180/cst.pi * 1e4; L1A.roll_sar * 180/cst.pi * 1e4; L1A.yaw_sar * 180/cst.pi * 1e4]);
pitch_bias = 0; % to CHD file
roll_bias = 0; % to CHD file
yaw_bias = 0; % to CHD file
% mispointing_bias = int32([pitch_bias * 180/cst.pi * 1e4, roll_bias * 180/cst.pi * 1e4, yaw_bias * 180/cst.pi * 1e4]);

%----------C. Configuration and quality variables -------------------------
l1_instrument_configuration = 0;
l1a_mcd = 0; % TBD


%----------D. Altimeter range variables -----------------------------------  
altimeter_range_calibrated = int32((L1A.win_delay * cst.c/2 - 0.7e6) * 1e4);
range_corr_internal_delay = 0;
z_cog_ant_corr=0; % TBDefined
range_corr_com = int16(-(z_cog_ant_corr) * 1e4);

%----------E. Altimeter power variables -----------------------------------
% attenuator_calibrated =int16(L1A.att_isp * 1e2);
altimeter_power_drift = 0; % int16(power_var_cal1 * 1e2);
power_corr_digital_processing = 0; %= int16((onboard_proc_sar)*1e2);
power_scaling_to_antenna = 0; %= int16(gain_corr_instr_sar * 1e2);

%----------F. Altimeter engineering variables -----------------------------
altimeter_clock = int32(1./chd.T0_nom - 3.95e8)* 1e9;
tm_h0 = int32(L1A.h0_comp_sar_isp);
tm_cor2 = int16(L1A.cor2_comp_sar_isp);
% cai = int32(L1A.cai_sar_isp(:,1)).';
% fai = int32(L1A.fai_sar_isp(:,1)).';
tm_pri = double(L1A.pri_sar*1e12);
tm_bri = double(L1A.bri_sar*1e6);
tm_ambiguity_rank = int16(L1A.ambiguity_order_sar_isp);
% tm_nimp = int16(L1A.nimp_sar_isp);

%----------L. Waveform related variables -----------------------------------
% tm_burst_num = int8(L1A.burst);
i_wfm_cal_corrected = real(L1A.wfm_cal_gain_corrected.');
q_wfm_cal_corrected = imag(L1A.wfm_cal_gain_corrected.');
% i_scale_factor = zeros(chd.N_pulses_burst);
% q_scale_factor = zeros(chd.N_pulses_burst);
% i_samples = int8(zeros(chd.N_samples_sar,chd.N_pulses_burst));
% q_samples = int8(zeros(chd.N_samples_sar,chd.N_pulses_burst));

% i_scale_factor = double(max(max(abs(i_wfm_cal_corrected(:,:))))*10^0.3 / (2^7-1));
% q_scale_factor = double(max(max(abs(q_wfm_cal_corrected(:,:))))*10^0.3 / (2^7-1));
% i_samples = int8(round(i_wfm_cal_corrected(:,:)*10^0.3./i_scale_factor));
% q_samples = int8(round(q_wfm_cal_corrected(:,:)*10^0.3./q_scale_factor));

if strcmp(cnf.processing_mode,'SIN')
    i_wfm_cal_corrected_2 = real(L1A.wfm_cal_gain_corrected_2.');
    q_wfm_cal_corrected_2 = imag(L1A.wfm_cal_gain_corrected_2.');
    
%     i_scale_factor_2 = double(max(max(abs(i_wfm_cal_corrected_2(:,:))))*10^0.3 / (2^7-1));
%     q_scale_factor_2 = double(max(max(abs(q_wfm_cal_corrected_2(:,:))))*10^0.3 / (2^7-1));
%     i_samples = int8(round(i_wfm_cal_corrected(:,:)*10^0.3./i_scale_factor));
%     q_samples = int8(round(q_wfm_cal_corrected(:,:)*10^0.3./q_scale_factor));
    
end



snr_estimation = 0; 



%% PACKING L1A

%----------A. Time variables ----------------------------------------------
var_id = netcdf.inqVarID(files.ncid_L1AS,'time_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,L1A.time);


%----------D. Position/Velocity variables ------------------------------
var_id = netcdf.inqVarID(files.ncid_L1AS,'x_pos_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,L1A.x_sar_sat);
var_id = netcdf.inqVarID(files.ncid_L1AS,'y_pos_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,L1A.y_sar_sat);
var_id = netcdf.inqVarID(files.ncid_L1AS,'z_pos_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,L1A.z_sar_sat);
var_id = netcdf.inqVarID(files.ncid_L1AS,'x_vel_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,com_velocity_vector(1));
var_id = netcdf.inqVarID(files.ncid_L1AS,'y_vel_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,com_velocity_vector(2));
var_id = netcdf.inqVarID(files.ncid_L1AS,'z_vel_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,com_velocity_vector(3));
var_id = netcdf.inqVarID(files.ncid_L1AS,'alt_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,com_altitude);
var_id = netcdf.inqVarID(files.ncid_L1AS,'alt_rate_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,com_altitude_rate);
var_id = netcdf.inqVarID(files.ncid_L1AS,'roll_mispointing_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,satellite_mispointing(1));
var_id = netcdf.inqVarID(files.ncid_L1AS,'pitch_mispointing_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,satellite_mispointing(2));
var_id = netcdf.inqVarID(files.ncid_L1AS,'yaw_mispointing_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,satellite_mispointing(3));

%---------- I. Altimeter range and Corrections ----------------------------
var_id = netcdf.inqVarID(files.ncid_L1AS,'range_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,altimeter_range_calibrated);

var_id = netcdf.inqVarID(files.ncid_L1AS,'ambiguity_order_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,tm_ambiguity_rank);

var_id = netcdf.inqVarID(files.ncid_L1AS,'bri_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,tm_bri);

var_id = netcdf.inqVarID(files.ncid_L1AS,'pri_l1a_echo');
netcdf.putVar(files.ncid_L1AS,var_id,i_burst,tm_pri);


%------------- M. Burst related variables -----------------------------
var_id = netcdf.inqVarID(files.ncid_L1AS,'rx1_complex_waveforms_i_samples');
netcdf.putVar(files.ncid_L1AS,var_id,[0 0 i_burst],[chd.N_samples_sar chd.N_pulses_burst 1], i_wfm_cal_corrected); %256 x 64
var_id = netcdf.inqVarID(files.ncid_L1AS,'rx1_complex_waveforms_q_samples');
netcdf.putVar(files.ncid_L1AS,var_id,[0 0 i_burst],[chd.N_samples_sar chd.N_pulses_burst 1], q_wfm_cal_corrected); %256 x 64

% var_id = netcdf.inqVarID(files.ncid_L1AS,'rx1_i_scale_factor');
% netcdf.putVar(files.ncid_L1AS,var_id,i_burst, i_scale_factor);
% var_id = netcdf.inqVarID(files.ncid_L1AS,'rx1_q_scale_factor');
% netcdf.putVar(files.ncid_L1AS,var_id, i_burst,q_scale_factor); 

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(files.ncid_L1AS,'rx2_complex_waveforms_i_samples');
    netcdf.putVar(files.ncid_L1AS,var_id,[0 0 i_burst],[chd.N_samples_sar chd.N_pulses_burst 1],i_wfm_cal_corrected_2);
    var_id = netcdf.inqVarID(files.ncid_L1AS,'rx2_complex_waveforms_q_samples');
    netcdf.putVar(files.ncid_L1AS,var_id,[0 0 i_burst],[chd.N_samples_sar chd.N_pulses_burst 1],q_wfm_cal_corrected_2);
    
%     var_id = netcdf.inqVarID(files.ncid_L1AS,'rx2_i_scale_factor');
%     netcdf.putVar(files.ncid_L1AS,var_id,i_burst, i_scale_factor_2);
%     var_id = netcdf.inqVarID(files.ncid_L1AS,'rx2_q_scale_factor');
%     netcdf.putVar(files.ncid_L1AS,var_id,i_burst, q_scale_factor_2);
    
end





