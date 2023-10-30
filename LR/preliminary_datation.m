%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Aligned with CRISTAL v4a_DPM_L1__20230116
% ---------------------------------------------------------
% Objective: The computation of the waveforms preliminary datation.
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
% v1.1 2019/06/10 Computation of the datation using the SWST provided by Aresys

function [PRE_DAT_OUT] = preliminary_datation (ISP, cnf, chd, cst)

%% 1.Time stamp delay
t_tx = ISP.time_lrm + chd.tx1_lrm + chd.pulse_length / 2; 

%% 2.Shift to each transmitted pulse
pri_nom_lrm = ISP.pri_lrm * chd.pri_T0_unit_conv * chd.T0_nom;

if( strcmp(chd.meas_mode, 'OPEN_BURST') ) 
    Namb = ISP.ambiguity_order_lrm;
    bri_lrm = pri_nom_lrm * chd.N_pulses_burst;
    
elseif( strcmp(chd.meas_mode, 'CLOSED_BURST') )
    Namb = 0;
    bri_lrm = chd.bri_cb;

end

for i_burst = 0:chd.N_bursts_rc-1
    for i_pulse = 0:chd.N_pulses_burst-1
        tp(i_burst+1,i_pulse+1) = t_tx - (Namb - i_pulse) * pri_nom_lrm + i_burst * bri_lrm;
    end
end

%% 3.Locate at the surface
delta_prop_delay = chd.mean_sat_alt / cst.c;
time_lrm_pre_dat_pulse = tp(:) + delta_prop_delay;

%% 4.Datation averaging
time_lrm_pre_dat = mean(time_lrm_pre_dat_pulse);

%% Output
PRE_DAT_OUT.time_lrm_pre_dat = time_lrm_pre_dat;
PRE_DAT_OUT.band = ISP.band;
PRE_DAT_OUT.chd.meas_mode = chd.meas_mode;

%% Auxiliary data selection
% Update auxiliary parameters calling different selection algorithms
% (section 3) with the new time value 'time_lrm_pre_dat'


end