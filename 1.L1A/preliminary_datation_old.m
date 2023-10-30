%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
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

function [L1A, L1AP] = preliminary_datation (L1A,L1AP, cnf, chd, cst)

%% 0.Initialisation and Memory allocation

tp                          = zeros(1,chd.N_pulses_burst);
time_sar_ku_pre_dat_pulse   = zeros(1,chd.N_pulses_burst);


%%ARESYS aproach:
%  - the Rx time of the first sample of the first echo of the burst t_rx_first(k) is read from the NetCDF data_record_time variable;
%  - the bouncing time t_bouncing(k) is the mean between t_rx_first(k) and the Tx time of the first sample of the last echo of the burst t_tx_last(k) and can be computed as:
%    t_bouncing(k) = mean(t_tx_last(k),t_rx_first(k)) =  t_rx_first(k) - ((SWST+RANK*PRI) - (Nechoes-1)*PRI)/2
%    where:
%                  Sºhe time when the receiving window is opened with respect to the RANK*PRI time;
%                  Nechoes is the number of echoes in each burst (64 for CLOSED BURST acquisition mode and 65 for OPEN BURST);
%  
% Probably, in order to compute the bouncing time you need the SWST value:
% for CLOSED BURST mode SWST = 2.4745193075612004e-05 s;
% for OPEN BURST mode SWST = 2.6037973915238429e-05 s;





%% 1.Time stamp delay

% burst_prop_delay = double(chd.N_pulses_burst .* L1A.pri_sar * double(L1A.burst-1)); %this is the propagation of the time stamp along the radar cycle %DO NOT NEEDED AS WE HAVE A DIFFERENT TIME FOR EACH BURST
burst_prop_delay = 0;
t_rx = L1A.time_rx_1st + chd.tx1_sar + chd.pulse_length/2 + burst_prop_delay;

%% 2.Shift to each received pulse
for i_pulse = 1:chd.N_pulses_burst
    tp(i_pulse) = t_rx -...
        double(L1A.ambiguity_order_sar_isp-i_pulse+1) * L1A.pri_sar;
end

%% 3.Locate at the surface
% delta_prop_delay                        = cst.mean_sat_alt/cst.c;
% L1AP.time_sar_ku_pre_dat_pulse    = tp - delta_prop_delay;

% JPLZ: using current sat alt instead of a mean one, and adding prop delay instead of substracting it
delta_prop_delay                        = L1A.alt_sar_sat/cst.c;
L1AP.time_sar_ku_pre_dat_pulse    = tp + delta_prop_delay;

%% 4.Datation averaging
% L1A.time = mean(L1AP.time_sar_ku_pre_dat_pulse(1+chd.N_cal_pulses_burst+chd.N_c_pulses_burst:chd.N_pulses_burst));

%L1A.time = L1A.time_rx_1st - ((chd.SWST+L1A.ambiguity_order_sar_isp*L1A.pri_sar) - (chd.N_pulses_burst-1)*L1A.pri_sar)/2;
L1A.time = mean(L1AP.time_sar_ku_pre_dat_pulse);




end