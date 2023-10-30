%% FINAL DATATION AND GEOLOCATION ALGORITHM
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L.
% --------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [FINAL_DAT_OUT] = final_datation (T0_pre_dat, WIN_DELAY_OUT, ATTITUDE, ORBIT, ISP, cnf, chd, cst)

%% 1. Centre of Gravity to antenna correction

if(cnf.flag_cog_antenna_datation_correction == 1)
    x_cog_ant1_corr = chd.x_cog_ant * ( cos(ATTITUDE.pitch_pre_dat) + cos(ATTITUDE.yaw_pre_dat) ) - chd.y_cog_ant * sin(ATTITUDE.yaw_pre_dat) + ...
        chd.z_cog_ant * sin(ATTITUDE.pitch_pre_dat);
    x_cog_ant2_corr = chd.x_cog_ant2 * ( cos(ATTITUDE.pitch_pre_dat) + cos(ATTITUDE.yaw_pre_dat) ) + chd.y_cog_ant2 * sin(ATTITUDE.yaw_pre_dat) + ...
        chd.z_cog_ant2 * sin(ATTITUDE.pitch_pre_dat);
    
    t_CoG2ant_1_corr = x_cog_ant1_corr / vel_along;
    t_CoG2ant_2_corr = ( (x_cog_ant1_corr + x_cog_ant2_corr) / 2 ) / vel_along;
    
end

%% 1. Datation Correction
% PRI is computed with the real altimeter clock period
PRI_c = ISP.pri * chd.pri_T0_unit_conv * T0_pre_dat;

if( strcmp(chd.meas_mode, 'OPEN_BURST') )
    Namb = ISP.ambiguity_order;
    bri = PRI_c * chd.N_pulses_burst;

elseif( strcmp(chd.meas_mode, 'CLOSED_BURST') )
    Namb = 0;
    bri = chd.bri;
    
end

% Band selection
if( strcmp(ISP.band, 'Ku') )
    t_cog_ant_corr = t_CoG2ant_1_corr;
    
elseif( strcmp(ISP.band, 'Ku2') )
    t_cog_ant_corr = t_CoG2ant_2_corr;
    
elseif( strcmp(ISP.band, 'Ka') )
    t_cog_ant_corr = t_CoG2ant_1_corr;
    
end

% Update window delay with altitude rate
win_delay = WIN_DELAY_OUT.win_delay + t_cog_ant_corr * ORBIT.alt_rate*2/cst.c;

%The propagation delay computed with the window delay
delta_prop_delay = win_delay / 2;

% burst_prop_delay = 0;
t_tx =ISP.time0 + chd.tx1 + chd.pulse_length/2;

for i_burst = 0:chd.N_bursts_rc-1
    for i_pulse = 0:chd.N_pulses_burst-1
        time_pulse(i_burst+1,i_pulse+1) = t_tx - ...
            (Namb - i_pulse) * PRI_c + i_burst*bri +...
            delta_prop_delay + t_cog_ant_corr;
    end
end

%% Averaging
t_surf = mean(time_pulse,2);

%% Geolocation
% Using auxiliary function with time_lrm as inpu

%% Output
FINAL_DAT_OUT.t_surf = t_surf;
FINAL_DAT_OUT.band = ISP.band;
FINAL_DAT_OUT.chd.meas_mode = chd.meas_mode;

end
