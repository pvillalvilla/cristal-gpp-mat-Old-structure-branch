%% FINAL DATATION AND GEOLOCATION ALGORITHM
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% ---------------------------------------------------------
% Objective: Manage all the data along the SAR-Ku chain. 
% Author:    Roger Escol�       / isardSAT
%            Albert Garcia      / isardSAT
% Reviewer:  Monica Roca        / isardSAT
% Last rev.: Albert Garcia        / isardSAT (10/06/2019)
%
% v1.1 2019/06/10	- Change datation method to include the SWST (commented) 
%					- Added Geolocation
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L1A,L1AP] = final_datation (L1A,L1AP, cnf, chd, cst, ORBIT)

    
    %% 1. Datation Correction
    %The propagation delay computed with the window delay
    delta_prop_delay = L1A.win_delay / 2;

    %% 2. Centre of Gravity to antenna correction
    %PRI computed with the real altimeter period.
%     pri_sar_pre_dat = L1A.pri_sar * chd.pri_T0_unit_conv * chd.T0_sar_pre_dat;
%     L1A.bri_sar = N_pulses_burst * L1A.pri_sar;
    
    %The CoG to antenna correction is added. From the calling to the Attitude Selection algorithm, we
    %have pitch/yaw_pre_dat(0:N_total_bursts_sar_ku_isp(i_ISP)-1) and from the OSV Selection, the satellite's velocity.
    x_cog_ant_aux  = chd.x_cog_ant * cos(L1A.pitch_sar) + chd.z_cog_ant * sin(L1A.pitch_sar);
    x_cog_ant_corr = x_cog_ant_aux * cos(L1A.yaw_sar);
    t_cog_ant_corr = - x_cog_ant_corr / L1A.x_vel_sat_sar;
    
    %% NEED TO ADD ROLL EFFECT
    
    
    %Now the datation is computed again
%     burst_prop_delay = bri_sar_pre_dat * (burst_sar_isp(i_burst)-1); %this is the propagation of the time stamp along
    burst_prop_delay = 0;
    t_rx = L1A.time_rx_1st + chd.tx1_sar + chd.pulse_length/2 + burst_prop_delay;
    for i_pulse = 1:chd.N_pulses_burst
        L1AP.time_sar_ku_pre_dat_pulse(i_pulse) =...
            t_rx + i_pulse* L1A.pri_sar ...
            - delta_prop_delay + t_cog_ant_corr;
    end
    
    %% 3. Datation averaging
     %L1A.time = mean(L1AP.time_sar_ku_pre_dat_pulse(1+chd.N_cal_pulses_burst+chd.N_c_pulses_burst:chd.N_pulses_burst));
%      L1A.time = L1A.time_rx_1st - ((chd.SWST+L1A.ambiguity_order_sar_isp*L1A.pri_sar) - (chd.N_pulses_burst-1)*L1A.pri_sar)/2;

    %% 4. GEOLOCATION
    N_bursts_interpol = 8;
    
    start_burst = max(1,L1A.burst-N_bursts_interpol);
    end_burst = min(start_burst+N_bursts_interpol,length(ORBIT.time));
    L1A.x_sar_sat       = spline (ORBIT.time(start_burst:end_burst), ORBIT.com_position_vector(1,start_burst:end_burst), L1A.time);
    L1A.y_sar_sat       = spline (ORBIT.time(start_burst:end_burst), ORBIT.com_position_vector(2,start_burst:end_burst), L1A.time);
    L1A.z_sar_sat       = spline (ORBIT.time(start_burst:end_burst), ORBIT.com_position_vector(3,start_burst:end_burst), L1A.time);
    L1A.x_vel_sat_sar   = spline (ORBIT.time(start_burst:end_burst), ORBIT.com_velocity_vector(1,start_burst:end_burst), L1A.time);
    L1A.alt_rate_sar_sat= spline (ORBIT.time(start_burst:end_burst), ORBIT.alt_rate_sar_sat(start_burst:end_burst), L1A.time);

    L1A.y_vel_sat_sar   = spline (ORBIT.time(start_burst:end_burst), ORBIT.com_velocity_vector(2,start_burst:end_burst), L1A.time);
    L1A.z_vel_sat_sar   = spline (ORBIT.time(start_burst:end_burst), ORBIT.com_velocity_vector(3,start_burst:end_burst), L1A.time);
    L1A.roll_sar        = spline (ORBIT.time(start_burst:end_burst), ORBIT.roll(start_burst:end_burst), L1A.time);
    L1A.pitch_sar       = spline (ORBIT.time(start_burst:end_burst), ORBIT.pitch(start_burst:end_burst), L1A.time);
    L1A.yaw_sar         = spline (ORBIT.time(start_burst:end_burst), ORBIT.yaw(start_burst:end_burst), L1A.time);
 
 

end
