% PRELIMINARY WINDOW DELAY ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% JasonCS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: The computation of the waveforms preliminary window delay.
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

function [L1A, L1AP]...
      = window_delay (L1A, L1AP, cnf, chd, cst)

%% 0.Initialisation and Memory allocation
    
    
%     %% Computation of H0
%     
%     for i_pulse=1:chd.N_pulses_burst
%         cor2n(i_pulse) = ((L1A.burst-1)*chd.N_pulses_burst + (i_pulse-1)) * L1A.cor2_comp_sar_isp;
%         if cor2n(i_pulse) > 0
%             cor2_inc(i_pulse) = floor(cor2n(i_pulse)/N_pulses_rc);
%         else
%             cor2_inc(i_pulse) = ceil(cor2n(i_pulse)/N_pulses_rc);
%         end
%         h_sar(i_pulse) = L1A.h0_comp_sar_isp* chd.h0_cor2_unit_conv + cor2_inc(i_pulse);
% 
%     %% Separation of CAI and FAI
%         cai_namb_sar(i_pulse) = floor(h_sar(i_pulse) / chd.cai_cor2_unit_conv);
%         fai_sar(i_pulse) = h_sar(i_pulse) - cai_namb_sar(i_pulse)* chd.cai_cor2_unit_conv;
%         
%     %% Computation of Rx delay
%         %rx_delay(i_pulse) = h_sar(i_pulse)*T0_sar_pre_dat/chd.T0_h0_unit_conv/chd.h0_cor2_unit_conv;
%         
%         rx_delay_nom(i_pulse) = h_sar(i_pulse)*chd.T0_nom/chd.T0_h0_unit_conv/chd.h0_cor2_unit_conv;
%     end
    %% Attitude corrections
        z_cog_ant_corr_pitch    = - chd.x_cog_ant * sin(L1A.pitch_sar) - chd.z_cog_ant * cos(L1A.pitch_sar);
        delta_pitch             = 2 * z_cog_ant_corr_pitch / cst.c;
        z_cog_ant_corr_roll     = chd.y_cog_ant * sin(L1A.roll_sar);
        delta_roll              = 2 * z_cog_ant_corr_roll / cst.c;
        L1A.win_delay = L1A.win_delay + delta_pitch + delta_roll;
        
            
    %% Averaging
%    rx_delay_sar_ku_mean = mean(rx_delay_nom(1+chd.N_cal_pulses_burst+chd.N_c_pulses_burst:chd.N_pulses_burst));
    
%     mean_hn_sar_ku      = round(mean(h_sar(1+chd.N_cal_pulses_burst+chd.N_c_pulses_burst:chd.N_pulses_burst)));    
    
    %% Doppler correction
%     doppler_corr_sar    = -2 * L1A.alt_rate_sar_sat/ (chd.wv_length * chd.chirp_slope_ku);

%     %% Instrument corrections
%  	int_delay_att_sar= attenuator_table_delay_ku(1 + attcode_sar_ku_isp);
%  	int_delay_att_cal= attenuator_table_delay_ku(1 + attcode_cal1);
%     
%     int_delay_att=  - int_delay_att_cal;
% 	instr_delay= int_delay_cor_cal1+ int_delay_att+ chd.ext_delay_ground_ku;
%     instr_delay= chd.ext_delay_ground_ku;
%     %% Final window delay
%     win_delay_sar_ku_av= rx_delay_sar_ku_mean+ instr_delay;
%     win_delay_sar_ku_ref= rx_delay_sar(1+chd.N_cal_pulses_burst+chd.N_c_pulses_burst) + instr_delay;

    %% not used as there are no CAL or C pulses
    %     win_delay_sar_ku_ref_nom= rx_delay_sar_nom(1+chd.N_cal_pulses_burst+chd.N_c_pulses_burst) + instr_delay;
    



end
