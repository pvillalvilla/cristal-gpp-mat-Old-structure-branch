%% HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
% ---------------------------------------------------------
% Objective: The purpose of the surface & window delay alignment is to coregister the 
% different leading edge positions of the surfaces over the track. Need to
% align also the window delay.
% ----------------------------------------------------------
% Author:    Albert Garcia  / isardSAT
%            Eduard Makhoul / isardSAT
% ----------------------------------------------------------
% Version  record
% 1.0 2018/07/18 First version imported from Dedop rev 125
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [L1BS,L1B] = surface_win_delay_alignment (L1BS,L1B,win_delay_surf_ref,alt_sat_ref, cnf, chd, cst)


wd_shift = ((L1BS.win_delay_surf-win_delay_surf_ref)-(L1BS.alt_sat-alt_sat_ref)*2/cst.c) / chd.T0;
L1B.wfm_cor_i2q2_wdcorr = circshift(L1B.wfm_cor_i2q2,[0,round(wd_shift*cnf.zp_fact_range)]);
L1BS.win_delay_surf_aligned = L1BS.win_delay_surf-wd_shift*chd.T0;
if strcmp(cnf.processing_mode,'SIN') && strcmp(cnf.processing_mode,'SIN')
   L1B.wfm_cor_i2q2_2_wdcorr = circshift(L1B.wfm_cor_i2q2_2,[0,round(wd_shift*cnf.zp_fact_range)]);
   L1B.phase_difference_wdcorr = circshift(L1B.phase_difference,[0,round(wd_shift*cnf.zp_fact_range)]);
   L1B.coherence_wdcorr = circshift(L1B.coherence,[0,round(wd_shift*cnf.zp_fact_range)]);
end
end


