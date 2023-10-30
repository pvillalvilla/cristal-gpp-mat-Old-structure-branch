% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code runs an TCOG retracker based on Technical note CryoVal-LI
% -------------------------------------------------------------------------
% 
% Author:           Albert Mondejar / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT isardSAT V1 18/04/2019L1B.scaled_waveforms.'

function [fit_res] = leading_edge_detector (L1B,cnf_L2)

[mean_noise,~]     = compute_noise_floor(L1B.scaled_waveforms.',cnf_L2);










end
