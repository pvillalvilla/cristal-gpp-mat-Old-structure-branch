%% HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% ---------------------------------------------------------
% This code implements the ambiguity mask as described in the
% isardSAT_GPPICE_ATBD_v1a
%
% ---------------------------------------------------------
% Objective: The purpose of this function is to generate the ambiguity mask to be used 
% in the stack masking procedure.
% 
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Albert Garcia / isardSAT (18/05/2019)
% 
% Revisions:

% v1.0 2019/06/10 First Version
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ambiguity_mask] = ambiguity_mask_construction(L1BS, cnf, chd, cst)

ambiguity_mask_margin_hr_cnf = 5; %fixed from S6

 

%% ----------- Initialize variables ---------------------------------------
ambiguity_mask      = ones(L1BS.N_beams_stack,...
                           chd.N_samples_sar*cnf.zp_fact_range);
leftambiguity_mask  = ones(L1BS.N_beams_stack,...
                           chd.N_samples_sar*cnf.zp_fact_range);
rightambiguity_mask = ones(L1BS.N_beams_stack,...
                           chd.N_samples_sar*cnf.zp_fact_range);

%% -------- Compute rough estimation of leading edge position -------------
N_beams_selected = 10;

[~,doppler_central] = min(abs(cst.pi/2-L1BS.beam_ang_surf));
if(doppler_central> L1BS.N_beams_stack - N_beams_selected)
    right_beams = L1BS.N_beams_stack-doppler_central;
else
    right_beams = N_beams_selected;
end
if(doppler_central < N_beams_selected + 1)
    left_beams = doppler_central-1;
else
    left_beams = N_beams_selected;
end
[~,max_pos_beams] = (max(L1BS.beams_rng_cmpr,[],2));
leading_edge = round(min(max_pos_beams(doppler_central-left_beams:doppler_central+right_beams)));

%% -------- Construct the mask --------------------------------------------
% Compute earth curvature factor
alpha_R = 1+L1BS.alt_sat/cst.earth_radius;

%Compute the satellite ground velocity
v_s = norm([L1BS.x_vel_sat,L1BS.y_vel_sat,L1BS.z_vel_sat]);
v_g = v_s/alpha_R;

%Compute the equivalent/effective satellite velocity
v_e = sqrt(v_s*v_g);

%Compute the azimuth chirp rate
Ka = (2.0*v_e.^2)./(chd.wv_length* L1BS.alt_sat);

for i_beam = 1:L1BS.N_beams_stack
    %Compute the range sampling [meters]
    deltaRange = cst.c/(2.0*cnf.zp_fact_range)*L1BS.T0_sar_surf(i_beam);
            
    %Compute the PRF
    PRF = 1./L1BS.pri_sar_sat_beam(i_beam);
    
    %Compute the Doppler frequency
    fd = 2/chd.wv_length*v_s*cos(L1BS.beam_ang_surf(i_beam));
    
    
    %Compute the (hyperbolic approx.) non-ambiguous slant range
    %variation in Doppler domain
    deltaR_noamb = (1.0* L1BS.alt_rate_sat)/Ka*fd + ...
    chd.wv_length^2.*L1BS.alt_sat.*fd.^2./(8.0*v_e.^2);
    
    %Compute the (hyperbolic approx.) 1st left-ambiguous slant range
    %variation in Doppler domain
    deltaR_lamb = (1.0* L1BS.alt_rate_sat)/Ka*(fd+PRF) + ...
    chd.wv_length^2.*L1BS.alt_sat.*(fd+PRF).^2./(8.0*v_e.^2);
    
    %Compute the (hyperbolic approx.) 1st right-ambiguous slant range
    %variation in Doppler domain
    deltaR_ramb = (1.0* L1BS.alt_rate_sat)/Ka*(fd-PRF) + ...
    chd.wv_length^2.*L1BS.alt_sat.*(fd-PRF).^2./(8.0*v_e.^2);
    
    %Residual range cell migration on left-ambiguity [in samples]
    res_RCMC_lamb = ceil((deltaR_lamb-deltaR_noamb)./deltaRange);
    
    %Residual range cell migration on right-ambiguity [in samples]
    res_RCMC_ramb = ceil((deltaR_ramb-deltaR_noamb)./deltaRange);
    
    %First ambiguous range bin (left-ambiguity)
    range_bin_lamb = leading_edge+res_RCMC_lamb - ...
    ambiguity_mask_margin_hr_cnf*cnf.zp_fact_range;
    
    %First ambiguous range bin (right-ambiguity)
    range_bin_ramb = leading_edge+res_RCMC_ramb - ...
    ambiguity_mask_margin_hr_cnf*cnf.zp_fact_range;
    
    
    %Left Ambiguity mask
    if range_bin_lamb>1 && range_bin_lamb<=chd.N_samples_sar*cnf.zp_fact_range
        start_sample = range_bin_lamb;
        final_sample = chd.N_samples_sar*cnf.zp_fact_range;
    elseif range_bin_lamb <=1
        start_sample = 1;
        final_sample = chd.N_samples_sar*cnf.zp_fact_range;
    end
    if range_bin_lamb<=chd.N_samples_sar*cnf.zp_fact_range
        leftambiguity_mask(i_beam,start_sample:final_sample) = 0;
    end
    
    %Right Ambiguity mask
    if range_bin_ramb>1 && range_bin_ramb<=chd.N_samples_sar*cnf.zp_fact_range
        start_sample = range_bin_ramb;
        final_sample = chd.N_samples_sar*cnf.zp_fact_range;
    elseif range_bin_ramb<=1
        start_sample = 1;
        final_sample = chd.N_samples_sar*cnf.zp_fact_range;
    end
    if range_bin_ramb<=chd.N_samples_sar*cnf.zp_fact_range
        rightambiguity_mask(i_beam,start_sample:final_sample) = 0;
    end
    
    ambiguity_mask(i_beam,:)=leftambiguity_mask(i_beam,:).*...
                             rightambiguity_mask(i_beam,:);
end



end

