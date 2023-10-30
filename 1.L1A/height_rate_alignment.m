%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% HEIGHT RATE APPLICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this function is to perform a better alignment with the 
% height rate information once the applied CAI/FAI have been compensated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%____________________________For each stack________________________________

function [burst] = height_rate_alignment (burst, cnf, chd, cst)

i_samples               = (0:(chd.N_samples_sar-1));
%i_pulses                = ((0:(chd.N_ku_pulses_burst-1))).';
i_pulses                = ((0:(chd.N_pulses_burst-1))-(chd.N_pulses_burst-1)/2).'; 

% %% ------------------------------ CAI/FAI compensation --------------------
% % In CryoSAt-2 no FAI alignment has been performed
% % Compensate only for CAI
% burst.wfm_cal_gain_corrected  = burst.wfm_cal_gain_corrected.* ...
%                     exp(-1i*2*cst.pi*...
%                     (((CAI-CAI(1)).*chd.cai_cor2_unit_conv / chd.h0_cor2_unit_conv / chd.T0_h0_unit_conv)*ones(1,chd.N_samples_sar)).*...
%                     (ones(chd.N_ku_pulses_burst,1)*i_samples/chd.N_samples_sar));

%% ----------------------------- Height rate application ------------------
% burst.wfm_cal_gain_corrected  = burst.wfm_cal_gain_corrected.* ...
%                     exp(-1i*2*cst.pi*(burst.alt_rate_sar_sat.*((i_pulses)*ones(1,chd.N_samples_sar)).*burst.pri_sar * 2/cst.c/burst.T0_sar+(FAI(1))./ chd.T0_h0_unit_conv).*...
%                     (ones(chd.N_pulses_burst,1)*i_samples/chd.N_samples_sar));
burst.wfm_cal_gain_corrected  = burst.wfm_cal_gain_corrected.* ...
                    exp(-1i*2*cst.pi*(burst.alt_rate_sar_sat.*((i_pulses)*ones(1,chd.N_samples_sar)).*burst.pri_sar * 2/cst.c/burst.T0_sar).*...
                    (ones(chd.N_pulses_burst,1)*i_samples/chd.N_samples_sar));                
             
                
if strcmp(cnf.processing_mode,'SIN')

    burst.wfm_cal_gain_corrected_2  = burst.wfm_cal_gain_corrected_2.* ...
                    exp(-1i*2*cst.pi*(burst.alt_rate_sar_sat.*((i_pulses)*ones(1,chd.N_samples_sar)).*burst.pri_sar * 2/cst.c/burst.T0_sar).*...
                    (ones(chd.N_pulses_burst,1)*i_samples/chd.N_samples_sar));               
    
end
%% ----------------------------- Window delay update ----------------------
%alt_rate_wd_corr = burst.alt_rate_sar_sat.*...
%            (i_pulses) * burst.pri_sar * 2/cst.c; % in seconds
%burst.win_delay_sar_ku=burst.win_delay_sar_ku+mean(alt_rate_wd_corr)+FAI(1)./ chd.T0_h0_unit_conv.*burst.T0_sar;
%burst.win_delay_sar_ku=burst.win_delay_sar_ku+mean(alt_rate_wd_corr);

end