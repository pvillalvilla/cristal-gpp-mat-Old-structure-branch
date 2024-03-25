% ONBOARD REVERSION ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% ---------------------------------------------------------
% Objective: - reverse the onboard RMC for the SAR RMC data.
% 
%  Author:   Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
%            Eduard Makhoul / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (11/09/2013)
%
% v1.0 2019/06/10 First version simplified imported from the s6 GPP
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L1A, chd] = onboard_reversion(L1A,filesBulk, chd,cnf,cst)

rmc_matrix_file= dir([filesBulk.auxPath '*RMC_matrix.mat']);
load([filesBulk.auxPath rmc_matrix_file.name]);

chd.N_samples_sar   = size(L1A.wfm_cal_gain_corrected,2)*2;
truncated_zeros1    = zeros(chd.N_ku_pulses_burst,chd.RMC_start_sample-1);
truncated_zeros2    = zeros(chd.N_ku_pulses_burst,chd.N_samples_sar/2-chd.RMC_start_sample+1);

win_delay_sar_ku    = 0;
win_delay_sar_ku_nom = 0;
uso_drift = 0;

wfm_5               = zeros(chd.N_ku_pulses_burst,chd.N_samples_sar);
% wfm_6               = zeros(chd.N_ku_pulses_burst,chd.N_samples_sar);
% wfm_7              = zeros(chd.N_ku_pulses_burst,chd.N_samples_sar);
% wfm_8               = zeros(chd.N_ku_pulses_burst,chd.N_samples_sar);
wfm_9               = zeros(chd.N_ku_pulses_burst,chd.N_samples_sar);
wfm_10               = zeros(chd.N_ku_pulses_burst,chd.N_samples_sar);
wfm_11               = zeros(chd.N_ku_pulses_burst,chd.N_samples_sar);
wfm_sar_reversed    = zeros(chd.N_ku_pulses_burst,chd.N_samples_sar);
alt_rate_wd_corr    = zeros(chd.N_ku_pulses_burst);
alt_rate_wd_corr_nom = zeros(chd.N_ku_pulses_burst);

%% 1. TRUNCATION reversion

    if(chd.RMC_start_sample==1)
        wfm_1   = (cat(2,L1A.wfm_cal_gain_corrected(:,1:chd.N_samples_sar/2),truncated_zeros2));
    else
        wfm_1   = (cat(2,truncated_zeros1,L1A.wfm_cal_gain_corrected(:,1:chd.N_samples_sar/2),truncated_zeros2));
    end


    
        %% 2. FFT range dimension: frequency to time
        %wfm_2(:,:) = fftshift(fft(wfm_1_.').',2); 
        wfm_2(:,:) = fftshift(fft(wfm_1,chd.N_samples_sar,2),2);
%         for i_beam=1:size(rmc_matrix,2)
%             rmc_matrix_(:,i_beam) =  interp1(1: chd.N_samples_sar,   rmc_matrix(:,i_beam) , 1:2: chd.N_samples_sar);
%         end
        
        %% 3. RMC reversion
        wfm_3 = wfm_2.* (conj(rmc_matrix(1:2:end,:)).');
        %% 4. IFFT azimuth dimension: doppler beams to pulses
        wfm_4 = sqrt(chd.N_ku_pulses_burst)*(ifft(fftshift(wfm_3,1)));
    
%         for i_pulse = 1:chd.N_ku_pulses_burst
%             %% 5. Doppler centroid reversion
%             wfm_5(i_pulse,:) = wfm_4(i_pulse,:) * exp(...
%                 -1i*2*cst.pi*(i_pulse-1)*delta_tau_isp(i_burst)*...
%                 T0_sar_pre_dat(i_burst)/(2^22)*chd.freq);
%             
            %% 6. Radial speed correction reversion
%             for i_sample = 1:chd.N_samples_sar
%                 wfm_6(i_pulse,i_sample) = wfm_5(i_pulse,i_sample)...
%                     * exp(-1i*2*cst.pi*(i_pulse-1)*delta_tau_isp(i_burst)/(2^22)...
%                     * (i_sample-1-chd.N_samples_sar/2) / chd.N_samples_sar);
%             end
            
%             %% 7. CAI/FAI reversion
%             for i_sample = 1:chd.N_samples_sar
%                 wfm_7(i_pulse,i_sample) = wfm_6(i_pulse,i_sample) ...
%                     * exp(+1i*2*cst.pi...
%                     * ((cai_namb_sar(i_pulse)...
%                     -cai_namb_sar(1))*cai_cor2_unit_conv_chd...
%                     + fai_sar(i_pulse)...
%                     -fai_sar(1))...
%                     / chd.h0_cor2_unit_conv / chd.T0_h0_unit_conv * ...
%                     (i_sample-1-chd.N_samples_sar/2)/chd.N_samples_sar);
%             end
% 
%             %% 8. Undo azimuth weighting
%             wfm_8(i_pulse,:)=wfm_7(i_pulse,:)./rmc_weights(i_pulse);
%             
%             
%         end 
            %% 9. IFFT range dimension: time to frequency
            L1A.wfm_cal_gain_corrected=[];
            L1A.wfm_cal_gain_corrected(:,:)= fftshift(ifft(ifft(fftshift(wfm_4.',1),chd.N_samples_sar)),1).'; % aligned with P4 model 4.3
            %wfm_9(:,:)= (ifft(fftshift(ifft(fftshift(wfm_8.',1)),1))).'; % aligned with P4 model 4.0
            %wfm_9(:,:)= fftshift(fft(ifft(fftshift(wfm_8.',1))),1).';

%             figure; mesh(abs(fftshift(fft(L1A.wfm_cal_gain_corrected.'),2).')); % aligned with P4 model 4.3
% 
    
    
    %% 10. ALTITUDE RATE ALIGNMENT
%     for i_pulse = 1:chd.N_ku_pulses_burst
%         
%         %% A. CAI/FAI RE-REVERSION
%         for i_sample = 1:chd.N_samples_sar
%             wfm_10(i_pulse,i_sample) = wfm_9(i_pulse,i_sample) ...
%                 * exp(-1i*2*cst.pi...
%                 * ((cai_namb_sar(i_pulse)...
%                 -cai_namb_sar(1))*cai_cor2_unit_conv_chd + ...
%                 fai_sar(i_pulse)...
%                 -fai_sar(1))...
%                 / chd.h0_cor2_unit_conv / chd.T0_h0_unit_conv * ...
%                 (i_sample-1-chd.N_samples_sar/2)/chd.N_samples_sar);
%         end   
%         
%         %% A_prima CAL-2 application
%         if cnf.flag_cal2_correction
%             aux = fftshift(fft(wfm_10(i_pulse,:),chd.N_samples_sar,2),2); % fft along range samples
%             wfm_11(i_pulse,:)=ifft(fftshift(aux.*wfm_cal2_science_sar,2),chd.N_samples_sar,2);
%         else
%             wfm_11(i_pulse,:)=wfm_10(i_pulse,:);
%         end
        
        %% B. ALTITUDE RATE CORRECTION
%         for i_sample = 1:chd.N_samples_sar
%                 wfm_sar_reversed(i_pulse,i_sample) = wfm_11(i_pulse,i_sample) ...
%                     * exp(1i*2*cst.pi*alt_rate_sar_sat(i_burst)*(i_pulse-1)...
%                     *pri_sar_pre_dat(i_burst) * 2/c_cst / T0_sar_pre_dat(i_burst)...
%                     * (i_sample-1-chd.N_samples_sar/2)/chd.N_samples_sar);
%         end
%         %% 12. Window delay update
%         alt_rate_wd_corr(i_pulse) = alt_rate_sar_sat(i_burst) *...
%             (i_pulse-1) * pri_sar_pre_dat(i_burst) * 2/c_cst; % in seconds
% 
%         alt_rate_wd_corr_nom(i_pulse) = alt_rate_sar_sat(i_burst) *...
%             (i_pulse-1) * chd.pri_nom(i_burst) * 2/c_cst; % in seconds
%     end
% 
%     win_delay_sar_ku(i_burst) = win_delay_sar_ku_ref(i_burst) + mean(alt_rate_wd_corr(:));
%     win_delay_sar_ku_nom(i_burst) = win_delay_sar_ku_ref_nom(i_burst) + mean(alt_rate_wd_corr_nom(:));
%     uso_drift(i_burst) = win_delay_sar_ku(i_burst) - win_delay_sar_ku_nom(i_burst);
%    

end
