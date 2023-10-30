%% HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
% ---------------------------------------------------------
% Objective: The purpose of the range compression is to perform range 
% transformation of the input stacks (with an FFT) and then generate the 
% power waveforms
% ----------------------------------------------------------
% Author:    Albert Garcia  / isardSAT
%            Eduard Makhoul / isardSAT
% ----------------------------------------------------------
% Version  record
% 1.0 2018/07/18 First version imported from Dedop rev 125
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [L1BS] = range_transformation (L1BS, cnf, chd, cst)



% L1BS.beams_rng_cmpr = zeros (max(L1BS.N_beams_stack), chd.N_samples_sar * cnf.zp_fact_range*L1BS.N_windows);
% if(strcmp(cnf.processing_mode,'SIN'))
%     L1BS.beams_rng_cmpr_2 = zeros (max(L1BS.N_beams_stack), chd.N_samples_sar * cnf.zp_fact_range*L1BS.N_windows);
%     L1BS.phase_diff     = zeros (max(L1BS.N_beams_stack), chd.N_samples_sar *cnf.zp_fact_range*L1BS.N_windows);
%     L1BS.coherence      = zeros (max(L1BS.N_beams_stack), chd.N_samples_sar * cnf.zp_fact_range*L1BS.N_windows);
%     L1BS.coherence_mask= zeros (max(L1BS.N_beams_stack), chd.N_samples_sar * cnf.zp_fact_range*L1BS.N_windows);
% end
% beams_rng_cmpr_I = zeros (L1BS.N_total_surf_loc, max(L1BS.N_beams_stack), chd.N_samples_sar * cnf.zp_fact_range);
% beams_rng_cmpr_Q = zeros (L1BS.N_total_surf_loc, max(L1BS.N_beams_stack), chd.N_samples_sar * cnf.zp_fact_range);
% beams_rng_cmpr_aux = zeros (L1BS.N_total_surf_loc, max(L1BS.N_beams_stack), chd.N_samples_sar * cnf.zp_fact_range);
   
    

    %% 1.FFT
	
    
			if(cnf.window_rg)
				% Hamming        
				beams_zp_fft    = (fft((fftshift(squeeze(L1BS.beam_geo_corr(:,:)),2).*(ones(length(L1BS.beam_geo_corr(:,1)),1)*window_range)).', chd.N_samples_sar * cnf.zp_fact_range)).'/sqrt(chd.N_samples_sar* cnf.zp_fact_range);
				if(strcmp(cnf.processing_mode,'SIN')) && strcmp(cnf.processing_mode,'SIN') 
					beams_zp_fft_2  = (fft((fftshift(squeeze(L1BS.beam_geo_corr_2(:,:)),2).*(ones(length(L1BS.beam_geo_corr(:,1)),1)*window_range)).', chd.N_samples_sar * cnf.zp_fact_range)).'/sqrt(chd.N_samples_sar* cnf.zp_fact_range);
				end
            else
                
                    non_centered_spectra=fftshift(squeeze(L1BS.beam_geo_corr(:,:)),2);
                    beams_zp_fft = fft([non_centered_spectra(:,1:chd.N_samples_sar/2),...
                    zeros(size(non_centered_spectra,1),(cnf.zp_fact_range-1)*chd.N_samples_sar),...
                    non_centered_spectra(:,chd.N_samples_sar/2+1:chd.N_samples_sar)],...
                    chd.N_samples_sar * cnf.zp_fact_range,2);
% 				beams_zp_fft    = (fft(fftshift(squeeze(L1BS.beam_geo_corr(:,:).'),1), chd.N_samples_sar * cnf.zp_fact_range)).'/sqrt(chd.N_samples_sar* cnf.zp_fact_range);
				if strcmp(cnf.processing_mode,'SIN')
                    non_centered_spectra_2=fftshift(squeeze(L1BS.beam_geo_corr_2(:,:)),2);
                    beams_zp_fft_2 = fft([non_centered_spectra_2(:,1:chd.N_samples_sar/2),...
                    zeros(size(non_centered_spectra_2,1),(cnf.zp_fact_range-1)*chd.N_samples_sar),...
                    non_centered_spectra_2(:,chd.N_samples_sar/2+1:chd.N_samples_sar)],...
                    chd.N_samples_sar * cnf.zp_fact_range,2);
% 					beams_zp_fft_2  = (fft(fftshift(squeeze(L1BS.beam_geo_corr_2(:,:).'),1), chd.N_samples_sar * cnf.zp_fact_range)).'/sqrt(chd.N_samples_sar* cnf.zp_fact_range);
				end
			end

    
    	%% 2. AMPLITUDE to POWER
        
        L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,:)   = abs(beams_zp_fft(1:L1BS.N_beams_stack,:)).^2;
        L1BS.beams_rng_cmprIQ(1:L1BS.N_beams_stack,:) = (beams_zp_fft(1:L1BS.N_beams_stack,:));
        
        if(strcmp(cnf.processing_mode,'SIN')) && strcmp(cnf.processing_mode,'SIN')
            
            L1BS.beams_rng_cmpr_2(1:L1BS.N_beams_stack,:) = abs(beams_zp_fft_2(1:L1BS.N_beams_stack,:)).^2;
            L1BS.beams_rng_cmprIQ_2(1:L1BS.N_beams_stack,:)= (beams_zp_fft_2(1:L1BS.N_beams_stack,:));
            if(cnf.stack_interferometry_flag == 1)
                % substacking performed, otherwise coherence will be always 1
                for i_beam=1:L1BS.N_beams_stack
                    beam_margin_l=2;
                    beam_margin_r=2;
                    N_substack = beam_margin_l+ beam_margin_r+1;
                    if(i_beam-beam_margin_l<1)
                        beam_margin_l=i_beam-1;
                        beam_margin_r = N_substack-beam_margin_l-1;
                        continue; % do not compute anything in the edges
                        
                    elseif((i_beam+beam_margin_r)>L1BS.N_beams_stack)
                        beam_margin_r=L1BS.N_beams_stack-i_beam;
                        beam_margin_l=N_substack-beam_margin_r-1;
                        continue;
                    end
                    beams = (i_beam-beam_margin_l):(i_beam+beam_margin_r);
                    aux_1 = L1BS.beams_rng_cmprIQ(beams,:);
                    aux_2 = L1BS.beams_rng_cmprIQ_2(beams,:);
                    aux_1(find(aux_1==0)) = nan;
                    aux_2(find(aux_2==0)) = nan;
                    stack_power(i_beam,:) = nanmean(L1BS.beams_rng_cmprIQ(beams,:).*conj(L1BS.beams_rng_cmprIQ_2(beams,:)));
                    a1_power(i_beam,:) = (nanmean(L1BS.beams_rng_cmpr(beams,:)));
                    a2_power(i_beam,:) = (nanmean(L1BS.beams_rng_cmpr_2(beams,:)));
                    if(L1BS.Gap_flag(i_beam))
                        L1BS.coherence(i_beam,:) = (abs(stack_power(i_beam,:)))./sqrt(a1_power(i_beam,:).*a2_power(i_beam,:));
                        L1BS.coherence_mask(i_beam,find(L1BS.coherence(i_beam,:)>cnf.coherence_threshold)) = 1; % to be checked

                    end
%                     stack_power = beams_zp_fft(1:L1BS.N_beams_stack,:).*conj(beams_zp_fft_2(1:L1BS.N_beams_stack,:));
%                     a1_power = squeeze(L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,:));
%                     a2_power = squeeze(L1BS.beams_rng_cmpr_2(1:L1BS.N_beams_stack,:));
%                     L1BS.coherence(1:L1BS.N_beams_stack,:) = sqrt((abs(stack_power).^2)./(a1_power.*a2_power));
%                     L1BS.coherence_mask(1:L1BS.N_beams_stack,find(L1BS.coherence(1:L1BS.N_beams_stack,:)>cnf.coherence_threshold)) = 1; % to be checked
                end
                 phase1 = atan(imag(beams_zp_fft(1:L1BS.N_beams_stack,:))./real(beams_zp_fft(1:L1BS.N_beams_stack,:)))*180/pi;
                 phase2 = atan(imag(beams_zp_fft_2(1:L1BS.N_beams_stack,:))./real(beams_zp_fft_2(1:L1BS.N_beams_stack,:)))*180/pi;
                 L1BS.phase_diff(1:L1BS.N_beams_stack,:) = phase1 - phase2;
                    
            %Based on the conventional CryoSAt-2 SARIn processing:
            %coherencey and phase info is after multilooking
            
%             
            end
                
%         %% 3. FORCING WRAPPED SAMPLES TO 0 (applying the mask)
%         if abs(L1BS.shift_coarse(i_beam)) < chd.N_samples_sar
%             if L1BS.shift_coarse(i_beam) > 0
%                 limit1 = 1;
%                 limit2 = cnf.zp_fact_range * L1BS.shift_coarse(i_beam);
%             else
%                 limit1 = (chd.N_samples_sar - abs(L1BS.shift_coarse(i_beam)))*cnf.zp_fact_range;
%                 limit2 = chd.N_samples_sar*cnf.zp_fact_range;
%             end
%         else
%             limit1 = 1;
%             limit2 = chd.N_samples_sar*cnf.zp_fact_range;
%         end
%         L1BS.beams_rng_cmpr(i_beam,limit1:limit2) = 0;
    
   

end

