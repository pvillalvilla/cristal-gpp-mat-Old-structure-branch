%% HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
% ---------------------------------------------------------
% Objective: The purpose of the geometry corrections is to compute and 
% apply Doppler correction,Slant range correction and Window delay 
% Misalignments
% ----------------------------------------------------------
% Author:    Albert Garcia  / isardSAT
%            Eduard Makhoul / isardSAT
% ----------------------------------------------------------
% Version  record
% 1.0 2018/07/18 First version imported from Dedop rev 125
% 1.1 2019/06/10 Removed geometry mask computation as it is now done in the masking function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [L1BS]  = geometry_corrections (L1BS, cnf, chd, cst)
                                                                                             



% Pre-allocating memory
doppler_range           = zeros(1,chd.N_max_beams_stack);
i_samples               = (0:(chd.N_samples_sar-1));
win_delay_ref           = 0;
beam_power              = zeros(1,chd.N_max_beams_stack);
beam_power2             = zeros(1,chd.N_max_beams_stack);


% L1BS.doppler_corr       = zeros(1,max(L1BS.N_beams_stack));
% L1BS.range_sat_surf          = zeros(1,max(L1BS.N_beams_stack));
% L1BS.slant_range_corr_time   = zeros(1,max(L1BS.N_beams_stack));
% L1BS.slant_range_corr        = zeros(1,max(L1BS.N_beams_stack));
% L1BS.wd_corr                 = zeros(1,max(L1BS.N_beams_stack));
% L1BS.shift                   = zeros(1,max(L1BS.N_beams_stack));
% L1BS.shift_coarse            = zeros(1,max(L1BS.N_beams_stack));
%L1BS.good_samples            = zeros(max(L1BS.N_beams_stack),chd.N_samples_sar*cnf.zp_fact_range);
%              
    
   %0 = normal alignment
    %1 = optimal alignment (wrt max beam accumulated power)
    %2 = force the alignment to the closest of a given window_delay
    %3 = force to the first one.
    %4 = force to the minimum one.
   
    %% 1. BEAM POWER
    if cnf.win_delay_ref == 0
        win_delay_ref = L1BS.win_delay_surf;
		[~,L1BS.beam_ref] = min(abs(win_delay_ref-L1BS.win_delay_beam(1:L1BS.N_beams_stack)));
    elseif cnf.win_delay_ref ==1
        beams_surf_aux1_fft = (fft(fftshift(L1BS.beams_surf(:,:),2).')).';
        
        beam_power(1:L1BS.N_beams_stack) = sum(abs(beams_surf_aux1_fft(1:L1BS.N_beams_stack,:)).^2,2);
        
        [~,L1BS.beam_ref] = max(beam_power(1:L1BS.N_beams_stack));
        win_delay_ref = L1BS.win_delay_beam(L1BS.beam_ref);
        L1BS.win_delay_surf = win_delay_ref;
    elseif cnf.win_delay_ref==2
       %need testing, not working
%         win_delay_ref_index_beams(:) = (L1A.alt_sar_sat((1:L1BS.N_beams_stack))- cnf.elevation_ref) * 2 / cst.c;
        [~,L1BS.beam_ref] = min(abs(win_delay_ref-L1BS.win_delay_beam(1:L1BS.N_beams_stack)));
    elseif cnf.win_delay_ref==3 
		L1BS.beam_ref=1;
        win_delay_ref = L1BS.win_delay_beam(L1BS.beam_ref);
		L1BS.win_delay_surf = win_delay_ref;
    elseif cnf.win_delay_ref==4
        [win_delay_ref,L1BS.beam_ref] = min(L1BS.win_delay_beam(1:L1BS.N_beams_stack));
		L1BS.win_delay_surf = win_delay_ref;
	elseif cnf.win_delay_ref==5
        [win_delay_ref,L1BS.beam_ref] = max(L1BS.win_delay_beam(1:L1BS.N_beams_stack));
		L1BS.win_delay_surf = win_delay_ref;
    end
      
    
        %% 2. DOPPLER CORRECTION
       norm_vel_sat(1:L1BS.N_beams_stack) = sqrt(L1BS.x_vel_sat_beam(1:L1BS.N_beams_stack).^2+ L1BS.y_vel_sat_beam(1:L1BS.N_beams_stack).^2+ L1BS.z_vel_sat_beam(1:L1BS.N_beams_stack).^2);

        
        % Range correction computation
        doppler_range(1:L1BS.N_beams_stack) = (- cst.c / chd.wv_length .* norm_vel_sat(1:L1BS.N_beams_stack) .* cos(L1BS.beam_ang_surf(1:L1BS.N_beams_stack))) .* chd.pulse_length / chd.bw;
        % Doppler correction computation
        L1BS.doppler_corr(1:L1BS.N_beams_stack) = 2 / cst.c .* doppler_range(1:L1BS.N_beams_stack) ./ L1BS.T0_sar_surf(1:L1BS.N_beams_stack);
%         L1BS.doppler_corr(i_beam) = 0; %testing
        L1BS.doppler_corr(isnan(L1BS.doppler_corr))=0;
        %% 3. SLANT RANGE CORRECTION
        
%         % Selection of the observed surface locations
%         clear x_seen y_seen z_seen
%         % if i_burst == 822
%         %     surf_loc_index(1) = surf_loc_index(1) + 1;
%         % end
%         L1BS.surf_counter
%         % 20 is a first step, the number should be computed as
%         % ceil((1.35*pi/180*vel_sat_norm*2/wv_length_ku/143.4)-64)
%         start_beam_amb(i_burst) = max(1,surf_loc_index(i_burst,1)-20)-1;
%         stop_beam_amb(i_burst)  = min(max(max(surf_loc_index(:,:))),surf_loc_index(i_burst,N_beams_sar_ku(i_burst))+20)-1;
%         x_seen (1:N_beams_sar_ku(i_burst)) = x_surf (surf_loc_index(i_burst,1)-start_beam_amb : surf_loc_index(i_burst,N_beams_sar_ku(i_burst))-stop_beam_amb);
%         y_seen (1:N_beams_sar_ku(i_burst)) = y_surf (surf_loc_index(i_burst,1)-start_beam_amb : surf_loc_index(i_burst,N_beams_sar_ku(i_burst))-stop_beam_amb);
%         z_seen (1:N_beams_sar_ku(i_burst)) = z_surf (surf_loc_index(i_burst,1)-start_beam_amb : surf_loc_index(i_burst,N_beams_sar_ku(i_burst))-stop_beam_amb);
        
        
        % Range computation
        range_sat_surf(:,1:L1BS.N_beams_stack) = [L1BS.x_sar_sat_beam((1:L1BS.N_beams_stack)); ...
                                               L1BS.y_sar_sat_beam((1:L1BS.N_beams_stack)); ...
                                               L1BS.z_sar_sat_beam((1:L1BS.N_beams_stack))] - [L1BS.x_surf; L1BS.y_surf; L1BS.z_surf]*ones(1,L1BS.N_beams_stack);
        L1BS.range_sat_surf(1:L1BS.N_beams_stack) = sqrt((range_sat_surf(1,:)).^2+ (range_sat_surf(2,:)).^2+ (range_sat_surf(3,:)).^2);
        % Range delay computation
        L1BS.slant_range_corr_time(1:L1BS.N_beams_stack) = L1BS.win_delay_surf - (L1BS.range_sat_surf(1:L1BS.N_beams_stack) .* 2 / cst.c);
        L1BS.slant_range_corr(1:L1BS.N_beams_stack) = L1BS.slant_range_corr_time(1:L1BS.N_beams_stack) ./ L1BS.T0_sar_surf(1:L1BS.N_beams_stack);
    
        %% 4. WINDOW DELAY MISALIGNMENTS CORRECTION
        L1BS.wd_corr(1:L1BS.N_beams_stack) = - (win_delay_ref - L1BS.win_delay_beam((1:L1BS.N_beams_stack))) ./ L1BS.T0_sar_surf(1:L1BS.N_beams_stack);
    
    
        L1BS.shift = -L1BS.doppler_corr.*cnf.apply_doppler + L1BS.slant_range_corr.*cnf.apply_slant + L1BS.wd_corr.*cnf.apply_wd; 
        % the Doppler correction Sign is changed after teleconf with Michele as they are generating the sinusoidal with a negative carrier frequency
        


L1BS.shift_coarse = round(L1BS.shift);
L1BS.shift_fine = L1BS.shift-L1BS.shift_coarse;


%% increase the window size 
if(cnf.avoid_wrapps)
    L1BS.N_windows  = ceil(max(max(abs(L1BS.shift))/chd.N_samples_sar)+1);
    truncated_zeros = zeros(1,chd.N_samples_sar*cnf.zp_fact_range*(L1BS.N_windows-1));
else
    L1BS.N_windows  = 1;
    truncated_zeros =[];
end

% L1BS.beam_geo_corr          = zeros(1,max(L1BS.N_beams_stack),chd.N_samples_sar*L1BS.N_windows);
% 
%  if(strcmp(cnf.processing_mode,'SIN'))
%     L1BS.beam_geo_corr_2 	= zeros(1,max(L1BS.N_beams_stack),chd.N_samples_sar*L1BS.N_windows);
%  end



%% 5. APPLYING THE ALIGNMENTS IN FREQUENCY (PICE) 
    
    
    
    
      
        L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:) = L1BS.beams_surf(1:L1BS.N_beams_stack,:) .* exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples);        
        
%         figure;imagesc(abs(fft(L1BS.beams_surf(1:L1BS.N_beams_stack,:).').'))
%         hold all
%         plot(-L1BS.shift,1:length(L1BS.doppler_corr))
%         plot(-L1BS.slant_range_corr,1:length(L1BS.doppler_corr))
%         plot(-L1BS.wd_corr,1:length(L1BS.doppler_corr))
%         figure;imagesc(abs(fft(L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:).').'))
        
        if(strcmp(cnf.processing_mode,'SIN')) && strcmp(cnf.processing_mode,'SIN')
%          	wfm_geo_corr_aux_2(i_beam,:)    = squeeze(L1BS.beams_surf_2(i_beam,:)).' .* exp(2i*cst.pi/chd.N_samples_sar*L1BS.shift(i_beam).*i_samples);
%             beam_zp_fft                     = fftshift(fft(wfm_geo_corr_aux_2(i_beam,:).', chd.N_samples_sar * cnf.zp_fact_range),1).'/sqrt(chd.N_samples_sar);
%             beam_surf_aux1_fft_increased    = (cat(2,beam_zp_fft,truncated_zeros));
%             beam_surf_aux1_fft_aligned      = circshift(beam_surf_aux1_fft_increased,[0,L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
%             L1BS.beam_geo_corr_2(i_beam,:) = beam_surf_aux1_fft_aligned;
             L1BS.beam_geo_corr_2(1:L1BS.N_beams_stack,:)    = L1BS.beams_surf_2(1:L1BS.N_beams_stack,:) .* exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples);
        end     


    for i_beam = 1:L1BS.N_beams_stack
        if(L1BS.shift_coarse(i_beam)>0 && L1BS.shift_coarse(i_beam)<chd.N_samples_sar)
            start_sample=ceil(L1BS.shift_coarse(i_beam))*cnf.zp_fact_range+1;
            final_sample=chd.N_samples_sar*cnf.zp_fact_range;
            L1BS.good_samples(i_beam,1:start_sample-1)=0;
            L1BS.good_samples(i_beam,start_sample:final_sample)=1;
        elseif(L1BS.shift_coarse(i_beam)<=0 && L1BS.shift_coarse(i_beam)>-1*chd.N_samples_sar)
            start_sample=1;
            final_sample=(chd.N_samples_sar*cnf.zp_fact_range+floor((L1BS.shift_coarse(i_beam))*cnf.zp_fact_range));
            L1BS.good_samples(i_beam,start_sample:final_sample)=1;
            L1BS.good_samples(i_beam,final_sample+1:chd.N_samples_sar*cnf.zp_fact_range)=0;
        else
            L1BS.good_samples(i_beam,1:chd.N_samples_sar*cnf.zp_fact_range)=0;
        end
    end
    

% L1BS = rmfield(L1BS,'beams_surf');
% if(strcmp(cnf.processing_mode,'SIN'))
%     L1BS = rmfield(L1BS,'beams_surf_2');
% end
end

    


