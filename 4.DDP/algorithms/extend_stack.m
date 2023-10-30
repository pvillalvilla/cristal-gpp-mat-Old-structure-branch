%% HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
% ---------------------------------------------------------
% Objective: The purpose of the extend stack is to increase the window size
% taking into account the maximum range migration performed so the
% information is not lost
% ----------------------------------------------------------
% Author:    Albert Garcia  / isardSAT
%            Eduard Makhoul / isardSAT
% ----------------------------------------------------------
% Version  record
% 1.0 2018/07/18 First version imported from Dedop rev 125
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


function [L1BS] = extend_stack(L1BS, cnf, chd, cst)



% shift_extension = max(L1BS.shift)-min(L1BS.shift);
shift_foreward = max(L1BS.shift);
shift_backward = min(L1BS.shift);
sample_start(1:L1BS.N_beams_stack)= L1BS.shift_coarse(1:L1BS.N_beams_stack)*cnf.zp_fact_range-min(L1BS.shift_coarse(1:L1BS.N_beams_stack)*cnf.zp_fact_range)+1;
[~,L1BS.beam_ref2]=min(L1BS.shift_coarse(1:L1BS.N_beams_stack));
if(shift_backward>0)
    shift_backward=0;
else
    shift_backward=abs(shift_backward);
end
if(shift_foreward<0)
    shift_foreward=0;
end
% L1BS.N_windows  = ceil(shift_extension/chd.N_samples_sar)+1;
L1BS.N_windows_back  = ceil(shift_backward/chd.N_samples_sar);
L1BS.N_windows_fore  = ceil(shift_foreward/chd.N_samples_sar)+1;
L1BS.N_windows = L1BS.N_windows_fore+L1BS.N_windows_back;
truncated_zeros_forward = NaN(1,chd.N_samples_sar*cnf.zp_fact_range*(L1BS.N_windows_fore-1)); %additional windows to be added
truncated_zeros_back = NaN(1,chd.N_samples_sar*cnf.zp_fact_range*(L1BS.N_windows_back)); %additional windows to be added
if(isempty(truncated_zeros_forward))
    truncated_zeros_forward=[];
end
if(isempty(truncated_zeros_back))
    truncated_zeros_back=[];
end
beams_rng_cmpr_extended = NaN(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range*(L1BS.N_windows));
beams_rng_cmprIQ_extended = NaN(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range*(L1BS.N_windows));
if(strcmp(cnf.processing_mode,'SIN')) && strcmp(cnf.processing_mode,'SIN')
    
    beams_rng_cmpr_extended_2   = NaN(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range*(L1BS.N_windows));
    beams_rng_cmprIQ_extended_2 = NaN(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range*(L1BS.N_windows));
    phase_diff_extended_2       = NaN(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range*(L1BS.N_windows));
%     coherence_extended_2        = NaN(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range*(L1BS.N_windows));
%     coherence_mask_extended_2   = NaN(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range*(L1BS.N_windows));
end
    
for i_beam = 1:L1BS.N_beams_stack
%% Undo the aligning
beams_rng_cmpr_misaligned           = circshift(L1BS.beams_rng_cmpr(i_beam,:),[0,-L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
beams_rng_cmprIQ_misaligned         = circshift(L1BS.beams_rng_cmprIQ(i_beam,:),[0,-L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
%% Extend the window
beams_rng_cmpr_extended(i_beam,:)   = [truncated_zeros_back,beams_rng_cmpr_misaligned,truncated_zeros_forward];
beams_rng_cmprIQ_extended(i_beam,:) = [truncated_zeros_back,beams_rng_cmprIQ_misaligned,truncated_zeros_forward];
%% Reapply align the beams
beams_rng_cmpr_extended(i_beam,:)   = circshift(beams_rng_cmpr_extended(i_beam,:),[0,L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
beams_rng_cmprIQ_extended(i_beam,:) = circshift(beams_rng_cmprIQ_extended(i_beam,:),[0,L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);

if(strcmp(cnf.processing_mode,'SIN'))&& strcmp(cnf.processing_mode,'SIN')
    %% Undo the aligning
    beams_rng_cmpr_misaligned_2     = circshift(L1BS.beams_rng_cmpr_2(i_beam,:),[0,-L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
    beams_rng_cmprIQ_misaligned_2   = circshift(L1BS.beams_rng_cmprIQ_2(i_beam,:),[0,-L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
    phase_diff_misaligned_2         = circshift(L1BS.phase_diff(i_beam,:),[0,-L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
%     coherence_misaligned_2          = circshift(L1BS.coherence(i_beam,:),[0,-L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
%     coherence_mask_misaligned_2     = circshift(L1BS.coherence(i_beam,:),[0,-L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);

%% Extend the window
    beams_rng_cmpr_extended_2(i_beam,:)     = [truncated_zeros_back,beams_rng_cmpr_misaligned_2,truncated_zeros_forward];
    beams_rng_cmprIQ_extended_2(i_beam,:)	= [truncated_zeros_back,beams_rng_cmprIQ_misaligned_2,truncated_zeros_forward];
    phase_diff_extended_2(i_beam,:)         = [truncated_zeros_back,phase_diff_misaligned_2,truncated_zeros_forward];
%     coherence_extended_2(i_beam,:)          = (horzcat(2,truncated_zeros_back,coherence_misaligned_2,truncated_zeros_forward));
%     coherence_mask_extended_2(i_beam,:)     = (horzcat(2,truncated_zeros_back,coherence_mask_misaligned_2,truncated_zeros_forward));
    %% Reapply align the beams
    beams_rng_cmpr_extended_2(i_beam,:)     = circshift(beams_rng_cmpr_extended_2(i_beam,:),[0,L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
    beams_rng_cmprIQ_extended_2(i_beam,:)   = circshift(beams_rng_cmprIQ_extended_2(i_beam,:),[0,L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
    phase_diff_extended_2(i_beam,:)         = circshift(phase_diff_extended_2(i_beam,:),[0,L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
%     coherence_extended_2(i_beam,:)          = circshift(coherence_extended_2(i_beam,:),[0,L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
%     coherence_mask_extended_2(i_beam,:)     = circshift(coherence_mask_extended_2(i_beam,:),[0,L1BS.shift_coarse(i_beam)*cnf.zp_fact_range]);
end

end


%plot_extended_stacks(L1BS,beams_rng_cmpr_extended);

L1BS.beams_rng_cmpr=[];     L1BS.beams_rng_cmpr = beams_rng_cmpr_extended;
L1BS.beams_rng_cmprIQ=[];   L1BS.beams_rng_cmprIQ = beams_rng_cmprIQ_extended;

if(strcmp(cnf.processing_mode,'SIN'))&& strcmp(cnf.processing_mode,'SIN')
    L1BS.beams_rng_cmpr_2=[];   L1BS.beams_rng_cmpr_2 = beams_rng_cmpr_extended_2;
    L1BS.beams_rng_cmprIQ_2=[]; L1BS.beams_rng_cmprIQ_2 = beams_rng_cmprIQ_extended_2;
    L1BS.phase_diff=[];         L1BS.phase_diff = phase_diff_extended_2;
    L1BS.coherence=[];          L1BS.coherence = NaN(size(L1BS.beams_rng_cmpr_2));
    L1BS.coherence_mask=[];     L1BS.coherence_mask = NaN(size(L1BS.beams_rng_cmpr_2));
    stack_power = NaN(size(L1BS.beams_rng_cmpr_2));
    a1_power    = NaN(size(L1BS.beams_rng_cmpr_2));
    a2_power    = NaN(size(L1BS.beams_rng_cmpr_2));
    
    % substacking
    for i_beam=1:L1BS.N_beams_stack
        beam_margin_l=2;
        beam_margin_r=2;
        if(i_beam-beam_margin_l<1)
            continue;
%             beam_margin_l=i_beam-1;
        elseif((i_beam+beam_margin_r)>L1BS.N_beams_stack)
            continue;
%             beam_margin_r=L1BS.N_beams_stack-i_beam;
        end
        beams = (i_beam-beam_margin_l):(i_beam+beam_margin_r);
%         for i_sample = sample_start(i_beam):sample_start(i_beam)+ chd.N_samples_sar*cnf.zp_fact_range
%             range_margin_l=2;
%             range_margin_r=2;
%             if(i_sample-range_margin_l<sample_start(i_beam))
% %                 range_margin_l=0;
%                   continue;
%             elseif((i_sample+range_margin_r)>sample_start(i_beam)+ chd.N_samples_sar*cnf.zp_fact_range)
% %                 range_margin_r=0;
%                   continue;
%             end
%             samples = i_sample-range_margin_l:i_sample+range_margin_r;
%             stack_power(i_beam,i_sample) = nanmean(nanmean(L1BS.beams_rng_cmprIQ(beams,samples).*conj(L1BS.beams_rng_cmprIQ_2(beams,samples))));
%             a1_power(i_beam,i_sample) = nanmean(nanmean(L1BS.beams_rng_cmpr(beams,samples)));
%             a2_power(i_beam,i_sample) = nanmean(nanmean(L1BS.beams_rng_cmpr_2(beams,samples)));
%             if(L1BS.Gap_flag(i_beam))
%                 L1BS.coherence(i_beam,i_sample) = (abs(stack_power(i_beam,i_sample)))./sqrt(a1_power(i_beam,i_sample).*a2_power(i_beam,i_sample));
%                 if(L1BS.coherence(i_beam,i_sample)>cnf.coherence_threshold)
%                     L1BS.coherence_mask(i_beam,i_sample) = 1; % to be checked
%                 end
%             end
            stack_power(i_beam,:) = (nanmean(L1BS.beams_rng_cmprIQ(beams,:).*conj(L1BS.beams_rng_cmprIQ_2(beams,:))));
            a1_power(i_beam,:) = (nanmean(L1BS.beams_rng_cmpr(beams,:)));

            a2_power(i_beam,:) = (nanmean(L1BS.beams_rng_cmpr_2(beams,:)));
            if(L1BS.Gap_flag(i_beam))
                L1BS.coherence(i_beam,:) = (abs(stack_power(i_beam,:)))./sqrt(a1_power(i_beam,:).*a2_power(i_beam,:));
                L1BS.coherence_mask(i_beam,find(L1BS.coherence(i_beam,:)>cnf.coherence_threshold)) = 1; % to be checked
                
            end

%         end
    end

end










end