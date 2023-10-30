function [fit_res]=threshold_retracker(data,cnf_p)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code runs a simple threshold retracker: based on determining the
% leading edge as the % of the peak
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 01/07/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -data    =  data input structure for L2 processing
%       -cnf_p = configuration parameters structure for L2 processing
%      OPTIONAL
%       
% OUTPUT:
%       -fit_res        =   structure of fitting results {Epoch,sigma0}
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% - Need to optimize the code avoiding so many different data structures
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0:
% V1.1: 
%% --------------------- SEED GENERATION based on Pablo's Approach --------
% TBC and TBD

%% --------------------- RUN DETECTOR -------------------------------------
%--------------------------------------------------------------------------
% -------------------------------------------------------------------------
% IF MASK - Small aliasing neglection
% Due to the fact that the IFMask is not perfect and it does not go to zero
% in the stop band some aliasing is introduced in the first few samples.
% It' been agreed among the community that these samples ~12 first and last
% shall not be considered for the fitting of the waveform
% -------------------------------------------------------------------------
if strcmp(cnf_p.mode, 'SAR') || cnf_p.mode == 2 %RAW
    data.HRM.power_wav_filtered    =   data.HRM.power_wav(cnf_p.IFmask_N*cnf_p.ZP+1:end-(cnf_p.IFmask_N*cnf_p.ZP),:)';
elseif  cnf_p.mode == 3 %RMC
    data.HRM.power_wav_filtered    =   data.HRM.power_wav(cnf_p.IFmask_N*cnf_p.ZP+1:cnf_p.Nrmc,:)';
elseif strcmp(cnf_p.mode, 'SARin')
    data.SIN.power_wav_filtered    =   data.SIN.power_wav(cnf_p.IFmask_N*cnf_p.ZP+1:end-(cnf_p.IFmask_N*cnf_p.ZP),:)';
else
    error('not a valid mode');
end
%--------------------------------------------------------------------------
%-------------------DETECTOR-BASED ----------------------------------------
%--------------------------------------------------------------------------
[peak_pow,idx_max_peak]=max(data.HRM.power_wav_filtered,[],2);
for m = 1:data.N_records
    % added EM: 20.09.2016
    %checking waveform valid due to L1B processing related to noise
    %filtering outer beams
    if ~any(~isnan(data.HRM.power_wav_filtered(m,:)))
        fit_res.Epoch(m)    =   NaN;
        fit_res.Pu(m)       =   NaN;                
        fit_res.flag_validity_L1B(m) = 0; %due to L1B processing
        continue;
    end
    
    %dumm=find(find(data.HRM.power_wav_filtered(m,:)<=cnf_p.th_retracker.percentage_peak/100.0*peak_pow(m))<idx_max_peak(m), 1, 'last' );    
    dumm=find(data.HRM.power_wav_filtered(m,:)>cnf_p.th_retracker.percentage_peak/100.0*peak_pow(m), 1);
    if ~isempty(dumm)
        %position by simple linear interpolation
        if dumm>1
            idx_leading_edge_ISD=dumm-1+((cnf_p.th_retracker.percentage_peak/100.0*peak_pow(m)-data.HRM.power_wav_filtered(m,dumm-1))/(data.HRM.power_wav_filtered(m,dumm)-data.HRM.power_wav_filtered(m,dumm-1)));
        else
            idx_leading_edge_ISD=dumm;
        end
    else
        %if there is no leading edge or the waveform has displaced that
        %much from the window to the left select the peak as leading
        %edge
        idx_leading_edge_ISD=idx_max_peak(m);
    end
    fit_res.Epoch(m)    =   idx_leading_edge_ISD + cnf_p.IFmask_N*cnf_p.ZP - 1;
    fit_res.Pu(m)       =   10*log10(peak_pow(m)); %dB
    fit_res.flag_validity_L1B(m) = 1;
end

end

