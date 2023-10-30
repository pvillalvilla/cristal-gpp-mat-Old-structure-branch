function [fit_res]=OCOG_retracker(scaled_waveforms,L1B_alt,cnf_L2)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code runs an OCOG retracker based on Technical note CryoVal-LI
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Monica Roca / isardSAT
%
% Last revision:    Albert Garcia / isardSAT V1 01/07/2016
%                   Juan Pedro López Zaragoza / isardSAT V2 10/01/2023: Adapted to read FF data
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -L1B    =  data input structure for L2 processing
%       -cnf_L2 = configuration parameters structure for L2 processing
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
% - 
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0:
% V1.1: Updated struct and variable names
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

%--------------------------------------------------------------------------
%-------------------DETECTOR-BASED ----------------------------------------
%--------------------------------------------------------------------------
n1=ceil(cnf_L2.OCOG_retracker.n1/cnf_L2.ZP); %first sample to be used
n2=ceil(cnf_L2.OCOG_retracker.n2/cnf_L2.ZP); %last sample to be used
offset=cnf_L2.OCOG_retracker.offset/cnf_L2.ZP; %offset
vec_range=n1:n2;
for m = 1:length(L1B_alt)
    % added EM: 20.09.2016
    %checking waveform valid due to L1B processing related to noise
    %filtering outer beams
    if ~any(~isnan(scaled_waveforms(m,:)))
        fit_res.Epoch(m)    =   NaN;
        fit_res.Pu(m)       =   NaN;        
        fit_res.COG(m)      =   NaN;
        fit_res.A(m)        =   NaN;
        fit_res.W(m)        =   NaN;
        fit_res.flag_validity_L1B(m) = 0; %due to L1B processing
        continue;
    end
   
    switch cnf_L2.OCOG_retracker.param_comp_method
        case 0
            % using the squares of the power waveform a la Frappart
            %--------------------- COG ---------------------------------------------
            fit_res.COG(m)=(sum(vec_range.*(scaled_waveforms(m,n1:n2)).^2)/sum((scaled_waveforms(m,n1:n2)).^2));
            
            %-------------- Width --------------------------------------------------
            fit_res.W(m)= (sum((scaled_waveforms(m,n1:n2)).^2)).^2/sum((scaled_waveforms(m,n1:n2)).^4);
        case 1
             % using the power waveform a la Wingham definition
            %--------------------- COG ---------------------------------------------
            fit_res.COG(m)=(sum(vec_range.*(scaled_waveforms(m,n1:n2)))/sum((scaled_waveforms(m,n1:n2))));
            
            %-------------- Width --------------------------------------------------
            fit_res.W(m)= (sum(scaled_waveforms(m,n1:n2))).^2/sum((scaled_waveforms(m,n1:n2)).^2);
    end
   
   
   %----------- Computing the epoch ---------------------------------------
   switch cnf_L2.OCOG_retracker.implementation_method
       case 0
           %----------- thresholding the amplitude --------------------------------
           % amplitude computed as Frappart: from n1 to n2 as for COG
           switch cnf_L2.OCOG_retracker.param_comp_method
               case 0
                   fit_res.A(m) = sqrt(sum((scaled_waveforms(m,n1:n2)).^4)/sum((scaled_waveforms(m,n1:n2)).^2));
               case 1
                   fit_res.A(m) = (sum((scaled_waveforms(m,n1:n2)).^2)/sum((scaled_waveforms(m,n1:n2))));
           end
           threshold_pow=cnf_L2.OCOG_retracker.percentage_pow_OCOG/100.0*fit_res.A(m);
           dumm = find(scaled_waveforms(m,:) > threshold_pow,1);
           if ~isempty(dumm)
               %position by simple linear interpolation
               if dumm>1
                   idx_leading_edge_ISD=offset+dumm-1+((threshold_pow-scaled_waveforms(m,dumm-1))/(scaled_waveforms(m,dumm)-scaled_waveforms(m,dumm-1)));
               else
                   idx_leading_edge_ISD=offset+dumm;
               end
           else
               %if there is no leading edge or the waveform has displaced that
               %much from the window to the left select the peak as leading
               %edge
               [~,idx_max_peak]=max(scaled_waveforms(m,:),[],2);
               idx_leading_edge_ISD=idx_max_peak;
           end          
       case 1
           %----------- thresholding the amplitude --------------------------------
           % amplitude computed within the window centered ot OCOG
           n1_bis=round(fit_res.COG(m)-fit_res.W(m)/2);
           n2_bis=round(fit_res.COG(m)+fit_res.W(m)/2);
           switch cnf_L2.OCOG_retracker.param_comp_method
               case 0
                   fit_res.A(m) = sqrt(sum((scaled_waveforms(m,n1_bis:n2_bis)).^4)/sum((scaled_waveforms(m,n1_bis:n2_bis)).^2));
               case 1
                   fit_res.A(m) = (sum((scaled_waveforms(m,n1_bis:n2_bis)).^2)/sum((scaled_waveforms(m,n1_bis:n2_bis))));
           end
           threshold_pow=cnf_L2.OCOG_retracker.percentage_pow_OCOG/100.0*fit_res.A(m);
           dumm = find(scaled_waveforms(m,:) > threshold_pow,1);
           if ~isempty(dumm)
               %position by simple linear interpolation
               if dumm>1
                   idx_leading_edge_ISD=offset+dumm-1+((threshold_pow-scaled_waveforms(m,dumm-1))/(scaled_waveforms(m,dumm)-scaled_waveforms(m,dumm-1)));
               else
                   idx_leading_edge_ISD=offset+dumm;
               end
           else
               %if there is no leading edge or the waveform has displaced that
               %much from the window to the left select the peak as leading
               %edge
               [~,idx_max_peak]=max(scaled_waveforms(m,:),[],2);
               idx_leading_edge_ISD=idx_max_peak;
           end

       case 2
           %----------- Using OCOG only --------------------------------
           %epoch=offset+(COG-W/2)
           switch cnf_L2.OCOG_retracker.param_comp_method
               case 0
                   fit_res.A(m) = sqrt(sum((scaled_waveforms(m,n1:n2)).^4)/sum((scaled_waveforms(m,n1:n2)).^2));
               case 1
                   fit_res.A(m) = (sum((scaled_waveforms(m,n1:n2)).^2)/sum((scaled_waveforms(m,n1:n2))));
           end
           idx_leading_edge_ISD=offset+(fit_res.COG(m)-fit_res.W(m)/2);
   end
   fit_res.Epoch(m)    =   idx_leading_edge_ISD + cnf_L2.IFmask_N*cnf_L2.ZP - 1;
   fit_res.Pu(m)       =   10*log10(fit_res.A(m)); %dB using the amplitude of the OCOG to be discussed
   fit_res.A(m)=10*log10(fit_res.A(m));
   fit_res.flag_validity_L1B(m) = 1;
   fit_res.COG(m)= fit_res.COG(m)/cnf_L2.ZP -1/cnf_L2.ZP+1;
end

end

