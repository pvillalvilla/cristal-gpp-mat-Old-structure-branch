% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code runs an TCOG retracker based on Technical note CryoVal-LI
% -------------------------------------------------------------------------
% 
% Author:           Albert Mondejar / isardSAT
%
% Reviewer:         MÃ²nica Roca / isardSAT
%
% Last revision:    Albert Mondejar / isardSAT V1 18/04/2019
%                   Juan Pedro López Zaragoza / isardSAT V2 10/01/2023: Adapted to read FF data
%                                               
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [fit_res] = TCOG_retracker(scaled_waveforms, L1B_alt, cnf_L2, chd)
cnf_L2.power_threshold = 0.8;
n1=ceil(cnf_L2.OCOG_retracker.n1/cnf_L2.ZP); %first sample to be used
n2=ceil(cnf_L2.OCOG_retracker.n2/cnf_L2.ZP); %last sample to be used


for m = 1:length(L1B_alt)
   
    if ~any(~isnan(scaled_waveforms(m,:)))
        fit_res.Epoch(m)    =   NaN;
        fit_res.Pu(m)       =   NaN;        
        fit_res.COG(m)      =   NaN;
        fit_res.A(m)        =   NaN;
        fit_res.W(m)        =   NaN;
        fit_res.flag_validity_L1B(m) = 0; %due to L1B processing
        continue;
    end
    if(strcmp(chd.simulation,'OC'))
        fit_res.A(m)    = sqrt(sum((scaled_waveforms(m,n1:n2)).^4)/sum((scaled_waveforms(m,n1:n2)).^2));
    
    else
        fit_res.A(m)    = max(scaled_waveforms(m,n1:n2));
    end
    above_threshold=find(scaled_waveforms(m,2:end)>cnf_L2.power_threshold*fit_res.A(m));
    
    init_cross      = min(above_threshold)+1;
    
    fit_res.Epoch(m) = init_cross -1 + (cnf_L2.power_threshold*fit_res.A(m)-scaled_waveforms(m,init_cross-1))/(scaled_waveforms(m,init_cross)-scaled_waveforms(m,init_cross-1));
    fit_res.Epoch(m) = fit_res.Epoch(m)/cnf_L2.ZP - 1/cnf_L2.ZP +1;
    
    
end
end
