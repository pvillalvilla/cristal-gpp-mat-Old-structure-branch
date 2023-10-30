% RANGE COMPRESSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% JasonCS 
% This code implements the RANGE COMPRESSION 
%
% ---------------------------------------------------------
% Objective: Perform the SECONDARY RANGE COMPRESSION
% 
% INPUTs: 
%
%     Name                    Description  
%     wvfm_input:              Input waveform to be transformed in range time
%     w_rg:                   Range window for side-lobes control
%     
% OUTPUTs:  
%     Name                    Description
%     wvfm_out
%     
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (14/11/2017)

function [wvfm_out]...
      = RC_processing (wvfm_input,w_rg,zp_fact_range_cnf)

s = size(wvfm_input);
N_samples_rg = s(2);
N_samples_az = s(1);


%% 1. APPLICATION OF THE RANGE WINDOW
%Input data is in the range frequency and (azimuth time or Doppler domain)
non_centered_spectra=fftshift(wvfm_input.*(ones(N_samples_az,1)*w_rg),2);

%% 2. TRANSFORMATION TO TIME DOMAIN WITH ZERO-PADDING
wvfm_out = fft([non_centered_spectra(:,1:N_samples_rg/2),...
                       zeros(N_samples_az,(zp_fact_range_cnf-1)*N_samples_rg),...
                       non_centered_spectra(:,N_samples_rg/2+1:N_samples_rg)],...
                       N_samples_rg * zp_fact_range_cnf,2);

progressbar([],[], ...
    100/100,[]);
end