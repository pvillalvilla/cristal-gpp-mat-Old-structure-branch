% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% CryoSat 2 calibration over transponders
% 
% This code implements the algorithm as described in the
% ISARD_ESA_CR2_TRP_CAL_DPM_030 2.b of 26/05/2011
%
% ---------------------------------------------------------
% ECHOFIT_CRY: Code that finds the parameters of the Gaussian that better fits
% with the original signal input
% 
% Calling
%   E =echofit_CRY(x0, wvf)
%
% Inputs
%   x0 : initial gaussian parameters
%   wvf : input signal, to be fitted with a gaussian
%
% Output
%   E : rms
%
% ----------------------------------------------------------
% 
% Author:   Mercedes Reche / Pildo Labs
% Reviewer: Mònica Roca / isardSAT
%
%
%
% $Id: echofit_CRY.m 139 2005-12-20 07:48:18Z merche $
%
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function E =echofit_CRY(x0, wvf)

x = 1:length(wvf); 
y = wvf; 
sum=0;
for j=1:length(x)     
       sum=sum+ ( x0(1)*exp(-(x(j)-x0(3)).^2/x0(2).^2)-y(j)).^2;    
end

E=sum;



