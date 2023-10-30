function stack_com = stack_com_mask(stack, look, nf_p, cnf_p)

% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% The tracker on board limits the acquisition of data to the observation
% window with differs in meters from mode to mode in CryoSat. Eg. 60 m in
% SAR mode versus 240m in SARin. To build the stack we need to account for
% this limitation and we must not assume that after range cell migration
% data is infinit. Instead it is finit and zero values are given to the
% tails of those looks acquired at the edges of the observation window 
% -------------------------------------------------------------------------
% 
% Author:           Cristina Martin-Puig / isardSAT
%
% Reviewer:         Cristina Martin-Puig / isardSAT
%
% Last revision:    Cristina Martin-Puig / isardSAT v6 16/10/2013
%
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% INPUT:
%       stack     =   stack including infinite acquisition
%       Look        =   Look index
%       nf_p        =   Non-fit parameters
%       cnf_p       =   Configuration parameters
%
% OUTPUT:
%       stack_com =   stack compensated for finite acquisition
% -------------------------------------------------------------------------
    
%     kl              =   -(((nf_p.Lx/nf_p.Ly)*look).^2); % shift to be applied to begining of waveform in [mm]
%     kl_s            =   kl./nf_p.Lz;    % shift to be applied to begining of waveform in [range bin numbers]
%     stack_com       =   stack;
%     
%     %plot(nf_p.Doppler_mask)
%     
%     if cnf_p.L1proc == 1
%         if cnf_p.mode ==3
%             nf_p.Doppler_mask   =   nf_p.Doppler_mask;
%         else
%             nf_p.Doppler_mask   =   ((ceil(nf_p.Doppler_mask/cnf_p.ZP)-1)+1) * cnf_p.ZP;
%         end
%     elseif cnf_p.L1proc == 3
%         nf_p.Doppler_mask    = (nf_p.Doppler_mask+1) * cnf_p.ZP;
%     end
%     % added by EM 01.02.2016 to consider or not zeros in the averaging
%     if cnf_p.use_zeros_cnf == 1
%         value_mask=1*0;
%     else
%         value_mask=NaN;
%     end
    for m = 1:length(nf_p.Doppler_mask)
        if nf_p.Doppler_mask(m) == cnf_p.ZP
            stack_com(m,1:end) = value_mask;
        elseif nf_p.Doppler_mask(m) < size(stack,2)
            stack_com(m,nf_p.Doppler_mask(m)+1:end) = value_mask;      
        end
    end
%     for m = 1:length(nf_p.Doppler_mask)
%         if nf_p.Doppler_mask(m) == cnf_p.ZP
%             stack_com(m,1:end) = 1*0;
%         elseif nf_p.Doppler_mask(m) < size(stack,2)
%             stack_com(m,nf_p.Doppler_mask(m)+1:end) = 1*0;      
%         end
%     end
end
