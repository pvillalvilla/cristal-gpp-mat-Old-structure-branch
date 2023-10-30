function [ml_waveform,stack,stack_real] = ml_wave_gen(x,fit_params,nf_p,cnf_p,looks,func_f0,func_f1)%varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% Generation of multi-looked waveform in IEEE TGRS "SAR Altimeter Backscattered Waveform Model" 
% DOI:10.1109/TGRS.2014.23330423
% -------------------------------------------------------------------------
% 
% Author:           Cristina Martin-Puig / isardSAT
%
% Reviewer:         Cristina Martin-Puig / isardSAT
%
% Last revision:    Cristina Martin-Puig / isardSAT V9 4/7/2014
%
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       x           =   Range bin index 
%       fit_p       =   Parameters to be fitted 
%       nf_p        =   parameters of need for the waveform generation, but not to be
%                       fitted
%       cnf_p       =   configuration parameters of L2 processor
% OUTPUT:
%       ml_waveform =   multi-looked waveform
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - sl_wave_gen: in charge of genrating the single look waveform based on
% the model proposed by Chris in IEEE TGRS "SAR Altimeter Backscattered Waveform Model" 
% DOI:10.1109/TGRS.2014.23330423 
% - stack_com_mask
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Versions control:

%added by EM 01.03.2016
%$$$$$$$$$$$$$$$$$$$$$$$$$ Parse optional input parameters $$$$$$$$$$$$$$$$
% p = inputParser;
% defaultff0=0.0;
% addOptional(p,'func_f0',defaultff0,@isnumeric);
% addOptional(p,'func_f1',defaultff0,@isnumeric);
% parse(p,varargin{:});
% func_f0=p.Results.func_f0;
% func_f1=p.Results.func_f1;
% %added by EM 15.03.2016
% 
% clear p
if nargin < 7 
    func_f1=0; % initialize altitude to 0;
end
if nargin <6
    func_f0=0; % initialize altitude to 0;
end

% switch cnf_p.mission
%     % -------------------------------------------------------------------------
%     % Build STACK
%     % -------------------------------------------------------------------------
%     case 'CS2'
%         if strcmp(cnf_p.mode, 'SAR')
%             % if the number of looks is less than maximum
%             if nf_p.Neff < (nf_p.Npulses*nf_p.Nbcycle)
%                 % we need to find out which is the number of looks excluded
%                 N_Tlooks    =   round (nf_p.Neff/nf_p.Nbcycle);
%                 % if it is equal to all the illumination this is much easier
%             elseif nf.p.Neff == (nf_p.Npulses*nf_p.Nbcycle)
%                 N_Tlooks    =   nf_p.Npulses;
%             else
%                 error('invalid Neff value')
%             end
%         elseif strcmp(cnf_p.mode, 'SARin')
%             if nf_p.Neff < nf_p.Npulses
%                 % we need to find out which is the number of looks excluded
%                 N_Tlooks    =   nf_p.Neff;
%                 % if it is equal to all the illumination this is much easier
%             elseif nf.p.Neff == nf_p.Npulses*nf_p.Nbcycle
%                 N_Tlooks    =   nf_p.Npulses;
%             else
%                 error('invalid Neff value')
%             end
%         else
%             error('not a valid mode');
%         end
%         
%         N_Tlooks        =   N_Tlooks-1;
%         %looks           =   -round(N_Tlooks/2):round(N_Tlooks/2);
%         looks           =   (-round(N_Tlooks/2):round(N_Tlooks/2))/(nf_p.Npulses);
%         stack           =   zeros(length(looks),length(x));
%         
%         if cnf_p.L1proc == 2
%             for l = 1:length(looks)
%                 % stack (l,:)     =   fit_params(3) * sl_wave_gen (x,looks(l),[fit_params(1) fit_params(2)],nf_p,cnf_p)+ nf_p.ThN;  % note this is erroneous but we need to emulate CNES
%                 % stack (l,:)     =   fit_params(3) *sl_wave_gen (x,looks(l),[fit_params(1) fit_params(2)],nf_p,cnf_p) + nf_p.sinc2G;
%                 stack (l,:)     =   fit_params(3) *sl_wave_gen (x,looks(l),[fit_params(1) fit_params(2)],nf_p,cnf_p);
%             end
%         else
%             for l = 1:length(looks)
%                 stack (l,:)     =   sl_wave_gen (x,looks(l),[fit_params(1) fit_params(2)],nf_p,cnf_p) + nf_p.ThN;
%             end
%         end
%     case 'JCS'
        stack       =   zeros(length(looks),length(x));
        
        for l = 1:length(looks)
            %nf_p.ThN = 0;
            %modified by EM 01.03.2016 / modified by EM 15.03.2016
            if cnf_p.lut_flag
                switch cnf_p.power_wfm_model
                    case 'simple'
                        %a   = sl_wave_gen(x,looks(l),fit_params,nf_p,cnf_p,'func_f0',func_f0);
                        a   = sl_wave_gen(x,looks(l),fit_params,nf_p,cnf_p,func_f0);
                    case 'complete'
                        %a   = sl_wave_gen(x,looks(l),fit_params,nf_p,cnf_p,'func_f0',func_f0,'func_f1',func_f1);
                        a   = sl_wave_gen(x,looks(l),fit_params,nf_p,cnf_p,func_f0,func_f1);
                end
            else
                a   = sl_wave_gen(x,looks(l),fit_params,nf_p,cnf_p);
            end
            %added by EM 16.03.2016
            %including noise fitting estimation window
            
            stack(l,:)=a;
            %                 % added by EM 08.02.2016
            %                 switch cnf_p.multilook_option
            %                     case 'Cris_old'
            %                         stack(l,:)= fit_params(3).*a+nf_p.ThN(l);
            %                     case 'NormML_plus_noise'
            %                         stack(l,:)= a;
            %                 end
        end
%     otherwise
%         error('not a valid mission')
% end
    
%     if cnf_p.fit_noise
%         switch cnf_p.Thn_method
%             case 'ML'
%                 idx_ledge=[];
%                 temp_noise_thr=cnf_p.threshold_noise;
%                 while isempty(idx_ledge)
%                     idx_ledge=find(diff(mean(stack,1))>temp_noise_thr);
%                     %                 if isempty(idx_ledge)
%                     %                     'stop'
%                     %                 end
%                     temp_noise_thr=temp_noise_thr/2;
%                 end
%                 nf_p.ThN        =   mean(input_data(1:idx_ledge(1)-1))*ones(1,length(looks));
%                 %idx_ledge(1)
%         end
%     end
switch cnf_p.multilook_option
    case 'Cris_old'
        stack= fit_params(3).*stack+(nf_p.ThN).'*ones(1,length(x));
end
%fit_params(3)
% -------------------------------------------------------------------------
% compensate for Range Cell migration
% -------------------------------------------------------------------------
stack_real     =   stack.*nf_p.Doppler_mask; %stack_com_mask(stack, looks, nf_p, cnf_p);


% ml_waveform = zeros(1,size(stack_real,2));
% if cnf_p.L1proc == 2 % CNES
%     %ml_waveform     =   sum(stack_real)/length(looks);
%     ml_waveform     =   (sum(stack_real)/length(looks)) + nf_p.ThN;
% else
    % Modified by EM 01.03.2016
    %         for i_sample =1:size(stack_real,2)
    %             %ml_waveform(i_sample) = mean(stack_real(isfinite(stack_real(:,i_sample)),i_sample));
    %             %modified by EM to remove the impact of the zero-regions (masked ones)
    %             %ml_waveform(i_sample) = mean(stack_real(isfinite(stack_real(:,i_sample)) & stack_real(:,i_sample)~=0,i_sample));
    %             % added by EM 05.02.2016
    %             ml_waveform(i_sample) = mean(stack_real(isfinite(stack_real(1:length(looks),i_sample)),i_sample));
    %
    %         end
    
    ml_waveform = nanmean(stack_real(1:length(looks),:),1);
    
    % added by EM 08.02.2016
    switch cnf_p.multilook_option
        case 'NormML_plus_noise'
            max_ML_signal=max(ml_waveform);
            %ml_waveform = fit_params(3).*(ml_waveform/max_ML_signal+nf_p.ThN(1));
            ml_waveform = (ml_waveform/max_ML_signal+nf_p.ThN(1));
            ml_waveform = ml_waveform/max(ml_waveform);
            %stack=fit_params(3).*(stack+nf_p.ThN(1)/size(stack,1)*max_ML_signal);
            % try not to add noise to the samples forced to zero by the
            % mask or due to the fact that there are no data for other
            % beams as not
            %stack_real(stack_real~=0)=fit_params(3).*(stack_real(stack_real~=0)+nf_p.ThN(1)/size(stack,1)*max_ML_signal);
    end
    %         if ~cnf_p.Doppler
    %             %ml_waveform_old     =   fit_params(3)*sum(stack_real(101:end-100,:))/(length(looks)-200);
    % %             N_average   =   zeros(1, size(stack_real,2));
    % %                 for m = 1:size(stack_real,2)
    % %                     N_average(m)   =   length(find(stack_real(:,m)>0));
    % %                 end
    %             ml_waveform     =   sum(stack_real)/(length(looks)-1);%./N_average;% change this!!!
    %         else
    %             ml_waveform     =   sum(stack_real)/(length(looks)-1);
    %         end
%end

end

