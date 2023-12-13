function fit_res = fitting_EM_bis(data, cnf_p, cst,  chd, nf_p, ini_p,varargin)

% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% Routine to perform the fitting of the waveforms
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%                   Cristina Martin-Puig / isardSAT
%                   
%
% Reviewer:        Monica Roca / isardSAT
%
% Last revision:    Albert Garcia / isardSAT 2019/06/10
%
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% INPUT:
%   MANDATORY:
%       data        =   SAR data structure 
%       nf_p        =   Non-fit parameters of the fitting
%       cnf_p       =   Configuration parameters of the L2 processor
%       ini_p       =   initialization parameters of need for the fitting
%   OPTIONAL:
%       LUT_f0_file =   file containing the look up table for the f0
%                       function parametrization (full path)
%       LUT_f1_file =   file containing the look up table for the f1
%                       function parametrization (full path)
%
% OUTPUT:
%       fit_res     =   structure containing fitting results
%           .Epoch  =   Epoch not compensated for IF_Mask
%           .Hs     =   significant wave height
%           .Pu     =   Pu if to be fitted
%           .GOF    =   Goodness of fitting in %
%           .flag   =   fitting flag output = the same as for lsqcurve fit
%                       flags
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - ml_wave_gen_EM: genrates the multilook waveform based on Chris Ray
% analytical model in IEEE TGRS "SAR Altimeter Backscattered Waveform Model" 
% DOI:10.1109/TGRS.2014.23330423
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% - 
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% V7.1: Upgrading the look indexation using different methods based on
% normalized indexation (first rough approach), Doppler freq info, look
% angle information. The latter two methods have the exact and approximate approaches
% with the possibility to use the whole angular information (beam angles
% and look angles respectively) from the stack itself or using the start
% and stop angles. Vector of look indexation is passed to the ml waveform
% function.
% V7.2: Option to perform the estimation of the noise per look provided the
% whole stack is accessible to the L2 processing (to date such option will be 
% considered as no big differences have been observed)
% V7.3: It has been incorportaed the option to construct the Power waveform
% model as simple (using only f0 function: equation (43) in paper Chris
% Ray) or the complete option (including f0 and f1, as in equation (42)).
% It has been considered the option to implement sucgh f0 and f1 using
% their expansions in termns of bessel functions or using Look Up Tables to
% accelerate the fitting.
% V7.4: It has been included the option to perform the LSE minimization
% either using the lsqcurvefit or the fminsearch (minimizing the MSE)
% V7.5: A pre-processing optional stage has been incorporated to perform a
% rough threshold-based retracker to have a first initial guess of the epoch
%, which can be used instead of the one provided by the user 
%(some issues encountered in the fitting procedure depending on 
% the initial epoch guess provided by the user --> local minimums). In such
% pre-processing is it also possible to perform an adaptive estimation of
% the noise thermal noise based on the derivative of the real waveform
% without using a fixed set of initial sample and window size for the
% estimation of the noise floor.
% V7.6: Construction of the stack using matrix notation instead of for
% loops to accelerate the formation of the multilook waveform. Computation
% and generation of the stack mask based on the stack_mask vector provided
% as input in the L1B.
% V7.7: Provide an average estimation of the initial parameters for the
% next surface based on a sliding window or as an accumulated average value
% from the first surface for both the Hs and roughness. Specific filtering
% in the average computation is made based on the mean and std.
% v7.8: Incorporating the option to perform a two-step fitting first on a
% ocean-like waveform using SWH and the second depending on the correlation
% threshold based on roughness fitting leads-like surface. Keep the best
% fitting among the two of them
%
% v7.9 changed struct and variable names



%% ---------------- Handling input variables ------------------------------
% if(nargin<4 || nargin>(4+4*2))
%     error('Wrong number of input parameters');   
% end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('LUT_f0_file',{''},@(x)ischar(x));
p.addParamValue('LUT_f1_file',{''},@(x)ischar(x));
p.addParamValue('path_Results',{''},@(x)ischar(x));
% p.addParamValue('L1B_filename',{''},@(x)ischar(x));
p.parse(varargin{:});
LUT_f0_file=char(p.Results.LUT_f0_file);
LUT_f1_file=char(p.Results.LUT_f1_file);
path_Results=char(p.Results.path_Results);
% L1B_filename=char(p.Results.L1B_filename);
clear p;


%% --------------- GLOBAL VARIABLES ---------------------------------------
%% -------------- LOADING THE LUTS -----------------------------------------
%--------------------------------------------------------------------------
% added by EM 01.03.2016 include the LUT load
% modified by EM 15.03.2016 for func_f1
if (cnf_p.lut_flag)
    [~,~,LUT_f0_file_ext]=fileparts(LUT_f0_file);
    switch LUT_f0_file_ext
        case '.mat'
            load(LUT_f0_file,'func_f0');
        otherwise
            error('No valid LUT f0 file format');
    end
    switch cnf_p.power_wfm_model
        case 'complete'
            [~,~,LUT_f1_file_ext]=fileparts(LUT_f1_file);
            switch LUT_f1_file_ext
                case '.mat'
                    load(LUT_f1_file,'func_f1');
                otherwise
                    error('No valid LUT f1 file format');
            end
        case 'simple'    
            %
        otherwise
            error('No valid power waveform model')
    end
end


% -------------------------------------------------------------------------
% NON-FIT MODEL PARAMETERS SINGLE VARIABLE
% read non-fit parameters which are single value
% -------------------------------------------------------------------------
    nf_params.waveskew  =   nf_p.waveskew;
    nf_params.EMbias    =   nf_p.EMbias;
    nf_params.rou       =   nf_p.rou;
    nf_params.Npulses   =  64;% nf_p.Npulses;
    nf_p.Npulses= 64 ;
    nf_params.Lz        =   nf_p.Lz;
    nf_params.alphag_a  =   nf_p.alphag_a;
    nf_params.alphag_r  =   nf_p.alphag_r;
    nf_params.Nbcycle   =   nf_p.Nbcycle;
    nf_params.bw_Rx     =   nf_p.bw_Rx;
    nf_params.Nsamples  =   nf_p.Nsamples;
    
    nf_params.A_s2Ga= nf_p.A_s2Ga;
    nf_params.A_s2Gr= nf_p.A_s2Gr;

% -------------------------------------------------------------------------
% FITTING SEEDS
% Initialize parameters to be fitting according to flags
% -------------------------------------------------------------------------
    %EM: 13.09.2016
    if cnf_p.two_step_fitting && cnf_p.rou_flag
        cnf_p.rou_flag=0; %ensure first the SWH fitting is performed
    end
    
    if cnf_p.rou_flag
        if strcmp(cnf_p.multilook_option,'Cris_old')
            fit_params_ini      =   [ini_p.Epoch ini_p.rou ini_p.Pu];
        else
            fit_params_ini      =   [ini_p.Epoch ini_p.rou];
        end
        
    else
        if strcmp(cnf_p.multilook_option,'Cris_old')
            fit_params_ini      =   [ini_p.Epoch ini_p.sigmaz ini_p.Pu];
        else
            fit_params_ini      =   [ini_p.Epoch ini_p.sigmaz];
        end
        
    end
%     %added by Em 16.03.2016 to include fitting optimization of thermal
%     %noise floor
%     if cnf_p.fit_noise
%         fit_params_ini=[fit_params_ini,ini_p.Thn_1st,ini_p.Thn_wdth];
%     end
    

% from the above it is clear that maximum three parameters are to be
% fitted. 
% (1)- Epoch --- is always fitted
% (2)- mss/Hs   --- roughtness or significant wave height are to be fitted.
% Either one or the other
% (3)- Pu   --- amplitude is to be fitted

% -------------------------------------------------------------------------
% IF MASK - Small aliasing neglection
% Due to the fact that the IFMask is not perfect and it does not go to zero
% in the stop band some aliasing is introduced in the first few samples.
% It' been agreed among the community that these samples ~12 first and last
% shall not be considered for the fitting of the waveform
% -------------------------------------------------------------------------

    if strcmp(cnf_p.mode, 'SAR') %%|| cnf_p.mode == 2 %RAW
        data.scaled_waveforms_filtered    =   data.scaled_waveforms(cnf_p.IFmask_N*cnf_p.ZP+1:end-(cnf_p.IFmask_N*cnf_p.ZP),:)';
    elseif  strcmp(cnf_p.mode, 'RMC')
        data.scaled_waveforms_filtered    =   data.scaled_waveforms(:,cnf_p.IFmask_N*cnf_p.ZP+1:cnf_p.Nrmc)';  
    elseif strcmp(cnf_p.mode, 'SARin')
        data.scaled_waveforms_filtered    =   data.scaled_waveforms(cnf_p.IFmask_N*cnf_p.ZP+1:end-(cnf_p.IFmask_N*cnf_p.ZP),:)';
    else
        error('not a valid mode');
    end
%--------------------------------------------------------------------------
% --------------- Define output structure ---------------------------------
fit_res.Epoch=zeros(1,length(data.lat));
fit_res.rou=zeros(1,length(data.lat));
fit_res.Hs=zeros(1,length(data.lat));
fit_res.Pu=zeros(1,length(data.lat));                     
fit_res.flag=-1.0*ones(1,length(data.lat));
fit_res.COR=zeros(1,length(data.lat));

accumulated=0; %keep an accumulated value of the Hs for 
%parameter initialization for next surface
m_rou=0; %keep control on how many surfaces enter the second step-fitting on roughness
m_fitted=0; %keep number of fitted waveforms (avoiding those where L1B waveform is not valid)

%due to missing initialization in L1B processor for check continuty
idx_pri_non_zero=find(data.pri,1);

%--------------------------------------------------------------------------
% ----------------- INITIALIZATION OF INDEX RANGE VECTOR ------------------
%--------------------------------------------------------------------------
switch cnf_p.range_index_method
    case 'conventional'
        range_index=1:nf_p.Nsamples;
    case 'resampling'
        range_index=interp((1:nf_p.Nsamples/cnf_p.ZP),cnf_p.ZP);
end

%ThN_estimate_accumulated=0.0;
%%
% -------------------------------------------------------------------------
% Least Square Fitting
% -------------------------------------------------------------------------
disp(strcat('Total # wvfms ',num2str(length(data.lat))))
    for m = 1:length(data.lat)
        %t_fitting_init=tic;
        % added EM: 20.09.2016
        % check if the formed ML waveform is all NaN due to noise filtering
        % outer looks
        if ~any(~isnan(data.scaled_waveforms_filtered(:,m)))
            fit_res.Epoch(m)    =   NaN;
            if cnf_p.rou_flag
                fit_res.rou(m)  =   NaN; 
            else
                fit_res.Hs(m)   =   NaN;
                if cnf_p.two_step_fitting
                    fit_res.rou(m)  =   NaN;
                end
            end
            fit_res.Pu(m)       =   NaN;
            fit_res.flag(m)     =   NaN;           
            fit_res.COR(m)      =   NaN; 
            fit_res.flag_validity_L1B(m) = 0; %due to L1B processing
            continue;
        end
        m_fitted=m_fitted+1;

%             fprintf('-------\n');
%             fprintf('- Fitting wave num = %d\n',m-1);

% -------------------------------------------------------------------------
% NON-FIT MODEL PARAMETERS DEFINED AS ARRAY
% read non-fit parameters which are single value
% -------------------------------------------------------------------------
            nf_params.h         =   nf_p.h(m);
            nf_params.xp        =   nf_p.xp(m);
            nf_params.yp        =   nf_p.yp(m);
            nf_params.alphax    =   nf_p.alphax(m);
            nf_params.alphay 	=   nf_p.alphay(m);
            nf_params.Lx        =   nf_p.Lx(m);
            nf_params.Ly        =   nf_p.Ly(m);
            nf_params.Lgamma    =   nf_p.Lgamma(m);
            nf_params.Neff      =   nf_p.Neff(m);
            


% -------------------------------------------------------------------------
% DEFINE data to fit 
% -------------------------------------------------------------------------
            if strcmp(cnf_p.mode, 'SAR') || strcmp(cnf_p.mode, 'RMC')
                max_power_wav       =   max(data.scaled_waveforms_filtered(:,m));
                fit_data            =   data.scaled_waveforms_filtered(:,m)/max_power_wav;          
            else
                max_power_wav       =   max(data.scaled_waveforms_filtered(:,m));
                fit_data            =   data.scaled_waveforms_filtered(:,m)/max_power_wav; 
            end
            
% Added EM 21.01.2016            
%--------------------------------------------------------------------------
% Define the look indexation
%--------------------------------------------------------------------------
       switch cnf_p.looks_index_method
            case 'Cris_old'
                N_Tlooks    =   nf_p.Neff(m);
                looks       =   -floor((N_Tlooks-1)/2):floor(N_Tlooks/2);                
           case 'Norm_index'
                N_Tlooks    =   nf_p.Neff(m);
                looks       =   -floor((N_Tlooks-1)/2):floor(N_Tlooks/2);
                looks       =   (looks)/(nf_p.Npulses); % changed by EM: Model Chris l was per burst; not for stack (l should be normalized accordingly) 
           case 'Doppler_freq'
           %Generate the look index according to the relationship between
           %beam angles and Doppler
                switch cnf_p.fd_method
                    case 'exact'
%                         looks=2.0/chd.wv_length.*data.vel_stack(m,1:nf_p.Neff(m)).*cos(data.HRM.beam_ang_stack(m,1:nf_p.Neff(m))).*...
%                             data.HRM.pri_stack(m,1:nf_p.Neff(m)).*nf_p.Npulses;
                    case 'approximate'
                        %compute the max and min beam angles from 
                        % Look and Doppler Angles information
                        beam_angles=...
                         linspace(pi/2+data.doppler_angle_start(m)-data.look_ang_start(m),...
                                  pi/2+data.doppler_angle_stop(m)-data.look_ang_stop(m),nf_p.Neff(m)).';
                        looks=2.0/chd.wv_length.*data.vel(m).*cos(beam_angles).*data.pri(m).*nf_p.Npulses;
                end 
           case 'Look_angle'
               if data.pri(m)==0
                   data.pri(m)=data.pri(idx_pri_non_zero);
               end
               delta_look_angle=asin(chd.wv_length./(data.pri(m).*2.0*...
                                    nf_p.Npulses*data.vel(m))); 
               switch cnf_p.look_ang_method
                    case 'exact'
%                         looks=data.HRM.look_ang_stack(m,1:nf_p.Neff(m))/delta_look_angle;
                    case 'approximate'
                        %compute the max and min beam angles from 
                        % Look and Doppler Angles information
                        looks=(linspace(data.look_ang_start(m),...
                            data.look_ang_stop(m),nf_p.Neff(m))/delta_look_angle).';
                        
                end 
               
       end 
       
       %-----------------------------------------------------------
       %---------- Create the mask (as a matrix)-------------------
       %-----------------------------------------------------------
       %generated once avoiding to create it within the ML waveform
       %function for computational load issues
       if isfield(data,'Doppler_mask')
           nf_params.Doppler_mask=ones(nf_params.Neff,length(range_index));
           %to be in-line with L1B-processing
           if cnf_p.use_zeros_cnf == 1
               value_mask=1*0;
           else
               value_mask=NaN;
           end
           for i_look=1:nf_params.Neff
               nf_params.Doppler_mask(i_look,data.Doppler_mask(i_look,m):end)=value_mask;
           end
       else
           %create the Doppler mask based on the geometry information
           %case of reading L1B data directly from ESA
           % TBD
       end
       
       
% added by EM 29.03.2016: pre-processing stage to estimate the initial
% epoch as mid-point leading edge (some percent of peak power) & for
% estimation of the floor noise using derivative of the data itself
% detecting the leading edge (not fixing the window size initially)
% -------------------------------------------------------------------------
% Pre-processing
% -------------------------------------------------------------------------
if cnf_p.pre_processing
    %used to estimate the initial value of the epoch 
    %& eventually if desired the noise floor
    [peak_pow,idx_max_peak]=max(fit_data);
    idx_leading_edge=find(find(fit_data<=cnf_p.percent_leading_edge/100.0*peak_pow)<idx_max_peak, 1, 'last' );
    if isempty(idx_leading_edge)
        %if there is no leading edge or the waveform has displaced that
        %much from the window to the left select the peak as leading
        %edge
        idx_leading_edge=idx_max_peak;
    end
    fit_params_ini(1)=range_index(idx_leading_edge);
    
    %noise floor estimation
    if cnf_p.fit_noise
        switch cnf_p.Thn_method
            case 'ML'
                idx_noise=[];
                temp_noise_thr=cnf_p.threshold_noise;
                iter_noise=1;
                while isempty(idx_noise) && iter_noise<cnf_p.max_iter_noise
                    idx_noise=find(abs(diff(fit_data))<=temp_noise_thr);
                    idx_noise=idx_noise(idx_noise<idx_leading_edge);
                    temp_noise_thr=temp_noise_thr*1.5;
                    iter_noise=iter_noise+1;
                end
                if iter_noise<cnf_p.max_iter_noise
                    nf_params.ThN        =   mean(fit_data(idx_noise))*ones(1,length(looks));
                    clear idx_noise tempo_noise_thr;
                else
                    %take an average mean value of the previous Thermal
                    %noise estimation and keep it as reference for this
                    %record
                    if m~=1
                        nf_params.ThN         = mean(nf_params.ThN)*ones(1,length(looks));
                    else
                        %!!!need to be reviewed
                        %peak at the beginning of the window and so no
                        %space left for estimation of the noise floor
                        nf_params.ThN         = min(fit_data)*ones(1,length(looks));
                    end
                end
        end
        %ThN_estimate_accumulated=ThN_estimate_accumulated+nf_params.ThN(1);
        %disp(nf_params.ThN(1))
    end
    
    %clear idx_leading_edge peak_pow idx_max_peak;
end
       
% -------------------------------------------------------------------------
% ThN
% -------------------------------------------------------------------------

            if cnf_p.Thn_flag && ~cnf_p.fit_noise
                  % added by EM 02/12/2015 generate a vector of thermal noise for the
                   % different looks
                   switch cnf_p.Thn_method
                       case 'ML'
                           nf_params.ThN        =   mean(fit_data(cnf_p.Thn_w_first:cnf_p.Thn_w_first+cnf_p.Thn_w_width-1))*ones(1,nf_p.Neff(m));
                       case 'SL'
                           dumm=squeeze(data.HRM.beams_rng_cmpr_mask(m,1:nf_p.Neff(m),:));
                           dumm(~isfinite(dumm))=NaN;
                           dumm(dumm==0)=NaN;
                           max_mean_dumm=max(mean(dumm,1,'omitnan'));
                           for i_beam=1:nf_p.Neff(m)
                               nf_params.ThN(i_beam)        =   mean(dumm(i_beam,cnf_p.Thn_w_first:cnf_p.Thn_w_first+cnf_p.Thn_w_width-1)/max_mean_dumm,'omitnan');
                           end
                           
                           clear dumm;
                   end
            end
            
% -------------------------------------------------------------------------
% DEFINE FITTING MODEL & CALL FITTING ROUTINE 
% -------------------------------------------------------------------------                        
            % -------------------------------------------------------------
            % --------- DEFINE FITTING FUNCTION MODEL ---------------------
            % -------------------------------------------------------------            
            if cnf_p.lut_flag
                switch cnf_p.power_wfm_model
                    case 'simple'
                        switch cnf_p.fitting_fun_type
                            case 'lsq'
                                mpfun               =   @(fit_params,x)ml_wave_gen_EM(x, fit_params, nf_params, cnf_p,cst, chd,looks,func_f0);
                            case 'fmin'
                                fminfun               =   @(fit_params)sum((ml_wave_gen_EM(range_index, fit_params, nf_params, cnf_p,cst, chd,looks,func_f0)-fit_data).^2);
                        end
                        
                    case 'complete'
                        switch cnf_p.fitting_fun_type
                            case 'lsq'
                                mpfun               =   @(fit_params,x)ml_wave_gen_EM(x, fit_params, nf_params, cnf_p,cst, chd,looks, func_f0,func_f1);
                            case 'fmin'
                                fminfun               =   @(fit_params)sum((ml_wave_gen_EM(range_index, fit_params, nf_params, cnf_p,cst, chd,looks,func_f0,func_f1)-fit_data).^2);
                        end
                end
            else
                switch cnf_p.fitting_fun_type
                    case 'lsq'
                        mpfun               =   @(fit_params,x)ml_wave_gen_EM(x, fit_params, nf_params, cnf_p,cst, chd,looks);
                    case 'fmin'
                        fminfun               =   @(fit_params)sum((ml_wave_gen_EM(range_index, fit_params, nf_params, cnf_p, cst, chd,looks)-fit_data).^2);
                end
            end
            
            % -------------------------------------------------------------
            % --------- DEFINE TYPE OF MINIMIZATION PROCEDURE -------------
            % -------------------------------------------------------------
            switch cnf_p.mode
                case {'SAR','SARin','SARIN','HR','RMC','RAW'}
                    switch cnf_p.fitting_fun_type
                        case 'lsq'
                            %[fit_params,~,~,flag]     =   lsqcurvefit (mpfun,fit_params_ini,range_index,fit_data(range_index).',cnf_p.fitting_options_lb,cnf_p.fitting_options_ub,cnf_p.fitting_options);
                            % JPLZ: modify the input data format, I was obtaining errors
                            [fit_params,~,~,flag]     =   lsqcurvefit (mpfun,fit_params_ini,range_index,fit_data.',cnf_p.fitting_options_lb,cnf_p.fitting_options_ub,cnf_p.fitting_options);
                        case 'fmin'
                            [fit_params,~,flag]     =   fminsearchbnd (fminfun,fit_params_ini,zeros(1,length(fit_params_ini)),cnf_p.fitting_options_lb,cnf_p.fitting_options);
                    end
                otherwise
                    error('No valid mode')
                    
            end
            
            % -------------------------------------------------------------------------
            % DEFINE OUTPUT STRUCTURE
            % -------------------------------------------------------------------------        
            if cnf_p.lut_flag
                switch cnf_p.power_wfm_model
                    case 'simple'
                        ml_wav          =   ml_wave_gen_EM(range_index,fit_params,nf_params, cnf_p,cst, chd,looks,func_f0);
                        
                    case 'complete'
                        ml_wav          =   ml_wave_gen_EM(range_index,fit_params,nf_params, cnf_p,cst, chd,looks,func_f0,func_f1);
                end
            else
                ml_wav          =   ml_wave_gen_EM(range_index,fit_params,nf_params, cnf_p,cst, chd,looks);
            end
            
            %correlation_fit         =   corrcoef(fit_data(range_index),ml_wav);
            % JPLZ: changed the output format, was yielding errors
            correlation_fit         =   corrcoef(fit_data,ml_wav);

            fit_res.COR(m)          =   correlation_fit(1,2);
            
            fit_res.Epoch(m)    =   fit_params(1) + cnf_p.IFmask_N*cnf_p.ZP - 1; % include again the IF mask contribution
            %fit_res.Epoch(m)    =   fit_params(1) + cnf_p.IFmask_N*cnf_p.ZP; % include again the IF mask contribution
            
            % !!!!!!!!!!!!!!!!!! THE -1 VALUE MIGHT CAUSE SOME BIAS
            % CHECK IT !!!
            if cnf_p.rou_flag
                fit_res.rou(m)  =   fit_params(2); 
            else
                fit_res.Hs(m)   =   (fit_params(2)*4); % sigma_z = Hs/4
                if cnf_p.two_step_fitting
                    fit_res.rou(m)  =   NaN;
                end
            end
            if strcmp(cnf_p.multilook_option,'Cris_old')                
                fit_res.Pu(m)       =   10*log10(fit_params(3)*max_power_wav);
            else
                fit_res.Pu(m)       =   10*log10(max_power_wav);
            end            
            fit_res.flag(m)     =   flag;            
            if cnf_p.two_step_fitting
                ml_wav_1st_iter = ml_wav;
                COR_1st_iter=fit_res.COR(m);
                Epoch_1st_iter=fit_res.Epoch(m);
                SWH_1st_iter=fit_res.Hs(m);
                Pu_1st_iter=fit_res.Pu(m);
            end            
            fit_res.flag_validity_L1B(m) =1; %fiting valid 
            
            % -------------------------------------------------------------------------
            % SET NEXT SEEDS / INITIAL VALUES OF FITTING PARAMS
            % -------------------------------------------------------------------------
%             if fit_params(3) <0
%                 %Pu
%                 %force setting this parameter to 0
%                 fit_params(3)=0;
%             end
            % forcing the amplitude fitting to 1
            if cnf_p.initial_param_fit_feedback_flag
                fit_params_ini(3) = 1;
                
                % initial epoch fitting
                % it has much no sense as right now take one provided by the
                % peak retracker
                if cnf_p.seed && m<length(data.lat)-1
                    if (mean(data.MEA.seed)*0.75 <= data.MEA.seed(m+1)) && (data.MEA.seed(m+1)<  mean(data.MEA.seed)*1.25) % rank +-25% over the mean value
                        fit_params_ini(1)  =   data.MEA.seed(m+1);
                    end
                else
                    fit_params_ini(1)  = fit_params(1);
                end
                
                % initial SWH/roughness seed
                fit_params(2)=abs(fit_params(2));
                if m~=1
                    if cnf_p.rou_flag
                        average_total=nanmean(fit_res.rou(1:m-1));
                        std_total=nanstd(fit_res.rou(1:m-1));                        
                    else
                        average_total=nanmean(fit_res.Hs(1:m-1)/4);
                        std_total=nanstd(fit_res.Hs(1:m-1)/4);                                              
                    end
                    
                    if isnan(average_total)
                        average_total=fit_params(2);
                    end
                    if isnan(std_total)
                        std_total=0.0;
                    end
                    
                    if cnf_p.ini_Hs_rou_sliding_win_opt==1
                        %use sliding window to estimate the next initial seed
                        init_sample=max(m-(cnf_p.ini_Hs_rou_sliding_win_size-1),1);
                        if cnf_p.rou_flag
                            average_sliding=nanmean(fit_res.rou(init_sample:m-1));
                            std_sliding=nanstd(fit_res.rou(init_sample:m-1));                                                        
                        else
                            average_sliding=nanmean(fit_res.Hs(init_sample:m-1)/4);
                            std_sliding=nanstd(fit_res.Hs(init_sample:m-1)/4);
                        end
                        
                        if isnan(average_sliding)
                            average_sliding=fit_params(2);
                        end
                        if isnan(std_sliding)
                            std_sliding=0.0;
                        end
                        
                        %to discard effects of having a set of different
                        %consecutive non-valid SWH/rou estimates
                        if (average_sliding < average_total-cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total) || (average_sliding > average_total+cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total)
                            %take the average from the very beginning
                            fit_params_ini(2)=average_total;
                        else
                            current_value=fit_params(2);
                            if ((current_value >= average_sliding-cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding) && (current_value <= average_sliding+cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding))
                                %if within the limits +- std of sliding include
                                %them in average
                                if cnf_p.rou_flag
                                    fit_params_ini(2)=nanmean(fit_res.rou(init_sample:m));
                                else
                                    fit_params_ini(2)=nanmean(fit_res.Hs(init_sample:m))/4.0;
                                end
                            else
                                %if not wihtin limits +- std sliding take the
                                %sliding average
                                fit_params_ini(2)=average_sliding;
                            end
                        end
                    else
                        accumulated=accumulated+fit_params(2);
                        fit_params_ini(2) = (accumulated/m_fitted);
                    end
                else
                    fit_params_ini(2)=fit_params(2);
                end
            end
            
           
           
           
            
            % -------------------------------------------------------------
            % ----        TWO STEP FITTING --------------------------------
            % -------------------------------------------------------------
            if cnf_p.two_step_fitting && fit_res.COR(m)<cnf_p.two_step_fitting_COR_threshold_rou/100
                %run the roughness fitting
                cnf_p.rou_flag=1;
                m_rou=m_rou+1;
                
                if m_rou==1
                    if strcmp(cnf_p.multilook_option,'Cris_old')
                        fit_params_ini_rou      =   [fit_params_ini(1) ini_p.rou fit_params_ini(3)];
                    else
                        fit_params_ini_rou      =   [fit_params_ini(1) ini_p.rou];
                    end
                else
                    if cnf_p.pre_processing
                        fit_params_ini_rou(1)  = idx_leading_edge;
                    else
                        fit_params_ini_rou(1)  = fit_params(1);
                    end
                end
                
                % -------------------------------------------------------------
                % --------- DEFINE FITTING FUNCTION MODEL ---------------------
                % -------------------------------------------------------------
                if cnf_p.lut_flag
                    switch cnf_p.power_wfm_model
                        case 'simple'
                            switch cnf_p.fitting_fun_type
                                case 'lsq'
                                    mpfun               =   @(fit_params,x)ml_wave_gen_EM(x, fit_params, nf_params, cnf_p,cst, chd,looks,func_f0);
                                case 'fmin'
                                    fminfun               =   @(fit_params)sum((ml_wave_gen_EM(range_index, fit_params, nf_params, cnf_p,cst, chd,looks,func_f0)-fit_data).^2);
                            end
                            
                        case 'complete'
                            switch cnf_p.fitting_fun_type
                                case 'lsq'
                                    mpfun               =   @(fit_params,x)ml_wave_gen_EM(x, fit_params, nf_params, cnf_p,cst, chd,looks,func_f0,func_f1);
                                case 'fmin'
                                    fminfun               =   @(fit_params)sum((ml_wave_gen_EM(range_index, fit_params, nf_params, cnf_p,cst, chd,looks,func_f0,func_f1)-fit_data).^2);
                            end
                    end
                else
                    switch cnf_p.fitting_fun_type
                        case 'lsq'
                            mpfun               =   @(fit_params,x)ml_wave_gen_EM(x, fit_params, nf_params, cnf_p,looks);
                        case 'fmin'
                            fminfun               =   @(fit_params)sum((ml_wave_gen_EM(range_index, fit_params, nf_params, cnf_p,looks)-fit_data).^2);
                    end
                end
                
                % -------------------------------------------------------------
                % --------- DEFINE TYPE OF MINIMIZATION PROCEDURE -------------
                % -------------------------------------------------------------
                switch cnf_p.mode
                    case {'SAR','SARin','SARIN','HR','RMC','RAW'}
                        switch cnf_p.fitting_fun_type
                            case 'lsq'
                                [fit_params,~,~,flag]     =   lsqcurvefit (mpfun,fit_params_ini_rou,1:length(fit_data),fit_data.',cnf_p.fitting_options_lb,cnf_p.fitting_options_ub,cnf_p.fitting_options);
                            case 'fmin'
                                [fit_params,~,flag]     =   fminsearchbnd (fminfun,fit_params_ini_rou,zeros(1,length(fit_params_ini_rou)),cnf_p.fitting_options_lb,cnf_p.fitting_options);
                        end
                    otherwise
                        error('No valid mode')
                        
                end
                
                % -------------------------------------------------------------------------
                % DEFINE OUTPUT STRUCTURE
                % -------------------------------------------------------------------------
                
                if cnf_p.lut_flag
                    switch cnf_p.power_wfm_model
                        case 'simple'
                            ml_wav          =   ml_wave_gen_EM(range_index,fit_params,nf_params, cnf_p, cst, chd, looks,func_f0);
                            
                        case 'complete'
                            ml_wav          =   ml_wave_gen_EM(range_index,fit_params,nf_params, cnf_p, cst, chd, looks,func_f0,func_f1);
                    end
                else
                    ml_wav          =   ml_wave_gen_EM(range_index,fit_params,nf_params, cnf_p, cst, chd, looks);
                end
                ml_wav_2nd_iter=ml_wav;
                correlation_fit         =   corrcoef(fit_data,ml_wav);
                
                COR_2nd_iter          =   correlation_fit(1,2);
                Epoch_2nd_iter    =   fit_params(1) + cnf_p.IFmask_N*cnf_p.ZP - 1; % include again the IF mask contribution
                %Epoch_2nd_iter    =   fit_params(1) + cnf_p.IFmask_N*cnf_p.ZP; % include again the IF mask contribution
                rou_2nd_iter  =   fit_params(2);               
                if strcmp(cnf_p.multilook_option,'Cris_old')
                    Pu_2nd_iter       =   10*log10(fit_params(3)*max_power_wav);
                else
                    Pu_2nd_iter       =   10*log10(max_power_wav);
                end
                
                if correlation_fit(1,2)>COR_1st_iter                    
                    fit_res.COR(m)          =   COR_2nd_iter;
                    fit_res.Epoch(m)    =   Epoch_2nd_iter; % include again the IF mask contribution
                    fit_res.rou(m)  =   rou_2nd_iter;
                    fit_res.Hs(m) = NaN;                    
                    fit_res.Pu(m)       =   Pu_2nd_iter;                    
                    fit_res.flag(m)     =   flag;
                else
                    ml_wav=ml_wav_1st_iter;
                end
                
                % -------------------------------------------------------------------------
                % SET NEXT SEED IF NEEDED
                % -------------------------------------------------------------------------
                %             if fit_params(3) <0
                %                 %Pu
                %                 %force setting this parameter to 0
                %                 fit_params(3)=0;
                %             end
                if cnf_p.initial_param_fit_feedback_flag
                    % forcing the amplitude fitting to 1
                    fit_params_ini(3) = 1;
                    
                    % initial epoch fitting
                    % it has much no sense as right now take one provided by the
                    % peak retracker
                    if cnf_p.seed && m<length(data.lat)-1
                        if (mean(data.MEA.seed)*0.75 <= data.MEA.seed(m+1)) && (data.MEA.seed(m+1)<  mean(data.MEA.seed)*1.25) % rank +-25% over the mean value
                            fit_params_ini(1)  =   data.MEA.seed(m+1);
                        end
                    else
                        fit_params_ini_rou(1)  = fit_params(1);
                    end
                    
                    % initial roughness seed
                    fit_params(2)=abs(fit_params(2));
                    if m~=1
                        average_total=nanmean(fit_res.rou(1:m-1));
                        std_total=nanstd(fit_res.rou(1:m-1));
                        
                        if isnan(average_total)
                            average_total=fit_params(2);
                        end
                        if isnan(std_total)
                            std_total=0.0;
                        end
                        
                        if cnf_p.ini_Hs_rou_sliding_win_opt==1
                            %use sliding window to estimate the next initial seed
                            init_sample=max(m-(cnf_p.ini_Hs_rou_sliding_win_size-1),1);
                            
                            average_sliding=nanmean(fit_res.rou(init_sample:m-1));
                            std_sliding=nanstd(fit_res.rou(init_sample:m-1));
                            
                            if isnan(average_sliding)
                                average_sliding=fit_params(2);
                            end
                            if isnan(std_sliding)
                                std_sliding=0.0;
                            end
                            
                            %to discard effects of having a set of different
                            %consecutive non-valid SWH/rou estimates
                            if (average_sliding < average_total-cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total) || (average_sliding > average_total+cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total)
                                %take the average from the very beginning
                                fit_params_ini_rou(2)=average_total;
                            else
                                current_value=fit_params(2);
                                if ((current_value >= average_sliding-cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding) && (current_value <= average_sliding+cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding))
                                    %if within the limits +- std of sliding include
                                    %them in average
                                    fit_params_ini_rou(2)=nanmean(fit_res.rou(init_sample:m));
                                else
                                    %if not wihtin limits +- std sliding take the
                                    %sliding average
                                    fit_params_ini_rou(2)=average_sliding;
                                end
                            end
                        else
                            accumulated_rou=accumulated_rou+fit_params(2);
                            fit_params_ini_rou(2) = (accumulated_rou/m_rou);
                        end
                    else
                        fit_params_ini_rou(2)=fit_params(2);
                    end
                end                

            end

            

% -------------------------------------------------------------------------
% PLOT FITTING RESULTS
% -------------------------------------------------------------------------       
            if cnf_p.plot_fits_flag && ((mod(m,cnf_p.plot_fits_downsampling)==0) || (m==1)) && (data.lat(m)>=cnf_p.plot_fits_lat_range(1) && data.lat(m)<=cnf_p.plot_fits_lat_range(2))
                plot(fit_data,'-b');
                hold on                
                if cnf_p.two_step_fitting %&& COR_1st_iter<cnf_p.two_step_fitting_COR_threshold_rou/100
                    plot(ml_wav_1st_iter,'--r'); %ocean-like
                    plot(ml_wav_2nd_iter,'.-g'); %ice-like 
                    annotation('textbox', [0.725 0.25 0.18 0.35],...
                        'String',{['Epoch (1st)= ', num2str(Epoch_1st_iter,4), ' [r.b]'],...
                        ['Epoch (2nd)= ', num2str(Epoch_2nd_iter,4), ' [r.b]'],...
                        ['SWH = ' ,num2str(abs(SWH_1st_iter),4), ' [m]'],...
                        ['Rou = ' ,num2str(abs(rou_2nd_iter),'%10.5e'), ' [-]'],...
                        ['Pu (1st) = ', num2str(Pu_1st_iter,4), ' [ - ]'],...
                        ['Pu (2nd) = ', num2str(Pu_2nd_iter,4), ' [ - ]'],...
                        ['r (1st) = ', num2str(COR_1st_iter*100,5), '[%]'],...
                        ['r (2nd) = ', num2str(COR_2nd_iter*100,5), '[%]'] },...                        
                        'FontSize',18,...
                        'FontName','Helvetica',...
                        'BackgroundColor',[1 1 1]);
                    h_leg=legend('L1b-Waveform', 'Analytical fit (SWH)','Analytical fit (Rou)');
                    type_fit='2step';
                else
                    plot(ml_wav, '--r')
                    if cnf_p.rou_flag
                        annotation('textbox', [0.7 0.5 0.2 0.2],...
                            'String',{['Epoch = ', num2str(fit_res.Epoch(m),4), ' [r.b]'],['Rou = ' ,...
                            num2str(abs(fit_res.rou(m)),'%10.5e'), ' [m]'],['Pu = ', num2str(fit_res.Pu(m),4), ' [ - ]'],...
                            ['r = ', num2str(fit_res.COR(m)*100,5), '[%]'], ['LSFlag = ', num2str(fit_res.flag(m))],...
                            ['sigma0 = ', num2str(fit_res.Pu(m)+data.sigma0(m)), ' [ dB ]']},...
                            'FontSize',18,...
                            'FontName','Helvetica',...
                            'BackgroundColor',[1 1 1]);
                        type_fit='MSS';
                    else                      
                        annotation('textbox', [0.6 0.6 0.3 0.3],...
                            'String',{['Epoch = ', num2str(fit_res.Epoch(m),4), ' [r.b]'],['SWH = ' ,...
                            num2str(abs(fit_res.Hs(m)),4), ' [m]'],['Pu = ', num2str(fit_res.Pu(m),4), ' [ dB ]'], ...
                            ['r = ', num2str(fit_res.COR(m)*100,5), '[%]'], ['LSFlag = ', num2str(fit_res.flag(m))]...
                             },...
                            'FontSize',18,...
                            'FontName','Helvetica',...
                            'BackgroundColor',[1 1 1]);
                        type_fit='SWH';
                    end
                    h_leg=legend('L1b-Waveform', 'Analytical fit');
                end

                grid on
                xlabel('range bin');
                title(['wav # = ', num2str(m)]);
                axis([1 length(fit_data) 0 1.0]);                
                
                if cnf_p.optional_ext_file_flag
                    print('-dpng ',[path_Results,'plots/fitted_waveforms/','_',cnf_p.file_ext_string,'_L2_fitting_',type_fit,'_record_',num2str(m),'.png']);
                else
                    print('-dpng ',[path_Results,'plots/fitted_waveforms/','L2_fitting_',type_fit,'_record_',num2str(m),'.png']);
                end
                disp(strcat('Waveform: ',num2str(m)));
                save([path_Results,'plots/fitted_waveforms/','L2_fitting_record_',num2str(m),'.mat'],'ml_wav','fit_data','fit_res');
                set(h_leg,'FontName','Helvetica','FontSize',18)
                hold off                
            end % end ploting option
            %disp(strcat('Time wvfm ',num2str(m),' [sec]:',num2str(toc(t_fitting_init))));
            if cnf_p.two_step_fitting
                cnf_p.rou_flag=0;
            end
    end % end for  
    %disp(strcat('Normalized Thermal noise average estimation ',num2str(ThN_estimate_accumulated/m)))
end % end function

%             % initial SWH seed
%             fit_params(2)=abs(fit_params(2));
%             if m~=1
%                 Hs_average_total=mean(fit_res.Hs(1:m-1));
%                 Hs_std_total=std(fit_res.Hs(1:m-1));
%                 if cnf_p.ini_Hs_sliding_win_opt==1
%                     %use sliding window to estimate the next initial seed
%                     init_sample=max(m-(cnf_p.ini_Hs_sliding_win_size-1),1);
%                     Hs_average_sliding=mean(fit_res.Hs(init_sample:m-1));
%                     Hs_std_sliding=std(fit_res.Hs(init_sample:m-1));
%                     %to discard effects of having a set of different
%                     %consecutive non-valid SWH estimates
%                     if (Hs_average_sliding < Hs_average_total-cnf_p.ini_Hs_sliding_win_opt_discard_std_threshold*Hs_std_total) || (Hs_average_sliding > Hs_average_total+cnf_p.ini_Hs_sliding_win_opt_discard_std_threshold*Hs_std_total)
%                         %take the average from the very beginning
%                         fit_params_ini(2)=Hs_average_total/4.0;
%                     else
%                         Hs_current=fit_params(2)*4;
%                         if ((Hs_current >= Hs_average_sliding-cnf_p.ini_Hs_sliding_win_opt_discard_std_threshold*Hs_std_sliding) && (Hs_current <= Hs_average_sliding+cnf_p.ini_Hs_sliding_win_opt_discard_std_threshold*Hs_std_sliding))                                                             
%                             %if within the limits +- std of sliding include
%                             %them in average 
%                             fit_params_ini(2)=mean(fit_res.Hs(init_sample:m))/4.0;
%                         else                 
%                             %if not wihtin limits +- std sliding take the
%                             %sliding average
%                             fit_params_ini(2)=Hs_average_sliding/4.0;
%                         end                                                
%                     end                                        
%                 else                    
%                     accumulated_Hs=accumulated_Hs+fit_params(2)*4;
%                     fit_params_ini(2) = (accumulated_Hs/m)/4.0;
%                 end
%             else
%                 fit_params_ini(2)=fit_params(2);
%             end
