function fit_res = fitting(data, cnf_p, nf_p, ini_p,varargin)

% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% Routine to perform the fitting of the waveforms
% -------------------------------------------------------------------------
% 
% Author:           Cristina Martin-Puig / isardSAT
%                   Eduard Makhoul / isardSAT
%
% Reviewer:         Cristina Martin-Puig / isardSAT
%
% Last revision:    Cristina Martin-Puig / isardSAT V7 21/06/2016
%
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% INPUT:
%       data        =   SAR data structure 
%       nf_p        =   Non-fit parameters of the fitting
%       cnf_p       =   Configuration parameters of the L2 processor
%       ini_p       =   initialization parameters of need for the fitting
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
% - ml_wave_gen: genrates the multilook waveform based on Chris Ray
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
%
%
%



%% ---------------- Handling input variables ------------------------------
if(nargin<3 || nargin>6)
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
addParameter(p,'LUT_f0_file',{''},@(x)ischar(x));
addParameter(p,'LUT_f1_file',{''},@(x)ischar(x));
p.parse(varargin{:});
LUT_f0_file=char(p.Results.LUT_f0_file);
LUT_f1_file=char(p.Results.LUT_f1_file);
clear p;

%--------------------------------------------------------------------------
% -------------- LOADING THE LUTS -----------------------------------------
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
    nf_params.Npulses   =   nf_p.Npulses;
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

    if strcmp(cnf_p.mode, 'SAR') || cnf_p.mode == 2 %RAW
        data.HRM.power_fit    =   data.HRM.power_wav(cnf_p.IFmask_N*cnf_p.ZP+1:end-(cnf_p.IFmask_N*cnf_p.ZP),:)';
    elseif  cnf_p.mode == 3 %RMC
        data.HRM.power_fit    =   data.HRM.power_wav(cnf_p.IFmask_N*cnf_p.ZP+1:cnf_p.Nrmc,:)';  
    elseif strcmp(cnf_p.mode, 'SARin')
        data.SIN.power_fit    =   data.SIN.power_wav(cnf_p.IFmask_N*cnf_p.ZP+1:end-(cnf_p.IFmask_N*cnf_p.ZP),:)';
    else
        error('not a valid mode');
    end

% -------------------------------------------------------------------------
% Least Square Fitting
% -------------------------------------------------------------------------
    for m = 1:length(data.GEO.LAT)

            fprintf('-------\n');
            fprintf('- Fitting wave num = %d\n',m-1);

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
            if strcmp(cnf_p.mode, 'SAR') || cnf_p.mode == 2 || cnf_p.mode == 3
                fit_data            =   data.HRM.power_fit(m,:)/max(data.HRM.power_fit(m,:));          
            else
                fit_data            =   data.SIN.power_fit(m,:)/max(data.SIN.power_fit(m,:)); 
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
                        looks=2.0/nf_p.lambda.*data.GEO.V_beam_surf(m,1:nf_p.Neff(m)).*cos(data.HRM.beam_ang_surf(m,1:nf_p.Neff(m))).*...
                            data.HRM.pri_sar_surf(m,1:nf_p.Neff(m)).*nf_p.Npulses;
                    case 'approximate'
                        %compute the max and min beam angles from 
                        % Look and Doppler Angles information
                        beam_angles=pi/2+...
                         linspace(data.HRM.doppler_angle_start_ku(m),data.HRM.doppler_angle_stop_ku(m),nf_p.Neff(m))...  
                        -linspace(data.HRM.look_angle_start_ku(m),data.HRM.look_angle_stop_ku(m),nf_p.Neff(m));
                        looks=2.0/nf_p.lambda.*data.GEO.V(m).*cos(beam_angles).*data.HRM.pri_sar_isp_surf(m).*nf_p.Npulses;
                end 
           case 'Look_angle'
               delta_look_angle=asin(nf_p.lambda./(data.HRM.pri_sar_isp_surf(m).*2.0*...
                                    nf_p.Npulses*data.GEO.V(m))); 
               switch cnf_p.look_ang_method
                    case 'exact'
                        looks=data.HRM.look_ang_surf(m,1:nf_p.Neff(m))/delta_look_angle;
                    case 'approximate'
                        %compute the max and min beam angles from 
                        % Look and Doppler Angles information
                        looks=linspace(data.HRM.look_angle_start_ku(m),...
                            data.HRM.look_angle_stop_ku(m),nf_p.Neff(m))/delta_look_angle;
                end 
               
       end 

       %-----------------------------------------------------------
       %---------- Create the mask (as a matrix)-------------------
       %-----------------------------------------------------------
       %generated once avoiding to create it within the ML waveform
       %function for computational load issues
       if isfield(data.HRM,'Doppler_mask')
           nf_params.Doppler_mask=ones(nf_params.Neff,data.N_samples);
           %to be in-line with L1B-processing
           if cnf_p.use_zeros_cnf == 1
               value_mask=1*0;
           else
               value_mask=NaN;
           end
           for i_look=1:nf_params.Neff
               nf_params.Doppler_mask(i_look,data.HRM.Doppler_mask(m,i_look):end)=value_mask;
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
if cnf_p.pre_processing && m==1
    %used to estimate the initial value of the epoch 
    %& eventually if desired the noise floor
    [peak_pow,idx_max_peak]=max(fit_data);
    idx_leading_edge=find(find(fit_data<=cnf_p.percent_leading_edge/100.0*peak_pow)<idx_max_peak, 1, 'last' );
    fit_params_ini(1)=idx_leading_edge;
    
    %noise floor estimation
    if cnf_p.fit_noise
        switch cnf_p.Thn_method
            case 'ML'
                idx_noise=[];
                temp_noise_thr=cnf_p.threshold_noise;
                while isempty(idx_noise)
                    idx_noise=find(abs(diff(fit_data))<=temp_noise_thr);
                    idx_noise=idx_noise(idx_noise<idx_leading_edge);
                    temp_noise_thr=temp_noise_thr/2;                    
                end
                nf_params.ThN        =   mean(fit_data(idx_noise))*ones(1,length(looks));
                clear idx_noise tempo_noise_thr;
        end
    end
    
    clear idx_leading_edge peak_pow idx_max_peak;
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
            if cnf_p.lut_flag
                switch cnf_p.power_wfm_model
                    case 'simple'
                        switch cnf_p.fitting_fun_type
                                case 'lsq'
                                    mpfun               =   @(fit_params,x)ml_wave_gen_EM(x, fit_params, nf_params, cnf_p,looks,func_f0);
                                case 'fmin'
                                    fminfun               =   @(fit_params)sum((ml_wave_gen_EM(1:length(fit_data), fit_params, nf_params, cnf_p,looks,func_f0)-fit_data).^2);
                        end
                               
                    case 'complete'
                        switch cnf_p.fitting_fun_type
                            case 'lsq'
                                mpfun               =   @(fit_params,x)ml_wave_gen_EM(x, fit_params, nf_params, cnf_p,looks,func_f0,func_f1);
                            case 'fmin'
                                fminfun               =   @(fit_params)sum((ml_wave_gen_EM(1:length(fit_data), fit_params, nf_params, cnf_p,looks,func_f0,func_f1)-fit_data).^2);
                        end                                                       
                end                
            else
                switch cnf_p.fitting_fun_type
                    case 'lsq'
                        mpfun               =   @(fit_params,x)ml_wave_gen_EM(x, fit_params, nf_params, cnf_p,looks);
                    case 'fmin'
                        fminfun               =   @(fit_params)sum((ml_wave_gen_EM(1:length(fit_data), fit_params, nf_params, cnf_p,looks)-fit_data).^2);
                end                                
            end
            
            %tic            
            if strcmp(cnf_p.mode, 'SAR') || cnf_p.mode == 2 || cnf_p.mode == 3
                switch cnf_p.fitting_fun_type
                    case 'lsq'
                        [fit_params,~,~,flag]     =   lsqcurvefit (mpfun,fit_params_ini,1:length(fit_data),fit_data,[],[],cnf_p.fitting_options);
                    case 'fmin'
                        [fit_params,~,flag]     =   fminsearchbnd (fminfun,fit_params_ini,zeros(1,length(fit_params_ini)),[],cnf_p.fitting_options);
                end
            elseif strcmp(cnf_p.mode, 'SARin')
                switch cnf_p.fitting_fun_type
                    case 'lsq'
                        [fit_params,~,~,flag]     =   lsqcurvefit (mpfun,fit_params_ini,1:length(fit_data),fit_data,[],[],cnf_p.fitting_options);
                    case 'fmin'
                        [fit_params,~,flag]     =   fminsearchbnd (fminfun,fit_params_ini,zeros(1,length(fit_params_ini)),[],cnf_p.fitting_options);
                end
            else
                error('not a valid mode');
            end
            %toc
% -------------------------------------------------------------------------
% SET NEXT SEED IF NEEDED
% -------------------------------------------------------------------------        
            if cnf_p.seed && m<length(data.GEO.LAT)-1
                if (mean(data.MEA.seed)*0.75 <= data.MEA.seed(m+1)) && (data.MEA.seed(m+1)<  mean(data.MEA.seed)*1.25) % rank +-25% over the mean value
                    fit_params_ini(1)  =   data.MEA.seed(m+1);
                end
            else
                fit_params_ini      =   fit_params;
            end

% -------------------------------------------------------------------------
% DEFINE OUTPUT STRUCTURE
% -------------------------------------------------------------------------        
            fit_res.Epoch(m)    =   fit_params(1) + cnf_p.IFmask_N*cnf_p.ZP - 1; % include again the IF mask contribution
            if cnf_p.rou_flag
                fit_res.rou(m)  =   fit_params(2); 
            else
                fit_res.Hs(m)   =   fit_params(2)*4; % sigma_z = Hs/4
            end
            if strcmp(cnf_p.multilook_option,'Cris_old')
                fit_res.Pu(m)       =   10*log10(fit_params(3)*max(data.HRM.power_fit(m,:)));
            else
                fit_res.Pu(m)       =   10*log10(max(data.HRM.power_fit(m,:)));
            end
            
            fit_res.flag(m)     =   flag;
            
            if cnf_p.lut_flag
                switch cnf_p.power_wfm_model
                    case 'simple'
                        ml_wav          =   ml_wave_gen_EM(1:length(fit_data),fit_params,nf_params, cnf_p,looks,func_f0);
                        
                    case 'complete'
                        ml_wav          =   ml_wave_gen_EM(1:length(fit_data),fit_params,nf_params, cnf_p,looks,func_f0,func_f1);                        
                end                
            else
                ml_wav          =   ml_wave_gen_EM(1:length(fit_data),fit_params,nf_params, cnf_p,looks);
            end
            
            correlation_fit         =   corrcoef(fit_data,ml_wav);
            fit_res.COR(m)          =   correlation_fit(1,2); 
            
            profile viewer;
            
            if cnf_p.L1proc == 2
                correlation_fit_cnes    =   corrcoef(fit_data,data.RES.model(:,1)');
                fit_res.COR_cnes(m)     =   correlation_fit_cnes(1,2);
            end

% -------------------------------------------------------------------------
% PLOT FITTING RESULTS
% -------------------------------------------------------------------------       

            if cnf_p.plot_fits_flag
                plot(fit_data,'b.-','MarkerSize',8);
                hold on
                plot(ml_wav, 'r.-', 'MarkerSize',8)
                grid on
                xlabel('range bin','FontName','Helvetica','FontSize',11,'FontWeight', 'bold' );
                title(['wav # = ', num2str(m)],'FontName','Helvetica','FontSize',11,'FontWeight', 'bold' );

                if cnf_p.rou_flag
                    annotation('textbox', [0.7 0.5 0.2 0.2],...
                        'String',{['Epoch = ', num2str(fit_res.Epoch(m),4), ' [r.b]'],['Hs = ' ,...
                        num2str(abs(fit_res.Hs(m)),4), ' [m]'],['Pu = ', num2str(fit_res.Pu(m),4), ' [ - ]'],...
                        ['r = ', num2str(fit_res.COR(m)*100,5), '[%]'], ['LSFlag = ', num2str(fit_res.flag(m))] },...
                        'FontSize',18,...
                        'FontName','Helvetica',...
                        'BackgroundColor',[1 1 1]);
                else

                    annotation('textbox', [0.6 0.6 0.3 0.3],...
                        'String',{['Epoch = ', num2str(fit_res.Epoch(m),4), ' [r.b]'],['Hs = ' ,...
                        num2str(abs(fit_res.Hs(m)),4), ' [m]'],['Pu = ', num2str(fit_res.Pu(m),4), ' [ dB ]'], ...
                        ['r = ', num2str(fit_res.COR(m)*100,5), '[%]'], ['LSFlag = ', num2str(fit_res.flag(m))], ['sigma0 = ', num2str(fit_res.Pu(m)+data.HRM.s0_sf(m)), ' [ dB ]'] },...
                        'FontSize',18,...
                        'FontName','Helvetica',...
                        'BackgroundColor',[1 1 1]);

%                        annotation('textbox', [0.7 0.7 0.2 0.2],...
%                         'String',{['Epoch = ', num2str(fit_res.Epoch(m),4), ' [r.b]'],['Hs = ' ,...
%                         num2str(abs(fit_res.Hs(m)),4), ' [m]'],['Pu = ', num2str(fit_res.Pu(m),4), ' [ - ]'], ...
%                         ['r = ', num2str(fit_res.COR(m)*100,5), '[%]'], ['LSFlag = ', num2str(fit_res.flag(m))] ,...
%                         ['BIAS_Hs = ', num2str(abs(fit_res.Hs(m))-data.RES.SWH(m))]},'FontSize',11,...
%                         'FontName','Helvetica',...
%                         'BackgroundColor',[1 1 1]);           
                end

                if cnf_p.L1proc == 0 
                    h_leg=legend('esa L1B-BB wfm', 'isardSAT fit')
                    axis([1 length(fit_data) 0 250])
                    print('-depsc ',['HRM_figures/fit_L1b_esa_BB_',num2str(m),'.jpg']);
                elseif cnf_p.L1proc == 2 % CNES
                    plot(data.RES.model(:,m),'b.-','MarkerSize',8);
                    h_leg=legend('cnes L1B wfm', 'isardSAT fit', 'cnes fit')
                    axis([1 length(fit_data) 0 250])
                    print('-depsc ',['HRM_figures/fit_L1b_cnes_cpp_',num2str(m),'.jpg']);
                elseif cnf_p.L1proc == 1 || cnf_p.L1proc == 3
                     h_leg=legend('L1 S6 CPP', 'isardSAT FIT','Location', 'NorthWest')
                     if cnf_p.mode == 2
                         %Mask applied
                         axis([0 512-2*cnf_p.IFmask_N 0 1]) % before 1 was 180
                         %Mask not applied
                         axis([0 450 0 1]) % before 1 was 180
                     else
                         axis([0 cnf_p.Nrmc 0 1]) % before 1 was 180
                     end
                     print('-depsc ',['resultats/plots/fit_S6GPP_wfm_',num2str(m),'.jpg']);

                end
                set(h_leg,'FontName','Helvetica','FontSize',18)
                hold off
            end


    end % end for
    %save(strcat('resultats/error_fd_exact_approx_',strrep(cnf_p.FileName,'L1B','L2'),'.mat'),'error_fd_exact_approx');
end % end function

