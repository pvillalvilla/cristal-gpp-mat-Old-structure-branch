% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% CR2 
% ---------------------------------------------------------
% Objective: Define Configuration parameters for L2 processing
% 
% INPUTs : - 
% OUTPUTs: -
%
% ----------------------------------------------------------
% Author:    Eduard Makhoul / isardSAT
% Reviewer:  M�nica Roca   / isardSAT
% Last rev.: M�nica Roca   / isardSAT (11/05/2015)
%
%
% Version  record
% 1.0 2016/04/04 Based on Code of Sentinel-6
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------- MISSION & MODE -------------------------------
% -------------------------------------------------------------------------
cnf_L2.mission           =   'CS2';  % this shall be CS2 or JCS;
cnf_L2.mode              =   'SARin';      % this shall be LRM, SAR, SARin --> as in CS2 or LRM(1), HRM_raw(2) and HRM_rmc(3) for Jason-CS/S6
cnf_L2.Nrmc              =   240;    % This will only be used in HRM_rmc processing


%% ------------------------- L1 PROCESSOR ---------------------------------
% -------------------------------------------------------------------------
cnf_L2.L1proc            =   'GPPICE';   
% The L1 processor shall be 
%   ESA --> ESA processing baseline for the different type of missions
%   ISD --> isardSAT processing baseline based on Sentinel-6 GPP
cnf_L2.use_zeros_cnf     =   1;

%% ------------------------ L2 PROCESSOR ----------------------------------
cnf_L2.retracker_name={'THRESHOLD','OCOG','ANALYTICAL','ANALYTICAL'}; %vector to run the different retrackers on the same data: 
                                 %'SAMOSA' or 'ANALYTICAL'--> analytical retracker
                                 %'THRESHOLD' --> Simple percentage of maximum of peak                                
                                 %'OCOG_ice' --> Offset Center of Gravity
%Rough approach when calling the ANALYTICAL more than once to differentiate
%them whether fitting is on SWH, ROUGHNESS or both of them (still to be implemneted)
cnf_L2.analytical_type_of_fitting={'SWH','MSS'}; %common structure 

%-------------------------------------------------------------------------
%----------------------- REF. SAMPLES ------------------------------------
%-------------------------------------------------------------------------
cnf_L2.ref_sample_wd=128/2;

%--------------------------------------------------------------------------
%----------------------- SEED INFORMATION ---------------------------------
%--------------------------------------------------------------------------
cnf_L2.seed              =   0;              % 0 - no seed introduced; 1 - seed provided by phase information of SARin data

%--------------------------------------------------------------------------
% -------------------- FILTER BY GEOGRAPHICAL LOCATION --------------------
%--------------------------------------------------------------------------
cnf_L2.mask_ROI_flag               =   0;              % 0 processe all track; Otherwise specify a KML file

%--------------------------------------------------------------------------
% -------------------- FILTER BY # LOOKS/BEAMS ----------------------------
%--------------------------------------------------------------------------
% define the minimum number of looks available per waveform to consider it
% for the processing --> we exclude those waveforms not generated with 95% of the looks available in the mode theoretically
cnf_L2.mask_looks_flag=0;
if strcmp(cnf_L2.mission,'CS2') || strcmp(cnf_L2.mission,'CR2')
    cnf_L2.Neff_thres        =   256*0.95;    % we tolerate up to 5% of missing echoes max
    if cnf_L2.L1proc == 2
        cnf_L2.Neff_thres    =   212;
    end
else
    cnf_L2.Neff_thres    =   486;
end

%--------------------------------------------------------------------------
% -------------------- IF MASK --------------------------------------------
%--------------------------------------------------------------------------
cnf_L2.IFmask_N          =   0;              % First/last N samples of the waveform affected by the IF Mask filter which should be excluded from processing. Specify if known, 0 if not known

%--------------------------------------------------------------------------
%--------------------- L1B processing characteristics ---------------------
%--------------------------------------------------------------------------
cnf_L2.ZP                =   1;              % 1 if no ZP, 2 if ZP of 2, 4 if ZP of 4, etc...
cnf_L2.window_type_a     =   'Boxcar';       % ALONG TRACK window type previous to AFFT can be Hamming, Hanning or Boxcar (Note: use capital letters)
cnf_L2.window_type_r     =   'Boxcar';       % ACROSS TRACK window type previous to AFFT can be Hamming, Hanning or Boxcar (Note: use capital letters)

%% ----------------- THRESHOLD-BASED RETRACKER ----------------------------
if any(strcmp(cnf_L2.retracker_name,'THRESHOLD'))
 cnf_L2.th_retracker.percentage_peak=50;
end

%% ----------------- OCOG RETRACKER ----------------------------
if any(strcmp(cnf_L2.retracker_name,'OCOG'))
   cnf_L2.OCOG_retracker.percentage_pow_OCOG=30;
   cnf_L2.OCOG_retracker.n1=35; %first zero-padded sample to be used in leading edge estimation
   cnf_L2.OCOG_retracker.n2=246; %last zero-padded sample to be used in leading edge estimation
   cnf_L2.OCOG_retracker.offset=3.5; %additional offset
   cnf_L2.OCOG_retracker.param_comp_method=0; % 0: Using squares of power samples as per Frappart
                                       % 1: Using the power samples as per Wingham
   cnf_L2.OCOG_retracker.implementation_method=0; % 0: Using threshold-based approach computed within n1 and n2
                                   % 1: Computing Amplitude reference 
                                   % within limits defined by the window
                                   % centered at the COG
                                   % 2: Using epoch=offset+(COG-W/2)

end


%% ----------------- ANALYTICAL RETRACKER ---------------------------------
cnf_L2.wvfm_portion_selec=0;
cnf_L2.wvfm_portion_selec_type='peak_valley';
cnf_L2.wvfm_portion_selec_l_samples=50;
cnf_L2.wvfm_portion_selec_r_samples=100;
% cnf_L2.wvfm_discard_samples=0;
cnf_L2.Doppler_mask_cons_option='internal';
cnf_L2.Doppler_mask_cons_internal='l1b_like';
cnf_L2.Thn_estimation_method= 'adaptive';
cnf_L2.Thn_estimation_method= 'fixed_window';
% cnf_L2.Thn_estimation_method= 'external';
% cnf_L2.external_Thn_value=0.2;
cnf_L2.Thn_w_first=1;
cnf_L2.Thn_w_width=10;
cnf_L2.Thn_ML_SL_method= 'ML';
cnf_L2.factor_increase_noise_iter = 0.01;
cnf_L2.quality_flag_misfit_th = 100;
cnf_L2.wvfm_discard_samples = 1;
cnf_L2.wvfm_discard_samples_begin=5;
cnf_L2.wvfm_discard_samples_end=10; 
cnf_L2.wvfm_select_samples_around_peak = 0;
cnf_L2.wvfm_select_samples_leftofpeak=30;
cnf_L2.wvfm_select_samples_rightofpeak=10;  
cnf_L2.wvfm_select_samples_rightofpeak=70;

% SAMOSA based retracker
if any(strcmp(cnf_L2.retracker_name,'ANALYTICAL') | strcmp(cnf_L2.retracker_name,'SAMOSA'))       
    
    %--------------------------------------------------------------------------
    % ----------------------NOISE ESTIMATION ----------------------------------
    %--------------------------------------------------------------------------
    cnf_L2.Thn_flag          =   1;              % 1 - To account for ThN in estimation; 0 - Otherwise
    cnf_L2.Thn_w_first       =   1;%50;             % Gate number after IF Mask to start Thermal noise windowing; this is a subscript indice thus must be > 0
    cnf_L2.Thn_w_width       =   20;%20;             % Thermal noise window width in range bins
    cnf_L2.fit_noise         =   1; % 0 using the values in Thn_w_first and Thn_w_width or 1
    %to perform an adaptive estimation of the window width based
    %on pre-processing step estimation using the
    %multilooked waveform (derivative and leading edge detection)
    %PERFORMED WHENEVER PRE-PROCESSING FLAG IS
    %ACTIVE
    cnf_L2.threshold_noise   = 1e-3; %threshold used to estimate the samples used in the noise window estimation
    cnf_L2.max_iter_noise    = 100;
    cnf_L2.Thn_method        =  'ML'; % String:
    % -'ML': using the multilloked waveform
    %  'SL': estimate the noise per look from
    %  the real stack (require to pass the whole stack)

    %--------------------------------------------------------------------------
    %----------------------- LOOK INDEXATION METHOD ---------------------------
    %--------------------------------------------------------------------------
    cnf_L2.looks_index_method=   'Look_angle';   % String:
    % - Cris_old (old version),
    % - Norm_index (norm. #pulses),
    % - Doppler_freq (exploiting beam angle, velocity and PRI info)
    % - Look_angle angle method (for flat surface: equivalent to pitch + pointing angle)
    switch cnf_L2.looks_index_method
        case 'Look_angle'
            cnf_L2.look_ang_method='approximate'; % String:
            % - approximate: exploiting the min and max
            % values of the Dopp., and Look angles and
            % velocity satellite's over the surface with the given PRI information per stack.
            % - exact : using the beam angles,
            % satellite's velocity, PRI
            % information for each beam in stack
        case 'Doppler_freq'
            cnf_L2.fd_method='approximate'; % String:
            % - approximate: exploiting the min and max
            % values of the Dopp., and Look angles and
            % velocity satellite's over the surface with the given PRI information per stack.
            % - exact : using the beam angles,
            % satellite's velocity, PRI
            % information for each beam in stack
    end
    %--------------------------------------------------------------------------
    % ------------------- MODEL CONFIGURATION ---------------------------------
    %--------------------------------------------------------------------------
    % ------------------ Surface parameters -----------------------------------
    cnf_L2.rou_flag          =   0;              % NO TOQUIS If Gaussian PDF to be considered, then no need to fit roughness which for Ocean is equivalent to Skewness
    %if ~cnf_L2.rou_flag
        cnf_L2.rou           =   1e-2;           % Gaussian PDF show a roughness value in the order of 10-2
        cnf_L2.Hs            = 1e-6; % value of the Hs when fitting the MSS or roughness
    %end
    
    cnf_L2.sign_pitch=-1.0;
    cnf_L2.sign_roll=1.0;
    
    %------------------- Power waveform model ---------------------------------
    cnf_L2.power_wfm_model='complete'; %Define the model approximation of power wfm whether to compute
    % 'simple': Pkl=Bkl*sqrt(gl)*func_f0
    % 'complete':
    % Pkl=Bkl*sqrt(gl)*(func_f0+Tkl*gl*simga_s^2*func_f1);
    %-------------------- Multilooking option ---------------------------------
    cnf_L2.multilook_option='Cris_old'; %String:
    % -Cris_old: thermal noise is added to
    % each SL
    % -NormML_plus_noise: ML is normalized
    % and the thermal noise is added
    % outside
    % -------------------------------------------------------------------------
    % ---------------------- LOOK UP TABLES CONFIGURATION ---------------------
    % -------------------------------------------------------------------------
    % LUT tables in .mat should be in the processor folder
    cnf_L2.lut_flag          =   1;              % 0 or 1;if null do not use Look up tables
    cnf_L2.N_int             =   200; %not used            % number of intervals in the integration
    cnf_L2.lim_int           =   5.1; %not used           % the range of the integration
    cnf_L2.LUT_ximax         =   50; % maximum value of the gl*k parameter LUT
    cnf_L2.LUT_ximin         =   -10; % minimum value of the gl*k parameter LUT
    cnf_L2.LUT_step          =   1e-5; % step of gl*k parameter in the LUT
    
    %


    %--------------------------------------------------------------------------
    %----------------- PRE-PROCESSING STAGE------------------------------------
    %--------------------------------------------------------------------------
    %provide better initial estimate of epoch--
    %include a pre_processing stage to compute the initial epoch if desired
    %based on a mid-leading edge estimation as a percentatge of the peak power
    %in the actual waveform and eventually include the estimation of the noise
    %floor adjusting automatically
    cnf_L2.pre_processing=1;
    cnf_L2.percent_leading_edge=87.0; % percentage of peak detect to establish the mid-point leading edge
    
    % -------------------------------------------------------------------------
    % ------------- FITTING CONFIGURATION -------------------------------------
    %---- Initial values-------------------------------------------------------
    cnf_L2.ini_Epoch         =   35;             % 45 for BB, 145 for MEAS18
    cnf_L2.ini_Hs            =   2.0;
    cnf_L2.ini_Pu            =   1;
    cnf_L2.init_rou          =   1e-2;
    cnf_L2.initial_param_fit_feedback_flag = 0; %activate the updating of the initial parameters per fitted waveform either using the sliding window or with acumulation of previous info
    cnf_L2.ini_Hs_rou_sliding_win_opt = 1;  % activate the averaging of Hs/roughness for next initial seed 
    cnf_L2.ini_Hs_rou_sliding_win_opt_discard_std_threshold = 3; %number of times std for removal (it is also used for comparison with the accumulated one)
    cnf_L2.ini_Hs_rou_sliding_win_size = 50; 
    
    %------------------ Two-step fitting ----------------------------------
    cnf_L2.two_step_fitting=0; %used to perform first a fitting on SWH and if 
    cnf_L2.two_step_fitting_COR_threshold_rou=99.0; %threshold below which second fitting on roughness is used
    
    %------------------ Range indexation ----------------------------------
    cnf_L2.range_index_method='resampling'; % 'conventional': when creating the range vector used in fitting use directly 1:1:N_samples (even with ZP)
                                             % 'resampling': creating the
                                             % range vector as a resampling
                                             % of the original one (if there is a resampling)
    
    %------------------ Fitting function --------------------------------------
    cnf_L2.fitting_fun_type='lsq'; % Define the fitting procedure used
    
    %String:
    %-'lsq': Using the lsqcurvfit based on levegender-Marquard
    %-'fmin': Using the fminsearch
    %algorithm
    switch cnf_L2.fitting_fun_type
        case 'lsq'
            %trust-region-reflective
            %levenberg-marquardt
            %'MaxIter',400
            cnf_L2.lsq_algorithm      = 'trust-region-reflective';
            %cnf_L2.fitting_options    = optimset('Algorithm',cnf_L2.lsq_algorithm,'FinDiffType','central','TolFun',1e-8,'TolX',1e-8,'Display','off','UseParallel','always');
            cnf_L2.fitting_options    = optimset('Algorithm',cnf_L2.lsq_algorithm,'Display','off','UseParallel',0);
            cnf_L2.fitting_options_lb = [0,0,0]; %lower boundaries
            cnf_L2.fitting_options_ub = []; %upper boundaries
        case 'fmin'
            cnf_L2.fitting_options     =   optimset('Display','off','TolFun',1e-7,'TolX',1e-7);
            cnf_L2.fitting_options_lb = []; %lower boundaries
            cnf_L2.fitting_options_ub = []; %upper boundaries
    end
    
    %--------------------------------------------------------------------------
    %------------------ VERBOSE PARAMETERS-------------------------------------
    %--------------------------------------------------------------------------
    cnf_L2.plot_fits_flag    =   1;              % 0 - not plot; 1 - plot fit results
    cnf_L2.plot_fits_lat_range = [-91,91]; %latitude range of waveforms fit to be plotted (by default between -91 and 91 to force all of them to be plotted)
    cnf_L2.plot_fits_downsampling =10;
    
end

%% -------------------- GEOPHYSICAL CORRECTIONS ----------------------------
% -------------------------------------------------------------------------
cnf_L2.geo_corr_application_flag=1; 
cnf_L2.force_geocorr_surf_type=1; % 0: use the geophysical corrections depending on the flag surface type included in the L1B product                          
                                 % 1: use the same group of geophysical corrections to all the records independently of the flag of the surface type: forcing the type surface
cnf_L2.product_type_surface='open_ocean'; % open_ocean, sea_ice, land_ice ...                                 


%% -------------------- OUTPUT PRODUCT -------------------------------------
% -------------------------------------------------------------------------
%Inherete from Sentinel-3 file formating
cnf_L2.Product_type='LAN'; % String: 'LAN' for land products
                          %'WAT' for water products
cnf_L2.write_output=1;
cnf_L2.output_product_format='nc'; % nc: NetCDF, mat: Matlab
cnf_L2.nc_name_surface_height='shh'; %string indicating which name to be used in the netcdf for the surface height: SSH (sea surface height) or SH (surface height)

cnf_L2.optional_ext_file_flag =0;
%cnf_L2.file_ext_string='PTR_sigmaa_new_sigmar_cris_Lz_Ly_ZP_rg_idx_conv';
cnf_L2.file_ext_string='L2_looks_by_Doppler';


%% ---------------- VERBOSE OPTIONS ---------------------------------------
cnf_L2.visible_figures=0;


switch cnf_L2.window_type_a
    case 'Boxcar'  % the Gaussian approximation for a boxcar window filtered waveform previous to FFT
        %             alpha     =   1;
        %             [A_s2Ga, alpha_ga]    =   sinc2Gauss(alpha, [-2+1/64:4/64:2]);    % note that alpha_g varies with the number of samples not the rank of values
        %             alpha_ga=1/(2*0.513^2);
        chd.A_s2Ga=1.0196;
        chd.alpha_ga=1.0/(2*(0.36012).^2);
    case 'Hanning'  % the Gaussian approximation for a boxcar window filtered waveform previous to FFT
        %             alpha     =   0.5;
        %             [A_s2Ga, alpha_ga]    =   sinc2Gauss(alpha, [-2+1/64:4/64:2]);    % note that alpha_g varies with the number of samples not the rank of values
        chd.A_s2Ga=1.0101;
        chd.alpha_ga=1.0/(2*(0.59824).^2);
    case 'Hamming' % the Gaussian approximation for a hamming window filtered waveform previous to FFT
        %             alpha     =   0.54;
        %             [A_s2Ga, alpha_ga]    =   sinc2Gauss(alpha, [-2+1/64:4/64:2]);    % note that alpha_g varies with the number of samples not the rank of values
        chd.A_s2Ga=1.0081;
        chd.alpha_ga=1.0/(2.0*(0.54351).^2);
    case 'forced'
        chd.alpha_ga=1.0/(2.0*(0.545333774157098).^2);
    otherwise
        error('invalid window type. Chose Hamming, Hanning or Boxcar')
end
disp(strcat('Azimuth PTR approx:',{''},num2str(sqrt(1.0/(2.0*chd.alpha_ga)))))

switch cnf_L2.window_type_r
    case 'Boxcar'  % the Gaussian approximation for a boxcar window filtered waveform previous to FFT
        %             alpha     =   1;
        %             [A_s2Gr, alpha_gr]    =   sinc2Gauss(alpha, [-2+1/64:4/64:2]);    % note that alpha_g varies with the number of samples not the rank of values
        %             alpha_gr =1/(2*0.513^2);
        chd.A_s2Gr=1.0196;
        chd.alpha_gr=1.0/(2*(0.36012).^2);
    case 'Hanning'  % the Gaussian approximation for a boxcar window filtered waveform previous to FFT
        chd.alpha     =   0.5;
        [chd.A_s2Gr, chd.alpha_gr]    =   sinc2Gauss(chd.alpha, [-2+1/64:4/64:2]);    % note that alpha_g varies with the number of samples not the rank of values
        A_s2Gr=1.0101;
        chd.alpha_gr=1.0/(2*(0.59824).^2);
    case 'Hamming' % the Gaussian approximation for a hamming window filtered waveform previous to FFT
        chd.alpha     =   0.54;
        [chd.A_s2Gr, chd.alpha_gr]    =   sinc2Gauss(chd.alpha, [-2+1/64:4/64:2]);    % note that alpha_g varies with the number of samples not the rank of values
        chd.A_s2Gr=1.0081;
        chd.alpha_gr=1.0/(2.0*(0.54351).^2);
    case 'forced'
        chd.alpha_gr=1.0/(2.0*(0.545333774157098).^2);
    otherwise
        error('invalid window type. Chose Hamming, Hanning or Boxcar');
end

%flag to proceeed to iceberg detection method
cnf_L2.iceberg_detection= 0;