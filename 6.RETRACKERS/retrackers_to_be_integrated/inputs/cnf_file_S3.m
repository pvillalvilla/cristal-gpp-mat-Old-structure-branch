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

global cnf_p
%% ------------------------- MISSION & MODE -------------------------------
% -------------------------------------------------------------------------
cnf_p.mission           =   'CS2';  % this shall be CS2 or JCS; to be updated when S3 is available
cnf_p.mode              =   'SAR';      % this shall be LRM, SAR, SARin --> as in CS2 or LRM(1), HRM_raw(2) and HRM_rmc(3) for Jason-CS/S6
cnf_p.Nrmc              =   240;    % This will only be used in HRM_rmc processing


%% ------------------------- L1 PROCESSOR ---------------------------------
% -------------------------------------------------------------------------
cnf_p.L1proc            =   'ISD';   
% The L1 processor shall be 
%   ESA --> ESA processing baseline for the different type of missions
%   ISD --> isardSAT processing baseline based on Sentinel-6 GPP
cnf_p.use_zeros_cnf     =   1;

%% ------------------------ L2 PROCESSOR ----------------------------------
cnf_p.retracker_name={'THRESHOLD','OCOG','ANALYTICAL','ANALYTICAL'}; %vector to run the different retrackers on the same data: 
                                 %'SAMOSA' or 'ANALYTICAL'--> analytical retracker
                                 %'THRESHOLD' --> Simple percentage of maximum of peak                                
                                 %'OCOG_ice' --> Offset Center of Gravity
%Rough approach when calling the ANALYTICAL more than once to differentiate
%them whether fitting is on SWH, ROUGHNESS or both of them (still to be implemneted)
cnf_p.analytical_type_of_fitting={'SWH','MSS'}; %common structure 

%-------------------------------------------------------------------------
%----------------------- REF. SAMPLES ------------------------------------
%-------------------------------------------------------------------------
cnf_p.ref_sample_wd=128/2;

%--------------------------------------------------------------------------
%----------------------- SEED INFORMATION ---------------------------------
%--------------------------------------------------------------------------
cnf_p.seed              =   0;              % 0 - no seed introduced; 1 - seed provided by phase information of SARin data

%--------------------------------------------------------------------------
% -------------------- FILTER BY GEOGRAPHICAL LOCATION --------------------
%--------------------------------------------------------------------------
cnf_p.mask_ROI_flag               =   0;              % 0 processe all track; Otherwise specify a KML file

%--------------------------------------------------------------------------
% -------------------- FILTER BY # LOOKS/BEAMS ----------------------------
%--------------------------------------------------------------------------
% define the minimum number of looks available per waveform to consider it
% for the processing --> we exclude those waveforms not generated with 95% of the looks available in the mode theoretically
cnf_p.mask_looks_flag=0;
if strcmp(cnf_p.mission,'CS2') || strcmp(cnf_p.mission,'CR2')
    cnf_p.Neff_thres        =   256*0.95;    % we tolerate up to 5% of missing echoes max
    if cnf_p.L1proc == 2
        cnf_p.Neff_thres    =   212;
    end
else
    cnf_p.Neff_thres    =   486;
end

%--------------------------------------------------------------------------
% -------------------- IF MASK --------------------------------------------
%--------------------------------------------------------------------------
cnf_p.IFmask_N          =   0;              % First/last N samples of the waveform affected by the IF Mask filter which should be excluded from processing. Specify if known, 0 if not known

%--------------------------------------------------------------------------
%--------------------- L1B processing characteristics ---------------------
%--------------------------------------------------------------------------
cnf_p.ZP                =   1;              % 1 if no ZP, 2 if ZP of 2, 4 if ZP of 4, etc...
cnf_p.window_type_a     =   'Boxcar';       % ALONG TRACK window type previous to AFFT can be Hamming, Hanning or Boxcar (Note: use capital letters)
cnf_p.window_type_r     =   'Boxcar';       % ACROSS TRACK window type previous to AFFT can be Hamming, Hanning or Boxcar (Note: use capital letters)

%% ----------------- THRESHOLD-BASED RETRACKER ----------------------------
if any(strcmp(cnf_p.retracker_name,'THRESHOLD'))
 cnf_p.th_retracker.percentage_peak=50;
end

%% ----------------- OCOG RETRACKER ----------------------------
if any(strcmp(cnf_p.retracker_name,'OCOG'))
   cnf_p.OCOG_retracker.percentage_pow_OCOG=30;
   cnf_p.OCOG_retracker.n1=35; %first zero-padded sample to be used in leading edge estimation
   cnf_p.OCOG_retracker.n2=128; %last zero-padded sample to be used in leading edge estimation
   cnf_p.OCOG_retracker.offset=3.5; %additional offset
   cnf_p.OCOG_retracker.param_comp_method=0; % 0: Using squares of power samples as per Frappart
                                       % 1: Using the power samples as per Wingham
   cnf_p.OCOG_retracker.implementation_method=0; % 0: Using threshold-based approach computed within n1 and n2
                                   % 1: Computing Amplitude reference 
                                   % within limits defined by the window
                                   % centered at the COG
                                   % 2: Using epoch=offset+(COG-W/2)

end


%% ----------------- ANALYTICAL RETRACKER ---------------------------------
% SAMOSA based retracker
if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))       
    
    %--------------------------------------------------------------------------
    % ----------------------NOISE ESTIMATION ----------------------------------
    %--------------------------------------------------------------------------
    cnf_p.Thn_flag          =   1;              % 1 - To account for ThN in estimation; 0 - Otherwise
    cnf_p.Thn_w_first       =   1;%50;             % Gate number after IF Mask to start Thermal noise windowing; this is a subscript indice thus must be > 0
    cnf_p.Thn_w_width       =   20;%20;             % Thermal noise window width in range bins
    cnf_p.fit_noise         =   1; % 0 using the values in Thn_w_first and Thn_w_width or 1
    %to perform an adaptive estimation of the window width based
    %on pre-processing step estimation using the
    %multilooked waveform (derivative and leading edge detection)
    %PERFORMED WHENEVER PRE-PROCESSING FLAG IS
    %ACTIVE
    cnf_p.threshold_noise   = 1e-3; %threshold used to estimate the samples used in the noise window estimation
    cnf_p.max_iter_noise    = 100;
    cnf_p.Thn_method        =  'ML'; % String:
    % -'ML': using the multilloked waveform
    %  'SL': estimate the noise per look from
    %  the real stack (require to pass the whole stack)

    %--------------------------------------------------------------------------
    %----------------------- LOOK INDEXATION METHOD ---------------------------
    %--------------------------------------------------------------------------
    cnf_p.looks_index_method=   'Look_angle';   % String:
    % - Cris_old (old version),
    % - Norm_index (norm. #pulses),
    % - Doppler_freq (exploiting beam angle, velocity and PRI info)
    % - Look_angle angle method (for flat surface: equivalent to pitch + pointing angle)
    switch cnf_p.looks_index_method
        case 'Look_angle'
            cnf_p.look_ang_method='approximate'; % String:
            % - approximate: exploiting the min and max
            % values of the Dopp., and Look angles and
            % velocity satellite's over the surface with the given PRI information per stack.
            % - exact : using the beam angles,
            % satellite's velocity, PRI
            % information for each beam in stack
        case 'Doppler_freq'
            cnf_p.fd_method='approximate'; % String:
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
    cnf_p.rou_flag          =   0;              % NO TOQUIS If Gaussian PDF to be considered, then no need to fit roughness which for Ocean is equivalent to Skewness
    %if ~cnf_p.rou_flag
        cnf_p.rou           =   1e-2;           % Gaussian PDF show a roughness value in the order of 10-2
        cnf_p.Hs            = 1e-6; % value of the Hs when fitting the MSS or roughness
    %end
    
    cnf_p.sign_pitch=-1.0;
    cnf_p.sign_roll=1.0;
    
    %------------------- Power waveform model ---------------------------------
    cnf_p.power_wfm_model='complete'; %Define the model approximation of power wfm whether to compute
    % 'simple': Pkl=Bkl*sqrt(gl)*func_f0
    % 'complete':
    % Pkl=Bkl*sqrt(gl)*(func_f0+Tkl*gl*simga_s^2*func_f1);
    %-------------------- Multilooking option ---------------------------------
    cnf_p.multilook_option='Cris_old'; %String:
    % -Cris_old: thermal noise is added to
    % each SL
    % -NormML_plus_noise: ML is normalized
    % and the thermal noise is added
    % outside
    % -------------------------------------------------------------------------
    % ---------------------- LOOK UP TABLES CONFIGURATION ---------------------
    % -------------------------------------------------------------------------
    % LUT tables in .mat should be in the processor folder
    cnf_p.lut_flag          =   1;              % 0 or 1;if null do not use Look up tables
    cnf_p.N_int             =   200; %not used            % number of intervals in the integration
    cnf_p.lim_int           =   5.1; %not used           % the range of the integration
    cnf_p.LUT_ximax         =   50; % maximum value of the gl*k parameter LUT
    cnf_p.LUT_ximin         =   -10; % minimum value of the gl*k parameter LUT
    cnf_p.LUT_step          =   1e-5; % step of gl*k parameter in the LUT
    
    %


    %--------------------------------------------------------------------------
    %----------------- PRE-PROCESSING STAGE------------------------------------
    %--------------------------------------------------------------------------
    %provide better initial estimate of epoch--
    %include a pre_processing stage to compute the initial epoch if desired
    %based on a mid-leading edge estimation as a percentatge of the peak power
    %in the actual waveform and eventually include the estimation of the noise
    %floor adjusting automatically
    cnf_p.pre_processing=1;
    cnf_p.percent_leading_edge=87.0; % percentage of peak detect to establish the mid-point leading edge
    
    % -------------------------------------------------------------------------
    % ------------- FITTING CONFIGURATION -------------------------------------
    %---- Initial values-------------------------------------------------------
    cnf_p.ini_Epoch         =   35;             % 45 for BB, 145 for MEAS18
    cnf_p.ini_Hs            =   2.0;
    cnf_p.ini_Pu            =   1;
    cnf_p.init_rou          =   1e-2;
    cnf_p.initial_param_fit_feedback_flag = 0; %activate the updating of the initial parameters per fitted waveform either using the sliding window or with acumulation of previous info
    cnf_p.ini_Hs_rou_sliding_win_opt = 1;  % activate the averaging of Hs/roughness for next initial seed 
    cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold = 3; %number of times std for removal (it is also used for comparison with the accumulated one)
    cnf_p.ini_Hs_rou_sliding_win_size = 50; 
    
    %------------------ Two-step fitting ----------------------------------
    cnf_p.two_step_fitting=0; %used to perform first a fitting on SWH and if 
    cnf_p.two_step_fitting_COR_threshold_rou=90.0; %threshold below which second fitting on roughness is used
    
    %------------------ Range indexation ----------------------------------
    cnf_p.range_index_method='resampling'; % 'conventional': when creating the range vector used in fitting use directly 1:1:N_samples (even with ZP)
                                             % 'resampling': creating the
                                             % range vector as a resampling
                                             % of the original one (if there is a resampling)
    
    %------------------ Fitting function --------------------------------------
    cnf_p.fitting_fun_type='lsq'; % Define the fitting procedure used
    
    %String:
    %-'lsq': Using the lsqcurvfit based on levegender-Marquard
    %-'fmin': Using the fminsearch
    %algorithm
    switch cnf_p.fitting_fun_type
        case 'lsq'
            %trust-region-reflective
            %levenberg-marquardt
            %'MaxIter',400
            cnf_p.lsq_algorithm      = 'trust-region-reflective';
            %cnf_p.fitting_options    = optimset('Algorithm',cnf_p.lsq_algorithm,'FinDiffType','central','TolFun',1e-8,'TolX',1e-8,'Display','off','UseParallel','always');
            cnf_p.fitting_options    = optimset('Algorithm',cnf_p.lsq_algorithm,'Display','off','UseParallel','always');
            cnf_p.fitting_options_lb = [0,0,0]; %lower boundaries
            cnf_p.fitting_options_ub = []; %upper boundaries
        case 'fmin'
            cnf_p.fitting_options     =   optimset('Display','off','TolFun',1e-7,'TolX',1e-7);
            cnf_p.fitting_options_lb = []; %lower boundaries
            cnf_p.fitting_options_ub = []; %upper boundaries
    end
    
    %--------------------------------------------------------------------------
    %------------------ VERBOSE PARAMETERS-------------------------------------
    %--------------------------------------------------------------------------
    cnf_p.plot_fits_flag    =   0;              % 0 - not plot; 1 - plot fit results
    cnf_p.plot_fits_lat_range = [-91,91]; %latitude range of waveforms fit to be plotted (by default between -91 and 91 to force all of them to be plotted)
    cnf_p.plot_fits_downsampling =100;
    
end

%% -------------------- GEOPHYSICAL CORRECTIONS ----------------------------
% -------------------------------------------------------------------------
cnf_p.geo_corr_application_flag=1; 
cnf_p.force_geocorr_surf_type=1; % 0: use the geophysical corrections depending on the flag surface type included in the L1B product                          
                                 % 1: use the same group of geophysical corrections to all the records independently of the flag of the surface type: forcing the type surface
cnf_p.product_type_surface='open_ocean'; % open_ocean, sea_ice, land_ice ...                                 


%% -------------------- OUTPUT PRODUCT -------------------------------------
% -------------------------------------------------------------------------
%Inherete from Sentinel-3 file formating
cnf_p.Product_type='LAN'; % String: 'LAN' for land products
                          %'WAT' for water products
cnf_p.write_output=1;
cnf_p.output_product_format='nc'; % nc: NetCDF, mat: Matlab
cnf_p.nc_name_surface_height='shh'; %string indicating which name to be used in the netcdf for the surface height: SSH (sea surface height) or SH (surface height)

cnf_p.optional_ext_file_flag =0;
%cnf_p.file_ext_string='PTR_sigmaa_new_sigmar_cris_Lz_Ly_ZP_rg_idx_conv';
cnf_p.file_ext_string='L2_looks_by_Doppler';


%% ---------------- VERBOSE OPTIONS ---------------------------------------
cnf_p.visible_figures=0;


                                         