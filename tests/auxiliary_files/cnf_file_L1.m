% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% GGPICE 
% ---------------------------------------------------------
% Objective: Define Configuration parameters
% 
% INPUTs : - 
% OUTPUTs: -
%
% ----------------------------------------------------------
% Author:    Albert Garcia / isardSAT
%
% Version  record
% 1.0 2018/07/12 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnf.mode            ='SAR';
cnf.processing_mode ='SAR';
cnf.onboard_reversion_flag=1;

cnf.verify_L0   = 0;
cnf.verify_L1A  = 0;
cnf.verify_L1BS_internal  = 0;
cnf.verify_L1B = 0;
% CHAINS
cnf.writting_flag = [1 0 1 0 0 0 1]; % L1A, L1BS_HR, L1B_HR, L1B_VHR_time, L1B_VHR_frequency, LROS, KML
cnf.run_L1A     = 0;
cnf.run_DDP     = 0;
cnf.run_FFtime  = 0;
cnf.run_FFfreq  = 0;
cnf.run_LROS 	= 0;
cnf.run_L2      = 1;
cnf.run_swath   = 0;

%% Filtering 
cnf.mask_flag = 0;

%% Processing and testing

%% PRELIMINARY DATATION
cnf.move_CoG2ant_alongtrack = 0;
%% PRELIMINARY WINDOW DELAY
cnf.window_delay_source = 1; % 0. Compute it from H0, COR2. 1.Take it from L0
%% FINAL BURST DATATION
cnf.height_rate_application = 0;
%% ONBOARD REVERSION

cnf.height_rate_application=0;
cnf.FAI_application=0;

%% INSTRUMENT GAIN

%% WAVEFORMS CORRECTION

cnf.CAL2_flag=0;
cnf.CAL1p2p_flag=0;

%% SURFACE LOCATION


cnf.N_surf_samp_interpol = 20;
cnf.N_points_spline_surf = 10;
cnf.smooth_fact_surf_pos = 0;
cnf.N_points_spline_orbit = cnf.N_points_spline_surf;
cnf.smooth_fact_orbit_pos = 0;
cnf.smooth_fact_orbit_vel = 0;
cnf.accu_orbit_alt = 0;
cnf.accu_orbit_lat = 0.0001;
cnf.sigma_alt_surf_geoloc_th = 100 * cnf.N_surf_samp_interpol; %SURFACE LOCATIONS (HIGH OR LOW ROUGHNESS)

cnf.create_surfaces_backward = 0;
cnf.create_surfaces_foreward = 0;
cnf.backward_margin = 32+16;
cnf.foreward_margin = 32+16;

cnf.trp_flag            = 0; % TRP PROCESSING
cnf.trp_flag_method     = 1;

if (cnf.trp_flag)
    % Simulated PT PICE
    lla = ecef2lla([chd.x_trp, chd.y_trp, chd.z_trp],cst.flat_coeff,cst.semi_major_axis);
    
     chd.lat_trp = lla(1);
     chd.lon_trp = lla(2);
     chd.alt_trp = lla(3);

    
    % Crete TRP
%      cnf.lat_trp = 35.3379302808;
%      cnf.lon_trp = 23.7795182869;
%      cnf.alt_trp = 1048.8184;
    % Svalbard TRP
%      lat_trp = 78.23052306;
%      lon_trp = 15.39376997;
%      alt_trp = 492.772;
end

%% BEAM ANGLES

%% AZIMUTH PROCESSING

cnf.sigma_alt_surf_th   = 1000; %AZIMUTH PROCESSING (HIGH OR LOW VARIABILITY)
cnf.zp_fact_azimut      = 1;
cnf.hamming_window      = 0;
cnf.force_exact_method  = 0;
if(cnf.hamming_window)
    %compute the gain scaling of the hamming window as defined by dinardo
%     hamm = hamming(chd.N_pulses_burst+1); %we want the hamming centered in sample #0 (if the samples go from -32 to 31)
%     hamm_win = hamm(1:chd.N_pulses_burst);
%     clear hamm;
    %Hamming according to Dinardo
    c1=0.08;
    c2=0.92;
    cnf.hamm_win=(c1+c2.*(cos(pi.*(0:1:chd.N_pulses_burst-1)./(chd.N_pulses_burst-1)-pi/2)).^2).';
%     %Unitary energy window
    %gain_scale_hamm=-10*log10(sum(hamm_win.^2));
    %Dinardo option
    cnf.gain_scale_hamm=20*log10(mean(cnf.hamm_win));
    
else
    cnf.gain_scale_hamm=0;
    cnf.hamm_win=[];
end

%% GEOMETRY CORRECTION

cnf.win_delay_ref = 0;       %0 = normal alignment    
                             %1 = optimal alignment for ICE (wrt max beam accumulated power) CR2
                             %2 = force the alignment to the closest of a given window_delay: Lakes
                             %3 = force to the first one.
                             %4 = force to the minimum one.
                             if(cnf.win_delay_ref==2)
                                 cnf.elevation_ref = 1048.8184; %force the window delay to a reference elevation
                             end
cnf.apply_doppler = 1;
cnf.apply_slant = 1;
cnf.apply_wd = 1;

%% RANGE COMPRESSION
cnf.zp_fact_range = 1;

cnf.coherence_threshold = 10;
cnf.stack_interferometry_flag = 1;
cnf.window_rg = 0;


cnf.extend_window = 0;
if(cnf.window_rg)
    %compute the gain scaling of the hamming window as defined by dinardo
%     hamm = hamming(chd.N_pulses_burst+1); %we want the hamming centered in sample #0 (if the samples go from -32 to 31)
%     hamm_win = hamm(1:chd.N_pulses_burst);
%     clear hamm;
    %Hamming according to Dinardo
    c1=0.08;
    c2=0.92;
    cnf.window_range=(c1+c2.*(cos(pi.*(0:1:chd.N_samples-1)./(chd.N_samples-1)-pi/2)).^2);
%     %Unitary energy window
    %gain_scale_hamm=-10*log10(sum(hamm_win.^2));
    %Dinardo option
    cnf.gain_scale_win_rg=20*log10(mean(cnf.window_range));
    
else
    cnf.gain_scale_win_rg=0;
    cnf.window_range=[];
end

%% MULTILOOKING and RETRACKERS

cnf.compute_L1_stadistics = 1;
cnf.N_beams_sub_stack = 5;
cnf.range_index1 = 12;
cnf.range_index2 = 16;
cnf.k_noise_top = 7;
cnf.k_noise_floor = 0;
cnf.use_zeros = 1;
cnf.apply_stack_mask=1;
cnf.delete_ambiguities=0;
cnf.coherency_mask = 1;

cnf.avoid_wrapps = 0;

cnf.avoid_beams_mask_allzeros=1;

chd.mask_look_angles_CR2_max=0.6*pi/180;
chd.mask_look_angles_CR2_min=-0.6*pi/180;
cnf.mask_look_angles_CR2=0;

cnf.avoid_noisy_beams=0; % waveforms of outer noisy beams in the stack are discarded 
                         % following Sentinel-3 baseline: threshold=noise_mean+3*std_noise       
cnf.param_std_noise_thres=3.0; %number of times the standard deviation of the noise used in the thresholding
cnf.noise_estimation_window=[12 16]; %first and last non-zero padded noise samples
cnf.method_noise_est='after_geom_corr'; % String: 
                       % 'after_geom_corr': noise estimation is performed
                       % right after the geometry corrections as done in
                       % Sentinel-3
                       % 'before_geom_corr': noise estimation is performed
                       % right before the geometry corrections as done in
                       % Sentinel-3                      
% cnf.method_noise_est='Beam_based'; % String: 
%                        % 'Beam_based': for each beam its related noise
%                        % statistics are estimated and accordingly applied
%                        % 'Stack_based': noise statistics estimated using
%                        % all the samples over the stack for the defined
%                        % noise estimation window --> thresholding for
%                        % different beams based on a single mean and std
% cnf.method_noisy_thresholding='Average_based'; % String:
%                                 % 'Peak_based': Compare peak with threshold
%                                 % 'Average_based': Compare average with
%                                 % threshold

            

%% not used?
cnf.clock_unit_conv = 256;
cnf.cor2_bit_shift = 4;

%% Product writing
cnf.include_wfms_aligned=0;
% global netcdf_type
% netcdf_type='netcdf4'; %string indicating the type of netcdf 3 or 4: netcdf3 or netcdf4: 
% %Consequence on the uint not available in the netcdf3 definition need to
% %enlarge to int ouf double of bits to accomodate it --> enlarge the product

%% Verbosity: figures pop up or not
cnf.figures_visible=0;

%% Configuration parameters FFt
% computation grid
cnf.FFt.spacing_al_surfaces = 0.10; %spacing between the surfaces along-track
cnf.FFt.surf_dist_trp = 3.0; % +- distance around transponder used to filter out the surfaces to be used in the processing

% Kernel of backprojection
cnf.FFt.antenna_pattern_impact_flag = 1; % compensate the effect of antenna pattern (to be implemented when paralelizing)
cnf.FFt.T_integration = 1.95; % %if set to -1 take all the input data for processing
cnf.FFt.num_pools = 4; % to paralelize the processing of each surface
cnf.FFt.zp_ac_TRP = 8; % For transponder: zero-padding for across-track over the final image in TRP case
cnf.FFt.zp_al_TRP = 8; % For transponder: zero-padding for across-track over the final image in TRP case
