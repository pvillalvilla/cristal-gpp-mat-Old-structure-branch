%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: Performs the core processing of the Fully Focussed time-domain
% processor, backprojection algorithm
%
% Calling: 
%
% INPUTs:
%  wfm_cal_gain_corrected:  waveforms per pulse concatenated for all the bursts
%  T0_sar_pre_dat_pulse:    Sampling period per pulse
%  win_delay_sar_ku_pulse:  window delay of each pulse
%  x/y/z_sar_sat_pulse:     x/y/z positions of the satellite per pulse 
%  lat_sar_sat_pulse:       latitude of satellite per pulse 
%  lon_sar_sat_pulse:       longitude of satellite per pulse
%  alt_sar_sat_pulse:       altitude of the satellite per pulse
%  pitch_sar_pulse:         pitch of the satellite per pulse
%  roll_sar_pulse:          roll of the satellite per pulse
%  yaw_sar_pulse:           yaw of the satellite per pulse
%  x/y/z_vel_sat_sar_pulse: x/y/z satellite velocity comp. per pulse 
%  alt_rate_sar_sat_pulse:  height rate of the satellite per pulse
%  x/y/z_surf:              x/y/z positions of surfaces
%  win_delay_surf:          window delay associated to surfaces  
%  time_surf:               datation time of surfaces
%  N_samples_az:            Number of azimuth input samples or pulses
%  N_samples_rg:            Number of range samples
%  indexes_integration:     Indexes where to compute the algorithm
%  cnf:                     Configuration parameters structure
%  chd:                     Characterization parameters structure
%  cst:                     Constant parameters structure
%
% OUTPUTs:
%  wfm_AC:                  Fully focused or azimuth compressed waveforms
%                           for each surface (integrated)
%  wfm_RC_phase_corr_masked:Fully focused or azimuth compressed waveforms
%                           for each surface (non integrated)
%
% COMMENTS: Current version performs the backprojection process for each input
% surface and for the specified integration time cnf.FFt.T_integration
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%            Albert Garcia / isardSAT
%            Ferran Gibert / isardSAT
% ----------------------------------------------------------
% v 1.1 2019/06/10 Removed Process_ID 
% v 1.2 2019/06/20 - Enable interferometry processing:
%                   'indexes_integration' variable moved to input variable
%                   'wfm_RC_phase_corr_masked' variable set as output
%                   variable
%                   'time_sar_ku_pulse' removed as it is no longer required
%                   'time_surf' removed as it is no longer required
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [wfm_AC, wfm_RC_phase_corr_masked, N_av_SL] = backprojection_kernel_freq_domain(...
                       wfm_cal_gain_corrected,...
                       T0_sar_pre_dat_pulse,...
                       win_delay_sar_ku_pulse,...
                       x_sar_sat_pulse,y_sar_sat_pulse,z_sar_sat_pulse,...
                       lat_sar_sat_pulse,lon_sar_sat_pulse,alt_sar_sat_pulse,...
                       pitch_sar_pulse,roll_sar_pulse,yaw_sar_pulse,...
                       x_vel_sat_sar_pulse,y_vel_sat_sar_pulse,z_vel_sat_sar_pulse,...
                       x_surf,y_surf,z_surf,win_delay_surf,...
                       N_samples_az,N_samples_rg,...
                       indexes_integration, ...
                       cnf,cst,chd,...                                           
                       varargin)  
                   
    %% Optional arguments
    p = inputParser;
    %p.addParamValue('matrix_approach',1);
    p.addParamValue('ProcessID_pulse',zeros(1,N_samples_az));
    p.addParamValue('i_sample_start_chd',-1);
    p.parse(varargin{:});
    ProcessID_pulse=p.Results.ProcessID_pulse;    
    i_sample_start_chd=p.Results.i_sample_start_chd;    
    clear p;
    
    
    %% CONSIDERING ONLY THE PORTIONS OF INTEREST DEFINED BY INTEGRATION TIME       
    % done in this way to allow paralelization of the different surfaces
    
    if any(ProcessID_pulse(indexes_integration)==1)
        ProcessID = 1; %RAW
    else
        ProcessID = 0; %RAW
    end
    
    wfm_cal_gain_corrected = wfm_cal_gain_corrected(indexes_integration,:);
    T0_sar_pre_dat_pulse = T0_sar_pre_dat_pulse(indexes_integration);
    win_delay_sar_ku_pulse =win_delay_sar_ku_pulse(indexes_integration);
    x_sar_sat_pulse = x_sar_sat_pulse(indexes_integration);
    y_sar_sat_pulse = y_sar_sat_pulse(indexes_integration);
    z_sar_sat_pulse = z_sar_sat_pulse(indexes_integration);
    lat_sar_sat_pulse = lat_sar_sat_pulse(indexes_integration);
    lon_sar_sat_pulse = lon_sar_sat_pulse(indexes_integration);
    alt_sar_sat_pulse = alt_sar_sat_pulse(indexes_integration);
    pitch_sar_pulse = pitch_sar_pulse(indexes_integration);
    roll_sar_pulse = roll_sar_pulse(indexes_integration);
    yaw_sar_pulse = yaw_sar_pulse(indexes_integration);
    x_vel_sat_sar_pulse = x_vel_sat_sar_pulse(indexes_integration);
    y_vel_sat_sar_pulse = y_vel_sat_sar_pulse(indexes_integration);
    z_vel_sat_sar_pulse = z_vel_sat_sar_pulse(indexes_integration);
    
    N_samples_az = length(find(indexes_integration));
    
    tau_delay_pixel_time_apert  = zeros(1,N_samples_az);
    delay_Doppler               = zeros(1,N_samples_az);
    beam_angle                  = zeros(1,N_samples_az);
    f_d                         = zeros(1,N_samples_az);
    slant_range_corr            = zeros(1,N_samples_az);
    wd_corr                     = zeros(1,N_samples_az);
    doppler_corr                = zeros(1,N_samples_az);
    
    wfm_cal_gain_corrected_corrected = zeros(N_samples_az,N_samples_rg);

    %% BACKPROJECTION 
    cnf.FFt.matrix_approach = 0; %Forced to call looping per pulse within integration time as Antenna pattern compensation in matrix formulation not available    
    if cnf.FFt.matrix_approach
        %Matrix notation
        % ------------------ Range frequency ----------------------------------
        f_rg_norm = ones(N_samples_az,1)*((0:1:N_samples_rg-1)/N_samples_rg-1/2d0);
        T0_matrix = (T0_sar_pre_dat_pulse.'*ones(1,N_samples_rg));
        
        
        %% SLANT RANGE/DELAY COMPUTATION
        vector_to_surface = ones(N_samples_az,1)*[x_surf, y_surf, z_surf] - ...
            [x_sar_sat_pulse.', y_sar_sat_pulse.', z_sar_sat_pulse.'];
        norm_vector_to_surface = sqrt(sum(vector_to_surface.^2,2));
        tau_delay_pixel_time_apert = (norm_vector_to_surface*ones(1,N_samples_rg))*2/cst.c;
        
        % slant range correction in samples
        slant_range_corr = (win_delay_surf-tau_delay_pixel_time_apert)./T0_matrix; 
        
        % for AZ phase compensation
        tau_delay_pixel_time_apert = (tau_delay_pixel_time_apert(:,1)).';


        %% WINDOW DELAY MISALIGNMENT
        % correction in samples
        wd_corr = (-(win_delay_surf-(win_delay_sar_ku_pulse.')*ones(1,N_samples_rg)))./T0_matrix;

        %% STOP & GO CORRECTION: FAST-TIME EFFECT: DOPPLER
        % Compute beam angle for each aperture time and surface
        vel_sat = [x_vel_sat_sar_pulse.',y_vel_sat_sar_pulse.',z_vel_sat_sar_pulse.'];
        norm_vel_sat = sqrt(sum(vel_sat.^2,2));
        beam_angle = acos (dot(vector_to_surface,vel_sat,2)./(norm_vector_to_surface.*norm_vel_sat));
        
        % Doppler related
        f_d = (2/chd.wv_length.*norm_vel_sat.*cos(beam_angle))*ones(1,N_samples_rg);
        
        % Delay related change
        delay_Doppler = f_d./(chd.bw/chd.pulse_length);
                
        % correction in samples
        doppler_corr = delay_Doppler./T0_matrix;
        
        %% APPLICATION OF CORRECTIONS
        wfm_cal_gain_corrected_corrected(:,:) = wfm_cal_gain_corrected(:,:).*...
            exp(1i*2*cst.pi*(slant_range_corr+doppler_corr+wd_corr).*f_rg_norm); %
        
        %from matrix to vector for mask construction
        slant_range_corr = slant_range_corr(:,1);
        wd_corr          = wd_corr(:,1);
        doppler_corr     = doppler_corr(:,1);
        %% COMPENSATION OF THE ANTENNA PATTERN
        % To be implemented (not an easy translation to matrix notation)
    else        
        for i_aperture=1:N_samples_az
            
            % ------------------ Norm. Range frequency ----------------------------------
            f_rg_norm = ((0:1:N_samples_rg-1)*1/N_samples_rg-1/2d0);
            
            
            %% SLANT RANGE/DELAY COMPUTATION TO TRP LOCATION
            vector_to_surface = [x_surf, y_surf, z_surf] - ...
                [x_sar_sat_pulse(i_aperture), y_sar_sat_pulse(i_aperture), z_sar_sat_pulse(i_aperture)];
            norm_vector_to_surface = norm(vector_to_surface);
            tau_delay_pixel_time_apert(i_aperture) = norm_vector_to_surface*2/cst.c;
            
            slant_range_corr(i_aperture)  = (win_delay_surf-tau_delay_pixel_time_apert(i_aperture))./T0_sar_pre_dat_pulse(i_aperture) ;
            

            %% WINDOW DELAY MISALIGNMENT            
            wd_corr(i_aperture) = (-(win_delay_surf-win_delay_sar_ku_pulse(i_aperture)))./T0_sar_pre_dat_pulse(i_aperture);
            
            
            %% STOP & GO CORRECTION: FAST-TIME EFFECT: DOPPLER
            % Compute beam angle for each aperture time and surface
            vel_sat = [x_vel_sat_sar_pulse(i_aperture),y_vel_sat_sar_pulse(i_aperture),z_vel_sat_sar_pulse(i_aperture)];
            norm_vel_sat = norm(vel_sat);
            beam_angle(i_aperture) = acos (dot(vector_to_surface,vel_sat) / (norm_vector_to_surface * norm_vel_sat));
            
            % Doppler related
            f_d(i_aperture) = 2/chd.wv_length.*norm(vel_sat).*cos(beam_angle(i_aperture));
            
            % Delay related change
            delay_Doppler(i_aperture) = f_d(i_aperture)/(chd.bw/chd.pulse_length);
            
            doppler_corr(i_aperture) = delay_Doppler(i_aperture)./T0_sar_pre_dat_pulse(i_aperture);
            
            %% APPLICATION OF CORRECTIONS
%             wfm_cal_gain_corrected_corrected(i_aperture,:) = wfm_cal_gain_corrected(i_aperture,:).*...
%                 exp(1i*2*cst.pi*(slant_range_corr(i_aperture)+wd_corr(i_aperture)+doppler_corr(i_aperture)).*f_rg_norm);
             wfm_cal_gain_corrected_corrected(i_aperture,:) = wfm_cal_gain_corrected(i_aperture,:).*... % testing not apply doppler correct to TPT.11
                 exp(1i*2*cst.pi*(slant_range_corr(i_aperture)+wd_corr(i_aperture)).*f_rg_norm);      
%              wfm_cal_gain_corrected_corrected(i_aperture,:) = wfm_cal_gain_corrected(i_aperture,:).*... % testing not apply doppler correct & slant to TPT.11
%                  exp(1i*2*cst.pi*(wd_corr(i_aperture)).*f_rg_norm);                
            
            %% COMPENSATION OF THE ANTENNA PATTERN
            if cnf.FFt.antenna_pattern_impact_flag
                % Compute the antenna weighting based on the position of TRP
                % and the antenna location and attitude
                
                % Conversion from ECEF to satellite framework (x || v_sat & z ||
                % nadir):
                % 1) Move location on-ground to satellite position: translation
                % NOTE: Actually just calculation of relative position vector
                x_trp_sat        = x_surf-x_sar_sat_pulse(i_aperture);
                y_trp_sat        = y_surf-y_sar_sat_pulse(i_aperture);
                z_trp_sat        = z_surf-z_sar_sat_pulse(i_aperture);
                
                % 2) Transformation to satellite system coordinate
                % 2.2) Compute the nadir vector in sat coordinates
                % 2.2.1) Subsatellite point in ECEF
                xyz_nadir_ECEF   = lla2ecef([lat_sar_sat_pulse(i_aperture),...
                                            lon_sar_sat_pulse(i_aperture),...
                                            alt_sar_sat_pulse(i_aperture)-win_delay_sar_ku_pulse(i_aperture)],...
                                            cst.flat_coeff,cst.semi_major_axis);%
%                 xyz_nadir_ECEF   = lla2ecef([lat_sar_sat_pulse(i_aperture),...
%                                             lon_sar_sat_pulse(i_aperture),...
%                                             alt_sar_sat_pulse(i_aperture)-win_delay_sar_ku_pulse(i_aperture)*cst.c/2],... % we were substracting win delay (s) to alt(m)??? Transformed to distance
%                                             cst.flat_coeff,cst.semi_major_axis);%
                
                % 2.2.2) Subsatellite point in satellite coordindates
                xyz_nadir_in_sat = [xyz_nadir_ECEF(1)-x_sar_sat_pulse(i_aperture),...
                                    xyz_nadir_ECEF(2)-y_sar_sat_pulse(i_aperture),...
                                    xyz_nadir_ECEF(3)-z_sar_sat_pulse(i_aperture)];
                
                z_new_in_S6      = xyz_nadir_in_sat./norm(xyz_nadir_in_sat);
                
                % Compute the normalized vector || satellite velocity
                sat_vel_norm     = norm([x_vel_sat_sar_pulse(i_aperture),y_vel_sat_sar_pulse(i_aperture),z_vel_sat_sar_pulse(i_aperture)]);
                v_velocity_norm      = [x_vel_sat_sar_pulse(i_aperture)/sat_vel_norm,...
                    y_vel_sat_sar_pulse(i_aperture)/sat_vel_norm,...
                    z_vel_sat_sar_pulse(i_aperture)/sat_vel_norm];
                
                % Vector perpendicular to velocity and nadir
                w = cross(z_new_in_S6,v_velocity_norm);
                w = w./norm(w);
                
                % Vector tangent velocitat (new x-axis)
                x_new_in_S6 = cross(w,z_new_in_S6);
                x_new_in_S6 = x_new_in_S6/norm(x_new_in_S6);
                
                % 2.2.3) Compute the y_new_in_S6 (the y axis of S6 coordinate system)
                y_new_in_S6      = cross(z_new_in_S6,x_new_in_S6);
                
                % 2.2.4) Transformation of point target location in S6
                % coordinates (antenna is in plane x-y)
                xyz_trp_new_S6   = ([x_new_in_S6; y_new_in_S6; z_new_in_S6]*...
                                    [x_trp_sat; y_trp_sat; z_trp_sat]);
                
                % Impact of pitch, roll and yaw of the antenna
                alpha            = yaw_sar_pulse(i_aperture); %yaw in degrees
                beta             = pitch_sar_pulse(i_aperture); %pitch in degrees
                gamma            = roll_sar_pulse(i_aperture); %roll in degrees
                
                rotation_matrix  = rotx(gamma)*roty(beta)*rotz(alpha); % FGG 2019/6/3: THIS LINE REQUIRES SPECIFIC FUNCTIONS FROM PHASED ARRAY SYSTEM TOOLBOX
                rotx_gamma =  [1 0 0; 0 cosd(gamma) -sind(gamma); 0 sind(gamma) cosd(gamma)];
                roty_beta  =  [cosd(beta) 0 sind(beta); 0 1 0; -sind(beta) 0 cosd(beta)];
                rotz_alpha =  [cosd(alpha) -sind(alpha) 0; sind(alpha) cosd(alpha) 0; 0 0 1];
                rotation_matrix = rotx_gamma*roty_beta*rotz_alpha;
                
                xyz_trp_new_S6   = (rotation_matrix*xyz_trp_new_S6).';
                
%                 %JPL: code above was from S6 backprojection_bias, made the changes in the no bias script
%                 rotx_gamma =  [1 0 0; 0 cosd(gamma) sind(gamma); 0 -sind(gamma) cosd(gamma)]; % Matrix for positive rotation of gamma around X in S6 base (roll)
%                 roty_beta  =  [cosd(beta) 0 -sind(beta); 0 1 0; sind(beta) 0 cosd(beta)];   % Matrix for positive rotation of beta around Y in S6 base (pitch)
%                 rotz_alpha =  [cosd(alpha) sind(alpha) 0; -sind(alpha) cosd(alpha) 0; 0 0 1]; % Matrix for positive rotation of alpha around Z in S6 base (yaw(
%                 rotation_matrix = rotz_alpha*roty_beta*rotx_gamma; % 1st Roll, 2nd Pitch, 3rd Yaw
%                 xyz_trp_new_S6   = (rotation_matrix*xyz_trp_new_S6);
                
                norm_xyz_trp_S6  = norm(xyz_trp_new_S6);
                
                % Conversion to u and v coordinates
                u = xyz_trp_new_S6(1)/norm_xyz_trp_S6;
                v = xyz_trp_new_S6(2)/norm_xyz_trp_S6;
                
                
                % Antenna pattern impact estimation
                antenna_pattern_2W_comp(i_aperture) = modelled_antenna_pattern(u,v,u,v,...
                    0,0,...
                    chd.antenna_beamwidth_along_track,...
                    chd.antenna_beamwidth_across_track,...
                    chd.antenna_beamwidth_along_track,...
                    chd.antenna_beamwidth_across_track);
                               
            else
                antenna_pattern_2W_comp(i_aperture) = 1.0;
            end
            
            wfm_cal_gain_corrected_corrected(i_aperture,:) = wfm_cal_gain_corrected_corrected(i_aperture,:)./sqrt(antenna_pattern_2W_comp(i_aperture));
            
        end % end of azimuth positions within synthetic aperture radar
    end
    
    %% RANGE COMPRESSION
    clear fs;
    fs = mean(1./T0_sar_pre_dat_pulse);
    f_rg=((1:1:N_samples_rg)*fs/N_samples_rg-fs/2d0);
    if cnf.window_rg
        f_centroid = 0.0;
        w_rg = window_computation(cnf.flag_range_win_type,...
                    cnf.range_BW_to_process,f_rg,f_centroid,fs,...
                    'win_param',cnf.flag_range_win_param,...
                    'normalize',cnf.flag_range_win_normalize,...
                    'normalize_type',cnf.flag_range_win_normalize_type);
    else
        w_rg = zeros(1, N_samples_rg) + 1;        
    end
    % ------------ 12.2 Perform the range compression ---------------------
    [wfm_RC] = RC_processing (wfm_cal_gain_corrected_corrected, w_rg,cnf.zp_fact_range);
    %clear wfm_SRC;
    
    %% AZIMUTH PHASE COMPENSATION        
    wfm_RC_phase_corr = wfm_RC.*((exp(1i*2*cst.pi*cst.c/chd.wv_length*(tau_delay_pixel_time_apert))).'*ones(1,N_samples_rg.*cnf.zp_fact_range));
    
    %% MASKING 
    wfm_RC_phase_corr_masked = mask_FFS(wfm_RC_phase_corr,doppler_corr,slant_range_corr,wd_corr,ProcessID,...
                                        N_samples_rg,cnf.zp_fact_range,i_sample_start_chd);
                                    
    %% INTEGRATION and normalization by number of non-nan pulses 
    N_av_SL = numel(wfm_RC_phase_corr_masked(:,1)) - sum(isnan(wfm_RC_phase_corr_masked(:,1)));
    wfm_AC = nansum(wfm_RC_phase_corr_masked)/N_av_SL';
    %wfm_AC = sum(wfm_RC_phase_corr_masked(1:end,:),1);
    %plot(abs((wfm_AC).^2));hold on
    
    end
    
    function [antenna_pattern] = modelled_antenna_pattern(u_TX,v_TX,u_RX,v_RX,...
        G_TX,G_RX,bw_al_TX,bw_ac_TX,bw_al_RX,bw_ac_RX)
        K = 0.5*log(2);
        % antenna pattern TX
        antenna_pattern_TX = 10^(G_TX/10).*...
                                    exp(-K.*(2.*asin(u_TX)/bw_al_TX).^2).*...
                                    exp(-K.*(2.*asin(v_TX)/bw_ac_TX).^2);

        % antenna pattern RX
        antenna_pattern_RX = 10^(G_RX/10).*...
                                    exp(-K.*(2.*asin(u_RX)/bw_al_RX).^2).*...
                                    exp(-K.*(2.*asin(v_RX)/bw_ac_RX).^2);

        % antenna pattern 2-W
        antenna_pattern = antenna_pattern_TX.*antenna_pattern_RX;
    
    
    end
