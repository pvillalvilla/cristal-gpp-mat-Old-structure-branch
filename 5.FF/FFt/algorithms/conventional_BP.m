%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Sentinel-6
%
% ---------------------------------------------------------
% Objective: Performs the core processing of the Fully Focussed time-domain
% processor based on simple conventional backprojection algorithm for each
% along-track surface input: every single surface consists of a set of
% surfaces (equal to number of range bins) to be focused across-track to generate a single waveform 
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
%  time_sar_ku_pulse:       time/datation of each pulse
%  alt_rate_sar_sat_pulse:  height rate of the satellite per pulse
%  x/y/z_surf:              x/y/z positions of surfaces
%  win_delay_surf:          window delay associated to surfaces  
%  time_surf:               datation time of surfaces
%  N_samples_az:            Number of azimuth input samples or pulses
%  N_samples_rg:            Number of range samples
%  cnf:                     Configuration parameters structure
%  chd:                     Characterization parameters structure
%  cst:                     Constant parameters structure
%
% OUTPUTs:
%  wfm_AC:                  Fully focused or azimuth compressed waveforms
%                           for each surface
%
%
% COMMENTS: Current version performs the backprojection process for each input
% surface and for the specified integration time integration_time
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%            Albert Garcia / isardSAT
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [wfm_AC,ProcessID] = conventional_BP ...
                                (wfm_RC,...
                                delay_matrix,...                                
                                win_delay_sar_ku_pulse,...
                                x_sar_sat_pulse, y_sar_sat_pulse, z_sar_sat_pulse,...
                                lat_sar_sat_pulse, lon_sar_sat_pulse, alt_sar_sat_pulse,...
                                pitch_sar_pulse, roll_sar_pulse, yaw_sar_pulse,...
                                x_vel_sat_sar_pulse, y_vel_sat_sar_pulse, z_vel_sat_sar_pulse,...
                                time_sar_ku_pulse,...
                                x_surf, y_surf, z_surf,...
                                time_surf,...                                
                                antenna_pattern_impact_flag,...
                                interpolation_method_backprojection_flag,...
                                c_cst_chd, pi_cst_chd, flat_coeff_cst_chd, semi_major_axis_cst_chd,...
                                wv_length_ku_chd, bw, Tp,...
                                antenna_beamwidth_along_track_ku, antenna_beamwidth_across_track_ku,...
                                integration_time,varargin)
   %% Optional arguments
    N_pulses = length(wfm_RC(:,1));
    p = inputParser;
    %p.addParamValue('matrix_approach',1);
    p.addParamValue('ProcessID_pulse',zeros(1,N_pulses));
    %p.addParamValue('i_sample_start_chd',-1);
    p.parse(varargin{:});
    ProcessID_pulse=p.Results.ProcessID_pulse;    
    %i_sample_start_chd=p.Results.i_sample_start_chd;    
    clear p;
    
    %% CONSIDERING ONLY THE PORTIONS OF INTEREST DEFINED BY INTEGRATION TIME           
    
    indexes_integration = logical(ones(1,length(wfm_RC(:,1))));
    
    if integration_time ~= -1
        indexes_integration = (time_sar_ku_pulse>=(time_surf-integration_time/2) & time_sar_ku_pulse<=(time_surf+integration_time/2));
    end
    
    if any(ProcessID_pulse(indexes_integration)==1)
        ProcessID = 1; %RMC
    else
        ProcessID = 0; %RAW
    end
    
    wfm_RC                      = wfm_RC(indexes_integration,:);
    delay_matrix                = delay_matrix(indexes_integration,:);
    win_delay_sar_ku_pulse      = win_delay_sar_ku_pulse(indexes_integration);
    x_sar_sat_pulse             = x_sar_sat_pulse(indexes_integration);
    y_sar_sat_pulse             = y_sar_sat_pulse(indexes_integration);
    z_sar_sat_pulse             = z_sar_sat_pulse(indexes_integration);
    lat_sar_sat_pulse           = lat_sar_sat_pulse(indexes_integration);
    lon_sar_sat_pulse           = lon_sar_sat_pulse(indexes_integration);
    alt_sar_sat_pulse           = alt_sar_sat_pulse(indexes_integration);
    pitch_sar_pulse             = pitch_sar_pulse(indexes_integration);
    roll_sar_pulse              = roll_sar_pulse(indexes_integration);
    yaw_sar_pulse               = yaw_sar_pulse(indexes_integration);
    x_vel_sat_sar_pulse         = x_vel_sat_sar_pulse(indexes_integration);
    y_vel_sat_sar_pulse         = y_vel_sat_sar_pulse(indexes_integration);
    z_vel_sat_sar_pulse         = z_vel_sat_sar_pulse(indexes_integration);
    
    %% Variables definition
    N_total_pulses              = length(wfm_RC(:,1));
    N_total_across_surf         = length(x_surf);
    %output waveform
    wfm_AC                      = zeros(1,N_total_across_surf);
    antenna_pattern_2W_comp     = zeros(1,N_total_pulses);
    

    %% BACKPROJECTION 
    for i_surf_across  = 1 : N_total_across_surf
        disp(strcat('Across_surf #=',num2str(i_surf_across)));
        for i_aperture = 1 : N_total_pulses
            % SLANT RANGE/DELAY COMPUTATION TO EACH PIXEL LOCATION
             vector_to_surface = [x_surf(i_surf_across), y_surf(i_surf_across), z_surf(i_surf_across)] - ...
                [x_sar_sat_pulse(i_aperture), y_sar_sat_pulse(i_aperture), z_sar_sat_pulse(i_aperture)];
            norm_vector_to_surface = norm(vector_to_surface);
            tau_slant_range = norm_vector_to_surface*2/c_cst_chd;
            
            
            % RANGE SHIFT DUE TO DOPPLER EFFECT (SATELLITE MOVING DURING TX/RX)
            vel_sat = [x_vel_sat_sar_pulse(i_aperture),y_vel_sat_sar_pulse(i_aperture),z_vel_sat_sar_pulse(i_aperture)];
            norm_vel_sat = norm(vel_sat);
            beam_angle = acos (dot(vector_to_surface,vel_sat) / (norm_vector_to_surface * norm_vel_sat));
            
            % Doppler related
            f_d = 2/wv_length_ku_chd.*norm(vel_sat).*cos(beam_angle);
            
            % Delay related change
            delay_Doppler = f_d/(bw/Tp);
                        
            
            % COMBINE BOTH DELAYS
            tau_delay_pixel_time_apert = tau_slant_range + delay_Doppler;
            
            
            % INTERPOLATION OF RANGE COMPRESSED DATA FOR FAST TIME DELAY
%             wfm_rf_comp_interp= interp1(delay_matrix(i_aperture,:),wfm_RC(i_aperture,:),tau_delay_pixel_time_apert,interpolation_method_backprojection_flag);
% 
%             if isnan(wfm_rf_comp_interp)
%                 wfm_rf_comp_interp=complex(0,0);
%             end
            
            
            %look for the closest range sample in the data
            if tau_delay_pixel_time_apert >= min(delay_matrix(i_aperture,:)) && tau_delay_pixel_time_apert <= max(delay_matrix(i_aperture,:)) 
                [~,rg_bin_int] = min(abs(delay_matrix(i_aperture,:)-tau_delay_pixel_time_apert)); 
                wfm_rf_comp_interp=wfm_RC(i_aperture,rg_bin_int);
            else
                wfm_rf_comp_interp=complex(0,0);
            end
            
            
            % COMPENSATION OF THE ANTENNA PATTERN
            if antenna_pattern_impact_flag
                % Compute the antenna weighting based on the position of TRP
                % and the antenna location and attitude
                
                % Conversion from ECEF to satellite framework (x || v_sat & z ||
                % nadir):
                % 1) Move location on-ground to satellite position: translation
                x_trp_sat        = x_surf(i_surf_across)-x_sar_sat_pulse(i_aperture);
                y_trp_sat        = y_surf(i_surf_across)-y_sar_sat_pulse(i_aperture);
                z_trp_sat        = z_surf(i_surf_across)-z_sar_sat_pulse(i_aperture);
                
                % 2) Transformation to satellite system coordinate
                % 2.2) Compute the nadir vector in sat coordinates
                % 2.2.1) Subsatellite point in ECEF
                xyz_nadir_ECEF   = lla2ecef([lat_sar_sat_pulse(i_aperture),...
                    lon_sar_sat_pulse(i_aperture),...
                    alt_sar_sat_pulse(i_aperture)-win_delay_sar_ku_pulse(i_aperture)],...
                    flat_coeff_cst_chd,semi_major_axis_cst_chd);%
                
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
                % coordindates (antenna is in plane x-y)
                xyz_trp_new_S6   = ([x_new_in_S6; y_new_in_S6; z_new_in_S6]*...
                    [x_trp_sat; y_trp_sat; z_trp_sat]);
                
                % Impact of pitch, roll and yaw of the antenna
                alpha            = yaw_sar_pulse(i_aperture); %yaw in degrees
                beta             = pitch_sar_pulse(i_aperture); %pitch in degrees
                gamma            = roll_sar_pulse(i_aperture); %roll in degrees
                rotation_matrix  = rotx(gamma)*roty(beta)*rotz(alpha);
                xyz_trp_new_S6   = (rotation_matrix*xyz_trp_new_S6).';
                
                norm_xyz_trp_S6  = norm(xyz_trp_new_S6);
                
                % Conversion to u and v coordinates
                u                = xyz_trp_new_S6(1)/norm_xyz_trp_S6;
                v                = xyz_trp_new_S6(2)/norm_xyz_trp_S6;
                
                
                % Antenna pattern impact estimation
                
                antenna_pattern_2W_comp(i_aperture) = modelled_antenna_pattern(u,v,u,v,...
                    0,0,...
                    antenna_beamwidth_along_track_ku,...
                    antenna_beamwidth_across_track_ku,...
                    antenna_beamwidth_along_track_ku,...
                    antenna_beamwidth_across_track_ku);
                
            else
                antenna_pattern_2W_comp(i_aperture) = 1.0;
            end
            
            wfm_rf_comp_interp = wfm_rf_comp_interp.*1/(antenna_pattern_2W_comp(i_aperture));
            
            %% PHASE COMPESNATION + COHERENT INTEGRATION
            wfm_AC(i_surf_across)=wfm_AC(i_surf_across)+...
                wfm_rf_comp_interp*exp(1i*2*pi_cst_chd*c_cst_chd/wv_length_ku_chd*tau_delay_pixel_time_apert);
                        
        end % end of azimuth positions within synthetic aperture radar

    end % different across-track surfaces
                
    end
    
	function [antenna_pattern] = modelled_antenna_pattern(u_TX,v_TX,u_RX,v_RX,...
        G_TX,G_RX,bw_al_TX,bw_ac_TX,bw_al_RX,bw_ac_RX)
    K = 0.5*log(2);
    % antenna pattern TX
    antenna_pattern_TX = 10^(G_TX/20).*exp(-K.*(2.*asin(u_TX)/bw_al_TX).^2).*...
        exp(-K.*(2.*asin(v_TX)/bw_ac_TX).^2);
    
    % antenna pattern RX
    antenna_pattern_RX = 10^(G_RX/20).*exp(-K.*(2.*asin(u_RX)/bw_al_RX).^2).*...
        exp(-K.*(2.*asin(v_RX)/bw_ac_RX).^2);
    
    % antenna pattern 2-W
    antenna_pattern = antenna_pattern_TX.*antenna_pattern_RX;
    
    
    end