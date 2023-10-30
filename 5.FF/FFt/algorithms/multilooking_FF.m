%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% CRISTAL
%
% ---------------------------------------------------------
% Objective: Perform the incoherent multilooking between the FF waveforms
% over a given final spatial resolution
%
% Calling: 
%
% INPUTs:
%  FF_wvfm:                 matrix with all waveforms from all surfaces:
%  each row contains the FF complex waveforms
%  lat_surf:                latitude of surfaces (at FF resolution)
%  lon_surf:                longitude of surfaces 
%  win_delay_surf:          window delay of surfaces
%  time_surf:               datation of the surface
%  x/y/z_vel_sat_surf:      z/y/z velocity components of the satellite for
%  each surface
%  alt_rate_sat_surf:       altitude rate of satellite for all surfaces
%  alt_sat_surf:            altitude of the satellite for all surfaces
%  pitch_sat_surf:          pitch of satellite for all surfaces
%  roll_sat_surf:           roll of satellite for all surfaces
%  yaw_sat_surf:            yaw of satellite for all surfaces
%  T0_surf:                 T0 (clock: inverse of sampling frequency) per surface
%  pri_surf:                Pulse repetition interval per surface

%
% OUTPUTs:
%
%  pow_ML_FF:                  matrix with all the ML power waveforms from all surfaces:
%                              each row contains the FF ML power waveform (waveforms have been aligned w.r.t first surface in the block before averaging)
%  pow_ML_FF_noalign:          matrix with all the ML power waveforms from all surfaces:
%                              each row contains the FF ML power waveform (waveforms have NOT been aligned w.r.t first surface in the block before averaging)
%  lat_surf_ML:                latitude of ML waveforms
%  lon_surf_ML:                longitude of ML waveforms
%  win_delay_surf_ML:          window delay of ML waveforms
%  time_surf_ML:               datation of ML waveforms
%  x/y/z_vel_sat_surf_ML:      z/y/z velocity components of the satellite for each ML waveform
%  alt_rate_sat_surf_ML:       altitude rate of satellite for ML waveforms
%  alt_sat_surf_ML:            altitude of the satellite for ML waveforms
%  pitch_sat_surf_ML:          pitch of satellite for ML waveforms
%  roll_sat_surf_ML:           roll of satellite for ML waveforms
%  yaw_sat_surf_ML:            yaw of satellite for ML waveforms
%  T0_surf_ML:                 T0 (clock: inverse of sampling frequency) per ML waveform
%  pri_surf_ML:                Pulse repetition interval per surface
%  N_looks:                    Number of averaged looks per surface
%
%
% COMMENTS: 
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%            Albert Garcia / isardSAT
%            Ferran Gibert / isardSAT 
%            Juan Pedro López-Zaragoza / isardSAT
% x.2 2023/03/28 FF_wvfm_comb pow_ML_FF_comb_noalign computation  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [pow_ML_FF, pow_ML_FF_noalign, pow_ML_FF_2, pow_ML_FF_2_noalign,pow_ML_FF_comb, pow_ML_FF_comb_noalign, phase_difference_ML_FF, wvfm_coh_ML_FF,... 
           lat_surf_ML,lon_surf_ML,win_delay_surf_ML, time_surf_ML,...
           x_vel_sat_surf_ML,y_vel_sat_surf_ML,z_vel_sat_surf_ML,...
           alt_rate_sat_surf_ML,alt_sat_surf_ML,...
           pitch_sat_surf_ML,roll_sat_surf_ML,yaw_sat_surf_ML,...
           num_aver_echoes_ML, ...
           T0_surf_ML,pri_surf_ML,N_looks,looks_flag]...
           ...
           = multilooking_FF ...
           (FF_wvfm, FF_wvfm_2,FF_wvfm_comb, phase_difference, wvfm_coh, ...
           lat_surf, lon_surf, win_delay_surf, time_surf,...        
           x_vel_sat_surf, y_vel_sat_surf, z_vel_sat_surf,...
           alt_rate_sat_surf, alt_sat_surf, ...       
           pitch_sat_surf, roll_sat_surf, yaw_sat_surf,...
           num_aver_echoes_surf, ...
           T0_surf, pri_surf,...
           out_spatial_res,...
           zp_fact_rg,...
           semi_major_axis_cst,semi_minor_axis_cst,c_cst,pi_cst,...
           processID,i_sample_start_chd,rmc_margin,processing_mode,ref_alignment)
        
    % Consider the multilooking based on actual distance on ground of
    % surfaces
    %% Initial variables
    N_surfaces_FF = length(FF_wvfm(:,1));%number of total FF surfaces separated with the theoretical expected resolution
    N_samples_rg  = length(FF_wvfm(1,:)); % number of samples of the range waveform (even if zero-padding has been considered)
    
    i_surf_ML_FF = 0; % index of the ML FF waveform
    i_ref_surf_FF = 1;
    i_sec_surf_FF = 2;
    %% Computation of the initial and final indexes of FF surfaces to be averaged     
    while i_ref_surf_FF <N_surfaces_FF
        %For a surface we compute the distance to a reference one until the
        %separation is larger or equal to the expected one over which to
        %average
        while i_sec_surf_FF <= N_surfaces_FF
            %disp(i_surf_ML_FF)
            [arclen,~]=distance(lat_surf(i_ref_surf_FF),lon_surf(i_ref_surf_FF),...
                lat_surf(i_sec_surf_FF),lon_surf(i_sec_surf_FF),...
                [semi_major_axis_cst sqrt(1-(semi_minor_axis_cst/semi_major_axis_cst).^2)]);
            
            if arclen>=out_spatial_res
                i_surf_ML_FF = i_surf_ML_FF+1;
                %keep with the surface providing the closest separation to
                %the expected output spatial resolution (between actual and previous one)                
                [arclen_prev,~]=distance(lat_surf(i_ref_surf_FF),lon_surf(i_ref_surf_FF),...
                    lat_surf(i_sec_surf_FF-1),lon_surf(i_sec_surf_FF-1),...
                    [semi_major_axis_cst sqrt(1-(semi_minor_axis_cst/semi_major_axis_cst).^2)]);
                
                if abs(arclen-out_spatial_res)<abs(arclen_prev-out_spatial_res)
                    %keep actual as last surface
                    last_surf_FF(i_surf_ML_FF)= i_sec_surf_FF;
                elseif abs(arclen-out_spatial_res)>abs(arclen_prev-out_spatial_res)
                    %keep the previous as last surface
                    last_surf_FF(i_surf_ML_FF)= i_sec_surf_FF-1;
                end
                init_surf_FF(i_surf_ML_FF)= i_ref_surf_FF;
                                

                %update the initial reference surface for the next averaging
                %block
                i_ref_surf_FF = last_surf_FF(i_surf_ML_FF)+1;
                i_sec_surf_FF = last_surf_FF(i_surf_ML_FF)+2;
                break;
            else
                if i_sec_surf_FF == N_surfaces_FF
                    %break both loops
                    i_surf_ML_FF = i_surf_ML_FF+1;
                    init_surf_FF(i_surf_ML_FF)= i_ref_surf_FF;
                    last_surf_FF(i_surf_ML_FF)= i_sec_surf_FF;
                    i_ref_surf_FF = i_sec_surf_FF+1;
                    i_sec_surf_FF = i_sec_surf_FF+2;
                    break;
                else
                    i_sec_surf_FF = i_sec_surf_FF+1;                
                end
                
            end            
        end
    end
    %In case there is only a single SLC surface to be averaged
    if i_ref_surf_FF == N_surfaces_FF
        i_surf_ML_FF = i_surf_ML_FF+1;
        init_surf_FF(i_surf_ML_FF)= N_surfaces_FF;
        last_surf_FF(i_surf_ML_FF)= N_surfaces_FF;
    end
    
    N_looks = last_surf_FF - init_surf_FF+1;
    N_total_ML_FF = i_surf_ML_FF; % total number of multilooked surfaces
    clear i_surf_ML_FF;
    
    % Flag to indicate if the last surface has less waveforms averaged than
    % the others because of the end of the dataset
    looks_flag(1:N_total_ML_FF-1)=1;
    if N_looks(N_total_ML_FF)<(N_total_ML_FF-1)
       looks_flag(N_total_ML_FF)=0;
    else
       looks_flag(N_total_ML_FF)=1;
    end
    
    %% Define output variables
    pow_ML_FF                 = zeros(N_total_ML_FF,N_samples_rg);
    pow_ML_FF_noalign         = zeros(N_total_ML_FF,N_samples_rg);
    %if(strcmp(processing_mode,'SIN'))
    pow_ML_FF_2                 = zeros(N_total_ML_FF,N_samples_rg);
    pow_ML_FF_2_noalign         = zeros(N_total_ML_FF,N_samples_rg);
    pow_ML_FF_comb              = zeros(N_total_ML_FF,N_samples_rg);
    pow_ML_FF_comb_noalign      = zeros(N_total_ML_FF,N_samples_rg);
    phase_difference_ML_FF      = zeros(N_total_ML_FF,N_samples_rg);
    wvfm_coh_ML_FF              = zeros(N_total_ML_FF,N_samples_rg);
    %end 
    lat_surf_ML               = zeros(1,N_total_ML_FF);
    lon_surf_ML               = zeros(1,N_total_ML_FF);
    win_delay_surf_ML         = zeros(1,N_total_ML_FF);
    time_surf_ML              = zeros(1,N_total_ML_FF);        
    x_vel_sat_surf_ML         = zeros(1,N_total_ML_FF);
    y_vel_sat_surf_ML         = zeros(1,N_total_ML_FF);
    z_vel_sat_surf_ML         = zeros(1,N_total_ML_FF);
    alt_rate_sat_surf_ML      = zeros(1,N_total_ML_FF);
    alt_sat_surf_ML           = zeros(1,N_total_ML_FF);       
    pitch_sat_surf_ML         = zeros(1,N_total_ML_FF);       
    roll_sat_surf_ML          = zeros(1,N_total_ML_FF);
    yaw_sat_surf_ML           = zeros(1,N_total_ML_FF);
    T0_surf_ML                = zeros(1,N_total_ML_FF);
    pri_surf_ML               = zeros(1,N_total_ML_FF);
    num_aver_echoes_ML        = zeros(1,N_total_ML_FF); 
    
    %% Compute the Multilooked waveforms & associated parameters
    %Assumed the reference surface is the first one in the block averaging
    %Alignment is based on phase ramp application of the specific delay
    
    %Normalized range frequency vector used for Phase ramp correction in
    %the range-frequency domain (to align waveforms)
    f_rg_norm = ((0:1:N_samples_rg-1)*1/N_samples_rg-1/2d0); %N_samples_rg includes any zero-padding in the data, 
                            %in the computation of the shift in samples
                            %zp_fact shall be considered to be aligned
    
    %loop over the total number of final ML FF waveforms
    for i_surf_ML_FF = 1: N_total_ML_FF
        
        %index to be used to identify the number of waveform within the
        %block of waveforms to be averaged
        i_surf_within_ML_block = 0;
        
        %---------- Move to the range-frequency domain via IFFT ------------
        %Consider only the data block of interest
        data_block = FF_wvfm(init_surf_FF(i_surf_ML_FF):last_surf_FF(i_surf_ML_FF),:);
        if(strcmp(processing_mode,'SIN'))
            data_block_2 = FF_wvfm_2(init_surf_FF(i_surf_ML_FF):last_surf_FF(i_surf_ML_FF),:);
            data_block_comb = FF_wvfm_comb(init_surf_FF(i_surf_ML_FF):last_surf_FF(i_surf_ML_FF),:);
            phase_difference_block = phase_difference(init_surf_FF(i_surf_ML_FF):last_surf_FF(i_surf_ML_FF),:);
            wvfm_coh_block = wvfm_coh(init_surf_FF(i_surf_ML_FF):last_surf_FF(i_surf_ML_FF),:);
        end
        %compute the ML waveform without any alignment between them (average power waveforms)
        pow_ML_FF_noalign(i_surf_ML_FF,:)=nanmean(abs(data_block).^2,1);
        if(strcmp(processing_mode,'SIN'))
            pow_ML_FF_2_noalign(i_surf_ML_FF,:)=nanmean(abs(data_block_2).^2,1);     
			pow_ML_FF_comb_noalign(i_surf_ML_FF,:)=nanmean(abs(data_block_comb).^2,1);
        end

        %move to the range-frequency domain
        data_block = fftshift(ifft(data_block,[],2),2);
        if(strcmp(processing_mode,'SIN'))
            data_block_2 = fftshift(ifft(data_block_2,[],2),2);
			data_block_comb = fftshift(ifft(data_block_comb,[],2),2);
        end
        
        % Computation of the reference window delay and altitude w.r.t
        % which perform the alignment: 0=1st SL surf; 1=central SL surf
        switch ref_alignment
            case 0
                % alignment w.r.t first SLC surface
                wd_ref            = win_delay_surf(init_surf_FF(i_surf_ML_FF));
                alt_ref           = alt_sat_surf(init_surf_FF(i_surf_ML_FF));
                
            case 1
                % alignment w.r.t central position over the different SLC surfaces to be averaged (in datation time & window delay)
                wd_ref            = mean(win_delay_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)));
                
                % to get the altitude of reference (interpolate at the central datation time)
                % computed as the average datation time over the block of
                % SLC waveforms
                time_datation_ref = mean(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)));
                
                %Interpolate at the central point using the altitude info
                %at that specific point based on time_datation_ref
                if((i_surf_ML_FF == N_total_ML_FF) && (i_ref_surf_FF == N_surfaces_FF))                   
                    alt_ref = alt_sat_surf(i_ref_surf_FF);
                else 
                    alt_ref           = spline(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),...
                                           alt_sat_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),time_datation_ref);   
                end                                     
        end

        %------------------ Compute/apply the phase correction ------------------
        for i_surf_FF=init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)
                        
            % shift due to differences between the window delays and altitude of
            % satellite for each surface w.r.t reference window delay and altitude 
            % two options foreseen: 
            % 0) the first SLC in the block "init_surf_FF(i_surf_ML_FF)"
            % 1) using the average window delay center of the collection of SLC surfaces
            % This is based on what is performed in LR-OS
            
            % JPLZ: wdshift same for both antennas when doing SARIn right? Maybe it should include modifications due to roll?
%            wdshift = ((win_delay_surf(i_surf_FF)-win_delay_surf(init_surf_FF(i_surf_ML_FF)))...
%                      -(alt_sat_surf(i_surf_FF)-alt_sat_surf(init_surf_FF(i_surf_ML_FF)))*2/c_cst)./(T0_surf(i_surf_ML_FF)/zp_fact_rg);
            wdshift(i_surf_FF-init_surf_FF(i_surf_ML_FF)+1) = ((win_delay_surf(i_surf_FF)-wd_ref)...
                     -(alt_sat_surf(i_surf_FF)-alt_ref)*2/c_cst)./(T0_surf(i_surf_FF)/zp_fact_rg);
            
            i_surf_within_ML_block = i_surf_within_ML_block+1;
            
            %apply the phase to the input data block as a phase ramp
            data_block(i_surf_within_ML_block,:) = data_block(i_surf_within_ML_block,:).*...
                exp(1i*2*pi_cst.*wdshift(i_surf_within_ML_block).*f_rg_norm);
            if(strcmp(processing_mode,'SIN'))
            data_block_2(i_surf_within_ML_block,:) = data_block_2(i_surf_within_ML_block,:).*...
                exp(1i*2*pi_cst.*wdshift(i_surf_within_ML_block).*f_rg_norm);
			data_block_comb(i_surf_within_ML_block,:) = data_block_comb(i_surf_within_ML_block,:).*...
                exp(1i*2*pi_cst.*wdshift(i_surf_within_ML_block).*f_rg_norm);
            end
            
        end% swapping all the FF non-ML surfaces of interest
        
        %-------- Move to the range domain via FFT ------------------------
        %move to power for incoherent averaging
        data_block = abs(fft(fftshift(data_block,2),[],2)).^2;
        if(strcmp(processing_mode,'SIN'))
            data_block_2 = abs(fft(fftshift(data_block_2,2),[],2)).^2;
			data_block_comb = abs(fft(fftshift(data_block_comb,2),[],2)).^2;
        end
        %--------- Average the waveforms in data block (and phase/coherence in case of SARIn mode) --------------------
        pow_ML_FF(i_surf_ML_FF,:)=nanmean(data_block,1);
        if(strcmp(processing_mode,'SIN'))
            pow_ML_FF_2(i_surf_ML_FF,:)=nanmean(data_block_2,1);
			pow_ML_FF_comb(i_surf_ML_FF,:)=nanmean(data_block_comb,1);
            phase_difference_ML_FF(i_surf_ML_FF,:)=nanmean(phase_difference_block,1);
            wvfm_coh_ML_FF(i_surf_ML_FF,:)=nanmean(wvfm_coh_block,1);      
        end
        %-------- Reapply the RMC mask ------------------------------------
        if any(processID(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF))==1)
            %compute the mask
            rmc_mask = ones(1,N_samples_rg);
            rmc_mask((N_samples_rg/zp_fact_rg/ 2 + i_sample_start_chd - 1 - rmc_margin)...
                * zp_fact_rg + 1 :N_samples_rg/zp_fact_rg * zp_fact_rg) = 0;
            
            if (i_sample_start_chd > 1)
                
                rmc_mask(1:(i_sample_start_chd - 1) * zp_fact_rg) = 0;
            end
            
            %apply the mask
            pow_ML_FF(i_surf_ML_FF,:) = pow_ML_FF(i_surf_ML_FF,:).*rmc_mask;
            if(strcmp(processing_mode,'SIN'))
                pow_ML_FF_2(i_surf_ML_FF,:) = pow_ML_FF_2(i_surf_ML_FF,:).*rmc_mask;
                pow_ML_FF_comb(i_surf_ML_FF,:) = pow_ML_FF_comb(i_surf_ML_FF,:).*rmc_mask;
            end
        end
                
        
        %---------- Set output variables for ML waveforms -----------------
        switch ref_alignment
            case 0
                %considering the first FF surface as the reference
                lat_surf_ML(i_surf_ML_FF)               = lat_surf(init_surf_FF(i_surf_ML_FF));
                lon_surf_ML(i_surf_ML_FF)               = lon_surf(init_surf_FF(i_surf_ML_FF));
                win_delay_surf_ML(i_surf_ML_FF)         = win_delay_surf(init_surf_FF(i_surf_ML_FF));
                time_surf_ML(i_surf_ML_FF)              = time_surf(init_surf_FF(i_surf_ML_FF));
                x_vel_sat_surf_ML(i_surf_ML_FF)         = x_vel_sat_surf(init_surf_FF(i_surf_ML_FF));
                y_vel_sat_surf_ML(i_surf_ML_FF)         = y_vel_sat_surf(init_surf_FF(i_surf_ML_FF));
                z_vel_sat_surf_ML(i_surf_ML_FF)         = z_vel_sat_surf(init_surf_FF(i_surf_ML_FF));
                alt_rate_sat_surf_ML(i_surf_ML_FF)      = alt_rate_sat_surf(init_surf_FF(i_surf_ML_FF));
                alt_sat_surf_ML(i_surf_ML_FF)           = alt_sat_surf(init_surf_FF(i_surf_ML_FF));
                pitch_sat_surf_ML(i_surf_ML_FF)         = pitch_sat_surf(init_surf_FF(i_surf_ML_FF));
                roll_sat_surf_ML(i_surf_ML_FF)          = roll_sat_surf(init_surf_FF(i_surf_ML_FF));
                yaw_sat_surf_ML(i_surf_ML_FF)           = yaw_sat_surf(init_surf_FF(i_surf_ML_FF));
                T0_surf_ML(i_surf_ML_FF)                = T0_surf(init_surf_FF(i_surf_ML_FF));
                pri_surf_ML(i_surf_ML_FF)               = pri_surf(init_surf_FF(i_surf_ML_FF));
                
            case 1
                %considering the central point in the block of SLC surfaces
                % Central datation time (time_datation_ref) to be used for interpolating the
                % geographical, velocities/attitude information
                if((i_surf_ML_FF == N_total_ML_FF) && (i_ref_surf_FF == N_surfaces_FF))
                    lat_surf_ML(i_surf_ML_FF)               = lat_surf(i_ref_surf_FF);
                    lon_surf_ML(i_surf_ML_FF)               = lon_surf(i_ref_surf_FF);
                    win_delay_surf_ML(i_surf_ML_FF)         = win_delay_surf(i_ref_surf_FF);
                    time_surf_ML(i_surf_ML_FF)              = time_surf(i_ref_surf_FF);
                    x_vel_sat_surf_ML(i_surf_ML_FF)         = x_vel_sat_surf(i_ref_surf_FF);
                    y_vel_sat_surf_ML(i_surf_ML_FF)         = y_vel_sat_surf(i_ref_surf_FF);
                    z_vel_sat_surf_ML(i_surf_ML_FF)         = z_vel_sat_surf(i_ref_surf_FF);
                    alt_rate_sat_surf_ML(i_surf_ML_FF)      = alt_rate_sat_surf(i_ref_surf_FF);
                    alt_sat_surf_ML(i_surf_ML_FF)           = alt_sat_surf(i_ref_surf_FF);
                    pitch_sat_surf_ML(i_surf_ML_FF)         = pitch_sat_surf(i_ref_surf_FF);
                    roll_sat_surf_ML(i_surf_ML_FF)          = roll_sat_surf(i_ref_surf_FF);
                    yaw_sat_surf_ML(i_surf_ML_FF)           = yaw_sat_surf(i_ref_surf_FF);
                    T0_surf_ML(i_surf_ML_FF)                = T0_surf(i_ref_surf_FF);
                    pri_surf_ML(i_surf_ML_FF)               = pri_surf(i_ref_surf_FF);
                    
                else
                    lat_surf_ML(i_surf_ML_FF)               = spline(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),...
                        lat_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),time_datation_ref);
                    
                    lon_surf_ML(i_surf_ML_FF)               = spline(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),...
                        lon_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),time_datation_ref);
                    
                    win_delay_surf_ML(i_surf_ML_FF)         = wd_ref;
                    time_surf_ML(i_surf_ML_FF)              = time_datation_ref;
                    
                    x_vel_sat_surf_ML(i_surf_ML_FF)         = spline(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),...
                        x_vel_sat_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),time_datation_ref);
                    
                    y_vel_sat_surf_ML(i_surf_ML_FF)         = spline(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),...
                        y_vel_sat_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),time_datation_ref);
                    
                    z_vel_sat_surf_ML(i_surf_ML_FF)         = spline(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),...
                        z_vel_sat_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),time_datation_ref);
                    
                    alt_rate_sat_surf_ML(i_surf_ML_FF)      = spline(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),...
                        alt_rate_sat_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),time_datation_ref);
                    
                    alt_sat_surf_ML(i_surf_ML_FF)           = alt_ref;
                    pitch_sat_surf_ML(i_surf_ML_FF)         = spline(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),...
                        pitch_sat_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),time_datation_ref);
                    
                    roll_sat_surf_ML(i_surf_ML_FF)          = spline(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),...
                        roll_sat_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),time_datation_ref);
                    
                    yaw_sat_surf_ML(i_surf_ML_FF)           = spline(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),...
                        yaw_sat_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),time_datation_ref);
                    
                    T0_surf_ML(i_surf_ML_FF)                = spline(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),...
                        T0_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),time_datation_ref);
                    
                    pri_surf_ML(i_surf_ML_FF)               = spline(time_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),...
                        pri_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)),time_datation_ref);
                end
                
        % Average also the sigma0_scale_factor:
        % sigma0_sc_factor_surf_ML(i_surf_ML_FF) = nanmean(sigma0_sc_factor_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)));
        % Average also the number of averaged echoes per SL waveform:
        num_aver_echoes_ML(i_surf_ML_FF) = nanmean(num_aver_echoes_surf(init_surf_FF(i_surf_ML_FF): last_surf_FF(i_surf_ML_FF)));
        
    end% total number of ML surfaces

end