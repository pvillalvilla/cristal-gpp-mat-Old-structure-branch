%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: Fully focused time domain processing from L1A to L1B
%
% Calling:
% INPUTs:
% filesBulk:  Structure of the input files including L0 and potentially L1A
% cnf:        Configuration parameters structure
% chd:        Characterization parameters structure
% cst:        Constant parameters structure
%
% OUTPUTs:
%
%
% COMMENTS:
% Current version the L0 from ARESYS is used directly, and the window
% delays of the pulses are computed within this processor
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%            Albert Garcia / isardSAT
%            Ferran Gibert / isardSAT
% ----------------------------------------------------------
% Version record:
% 1.1 2019/06/20 - Enable possibility of forcing processing over desired and
%                  coordinates (cnf.FFt.force_POI)
%                - Enable interferometry processing (swath processing +
%                  coherence)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filesBulk] = FFt_chain (filesBulk, chd, cnf, cst, options)
%% INIT VARIABLES

N_bursts            = get_num_record(filesBulk.filename_L0,'nb');
[L1A]               = create_L1A_struct(cnf, chd);
[L1A_buffer]        = create_L1A_struct(cnf, chd);
%[L1BS_buffer]        = create_L1BS_struct(cnf, chd);


%% READ THE L0 FROM ARESYS

 bursts_skipped=0;       
for i_burst = 1:N_bursts
    
     %Flag to skip burst 2 each RC 
    if cnf.skip_burst2_xRC && (strcmp(cnf.processing_mode,'SIN')) && (strcmp(chd.meas_mode,'OPEN_BURST')) && (i_burst==2||any(i_burst==2+12*int16(1:59)))
        %burst 2 in sarinob, SKIP
        i_burst=i_burst+1;
        bursts_skipped=bursts_skipped+1;
    end
    
    [L1A,filesBulk]  = read_adapt_L0_record(filesBulk,L1A,i_burst,N_bursts, cnf, chd, cst);
    % Set to zero the CAL4 pulses in pulse 10 of burst 2 each RC 
    if (cnf.CAL4_to_zero==1) && (i_burst==2||any(i_burst==2+12*int16(1:59)))
        L1A.wfm_cal_gain_corrected(10,:)=0;
        L1A.wfm_cal_gain_corrected_2(10,:)=0; 
    end
    L1A_buffer(i_burst-bursts_skipped) = L1A;
end
disp(strcat('Processing: ',filesBulk.filename_L0))


%% DATATION & WINDOW DELAY & PRI & T0 PER PULSE
[time_pulse,win_delay_pulse,pri_pulse,T0_pulse,Process_ID_pulse,...
    L1A_buffer] = ...
    datation_win_delay_pulse(L1A_buffer, cnf, chd, cst);


%% OSV & ATTITUDE SELECTION PER PULSE (MEASURED SIGNAL)
[x_sar_sat_pulse, y_sar_sat_pulse, z_sar_sat_pulse,...
    x_vel_sat_sar_pulse,  y_vel_sat_sar_pulse,  z_vel_sat_sar_pulse,...
    lat_sar_sat_pulse, lon_sar_sat_pulse,...
    alt_sar_sat_pulse,alt_rate_sar_sat_pulse,...
    roll_sar_pulse, pitch_sar_pulse, yaw_sar_pulse] = ...
    osv_attitude_pulse(L1A_buffer, time_pulse, cnf, chd, cst);


%% WAVEFORMS RE-ORGANIZATION
% Concatenate all the waveforms of the pulses from the different bursts
%%% On-board reversion for RMC assumed to be performed in L1A
N_total_pulses = length(time_pulse);
N_samples_rg   = length([L1A_buffer(1).wfm_cal_gain_corrected(1,:)]);
wfm_cal_gain_corrected = zeros(N_total_pulses,N_samples_rg);

%if(strcmp(cnf.processing_mode,'SIN')) % JPLZ: we define SARin variables even if doing just SAR for the parfor loop not to yield errors
    wfm_cal_gain_corrected_2 = zeros(N_total_pulses,N_samples_rg);
%end

for i_burst = 1:N_bursts-bursts_skipped
    
    wfm_cal_gain_corrected((i_burst-1)*chd.N_pulses_burst+1:(i_burst-1)*chd.N_pulses_burst+chd.N_pulses_burst,:)=...
        L1A_buffer(i_burst).wfm_cal_gain_corrected;
    
    if(strcmp(cnf.processing_mode,'SIN'))
        wfm_cal_gain_corrected_2((i_burst-1)*chd.N_pulses_burst+1:(i_burst-1)*chd.N_pulses_burst+chd.N_pulses_burst,:)=...
            L1A_buffer(i_burst).wfm_cal_gain_corrected_2;
    end
end

% TESTING JPL: In SARIN OB try not to take the CAL pulses
% for i_burst = 1:N_bursts
%     
%     if (strcmp(chd.meas_mode,'OPEN_BURST')) &&  (i_burst==1||any(i_burst==1+12*int16(1:59)))
%         
%     wfm_cal_gain_corrected((i_burst-1)*chd.N_pulses_burst+1:(i_burst-1)*chd.N_pulses_burst+chd.N_pulses_burst,:)=...
%         L1A_buffer(i_burst).wfm_cal_gain_corrected(5:64,:);
%     
%     else
%         wfm_cal_gain_corrected((i_burst-1)*chd.N_pulses_burst+1:(i_burst-1)*chd.N_pulses_burst+chd.N_pulses_burst,:)=...
%          L1A_buffer(i_burst).wfm_cal_gain_corrected;
%         
%     end
%     
%     if(strcmp(cnf.processing_mode,'SIN'))
%         wfm_cal_gain_corrected_2((i_burst-1)*chd.N_pulses_burst+1:(i_burst-1)*chd.N_pulses_burst+chd.N_pulses_burst,:)=...
%             L1A_buffer(i_burst).wfm_cal_gain_corrected_2;
%     end
% end


%% GRID (SURFACE) COMPUTATION
[x_surf, y_surf, z_surf,...
    lat_surf, lon_surf, alt_surf, win_delay_surf,...
    time_surf,...
    x_vel_sat_surf,y_vel_sat_surf,z_vel_sat_surf,...
    x_sat_surf,y_sat_surf,z_sat_surf,...
    alt_rate_sat_surf,alt_sat_surf,...
    pitch_sat_surf,roll_sat_surf,yaw_sat_surf,...
    TRP_surf_idx] = ...
    focusing_grid_computation(time_pulse,...
    alt_rate_sar_sat_pulse, ...
    x_vel_sat_sar_pulse, y_vel_sat_sar_pulse, z_vel_sat_sar_pulse,...
    roll_sar_pulse, pitch_sar_pulse, yaw_sar_pulse,...
    x_sar_sat_pulse, y_sar_sat_pulse, z_sar_sat_pulse,...
    lat_sar_sat_pulse,lon_sar_sat_pulse,alt_sar_sat_pulse,...
    win_delay_pulse,pri_pulse,T0_pulse,...
    cnf,chd,cst);

%% Forcing number of surfaces for testing

    if (cnf.FFt.force_max_N_surf)
        if cnf.FFt.force_max_N_surf
            N_surfaces                          = cnf.FFt.max_N_surf; 
        else
            N_surfaces                          = length(x_surf);
        end
    else
        N_surfaces                          = length(x_surf);
    end

%% Variables definition for paralelization

% write Intermediate product
%    [filesBulk] = create_NetCDF_L1ASS(filesBulk, N_bursts, N_surfaces, cnf, chd, cst);
%
%    for i_burst = 1:N_bursts
%         write_NetCDF_L1ASS(filesBulk, L1A_buffer(i_burst),i_burst, cnf,chd,cst);
%    end
%    write_NetCDF_L1ASS_surf(filesBulk, time_surf, x_surf, y_surf, z_surf, cnf,chd,cst);

%tau_delay_pixel_time_apert = zeros(N_surfaces,N_samples_az);

wfm_AC = zeros(N_surfaces,N_samples_rg.*cnf.zp_fact_range);

%if(strcmp(cnf.processing_mode,'SIN')) % JPLZ: we define SARin variables even if doing just SAR for the parfor loop not to yield errors
    wfm_AC_2 = zeros(N_surfaces,N_samples_rg.*cnf.zp_fact_range);
    wfm_AC_combined = zeros(N_surfaces,N_samples_rg.*cnf.zp_fact_range);
    
    cross_product     = zeros(N_surfaces, N_samples_rg.*cnf.zp_fact_range);
    phase_difference  = zeros(N_surfaces, N_samples_rg.*cnf.zp_fact_range);
    
    wfm_AC_avPower_1  = zeros(N_surfaces, N_samples_rg.*cnf.zp_fact_range);
    wfm_AC_avPower_2  = zeros(N_surfaces, N_samples_rg.*cnf.zp_fact_range);
    
    wfm_avPower_SARin = zeros(N_surfaces, N_samples_rg.*cnf.zp_fact_range);
    
    wvfm_coh          = zeros(N_surfaces, N_samples_rg.*cnf.zp_fact_range);
    
%end

%matlab version for paralelization issues
version_matlab=version;


%% FOCUSING PROCESS

if cnf.FFt.num_pools~=1
    if str2double(version_matlab(end-5:end-2))>2013
        parpool(cnf.FFt.num_pools);
    else
        matlabpool('open',cnf.FFt.num_pools);
    end
    parfor_progress(N_surfaces);
    parfor i_surf = 1: N_surfaces 
         %disp(strcat(cnf.processing_mode));
         disp(strcat('Surf: ',num2str(i_surf)));
        % Select indexes where to compute backprojection according to
        % integration time
        indexes_integration = logical(ones(1, N_total_pulses));
        if cnf.FFt.T_integration ~= -1
            indexes_integration = (time_pulse>=(time_surf(i_surf)-cnf.FFt.T_integration/2) & time_pulse<=(time_surf(i_surf)+cnf.FFt.T_integration/2));
        end
        
        % Compute backprojection on channel rx1:
        [wfm_AC(i_surf,:), wfm_RC_phase_corr_masked, N_averages_SL(i_surf)] = ...
            backprojection_kernel_freq_domain(...
            wfm_cal_gain_corrected,...
            T0_pulse,...
            win_delay_pulse,...
            x_sar_sat_pulse,y_sar_sat_pulse,z_sar_sat_pulse,...
            lat_sar_sat_pulse,lon_sar_sat_pulse,alt_sar_sat_pulse,...
            pitch_sar_pulse,roll_sar_pulse,yaw_sar_pulse,...
            x_vel_sat_sar_pulse,y_vel_sat_sar_pulse,z_vel_sat_sar_pulse,...
            x_surf(i_surf),y_surf(i_surf),z_surf(i_surf),win_delay_surf(i_surf),...
            N_total_pulses,N_samples_rg,...
            indexes_integration, ...
            cnf,cst,chd,...
            'ProcessID_pulse',Process_ID_pulse,...
            'i_sample_start_chd',chd.RMC_start_sample);
        
        % Eventually compute backprojection on channel rx2:
        if(strcmp(cnf.processing_mode,'SIN'))
            [wfm_AC_2(i_surf,:), wfm_RC_phase_corr_masked_2, N_averages_SL_2(i_surf)] = ...
                backprojection_kernel_freq_domain(...
                wfm_cal_gain_corrected_2,...
                T0_pulse,...
                win_delay_pulse, ...
                x_sar_sat_pulse,y_sar_sat_pulse,z_sar_sat_pulse,...
                lat_sar_sat_pulse,lon_sar_sat_pulse,alt_sar_sat_pulse,...
                pitch_sar_pulse,roll_sar_pulse,yaw_sar_pulse,...
                x_vel_sat_sar_pulse,y_vel_sat_sar_pulse,z_vel_sat_sar_pulse,...
                x_surf(i_surf),y_surf(i_surf),z_surf(i_surf),win_delay_surf(i_surf),...
                N_total_pulses,N_samples_rg,...
                indexes_integration, ...
                cnf,cst,chd,...
                'ProcessID_pulse',Process_ID_pulse,...
                'i_sample_start_chd',chd.RMC_start_sample);
            
            % Evaluate cross product
            cross_product(i_surf,:) = sum(wfm_RC_phase_corr_masked.*conj(wfm_RC_phase_corr_masked_2), 1)/numel(wfm_RC_phase_corr_masked(:,1));
			
			% Compute the combined waveform
            wfm_AC_combined(i_surf,:) =  nanmean(cross_product(i_surf,:).^2,1); 

            % Evaluate phase difference
            phase_difference(i_surf,:) = angle(cross_product(i_surf,:)); %recoverying phase between [-pi,pi] (wraped)
            
            % Average waveform power per surface
            wfm_AC_avPower_1(i_surf,:) = sum(wfm_RC_phase_corr_masked.*conj(wfm_RC_phase_corr_masked), 1)/numel(wfm_RC_phase_corr_masked(:,1));
            wfm_AC_avPower_2(i_surf,:) = sum(wfm_RC_phase_corr_masked_2.*conj(wfm_RC_phase_corr_masked_2), 1)/numel(wfm_RC_phase_corr_masked_2(:,1));
            
            % Compute SARin power
            wfm_avPower_SARin(i_surf,:) = sum(0.5*(wfm_RC_phase_corr_masked.*conj(wfm_RC_phase_corr_masked) + ...
                wfm_RC_phase_corr_masked_2.*conj(wfm_RC_phase_corr_masked_2)))/ ...
                numel(wfm_RC_phase_corr_masked(:,1));
            
            % Evaluate coherence
            wvfm_coh(i_surf,:) = sqrt(cross_product(i_surf,:).*conj(cross_product(i_surf,:))./(wfm_AC_avPower_1(i_surf,:).*wfm_AC_avPower_2(i_surf,:)));
            
        end
        
        % Pack data to write file:
        %L1B(i_surf).beam_ang_surf
        %L1B(i_surf).burst_index

        L1B(i_surf).T0_sar_surf    = T0_pulse(1);
        L1B(i_surf).time_surf      = time_surf(i_surf);

        L1B(i_surf).x_sat          = x_sat_surf(i_surf);
        L1B(i_surf).y_sat          = y_sat_surf(i_surf);
        L1B(i_surf).z_sat          = z_sat_surf(i_surf);

        L1B(i_surf).lat_sat        = lat_surf(i_surf);
        L1B(i_surf).lon_sat        = lon_surf(i_surf);
        L1B(i_surf).alt_sat        = alt_sat_surf(i_surf);
        L1B(i_surf).alt_rate_sat   = alt_rate_sat_surf(i_surf);

        L1B(i_surf).pitch_surf     = pitch_sat_surf(i_surf);
        L1B(i_surf).roll_surf      = roll_sat_surf(i_surf);
        L1B(i_surf).yaw_surf       = yaw_sat_surf(i_surf);

        L1B(i_surf).x_vel_sat      = x_vel_sat_surf(i_surf);
        L1B(i_surf).y_vel_sat      = y_vel_sat_surf(i_surf);
        L1B(i_surf).z_vel_sat      = z_vel_sat_surf(i_surf);

        L1B(i_surf).win_delay_surf = win_delay_surf(i_surf);

        L1B(i_surf).wfm_scaling_factor     = 1; % Set to one provionally

        L1B(i_surf).wfm_cor_i2q2        = wfm_AC(i_surf,:);
        
        if(strcmp(cnf.processing_mode,'SIN')) 
            L1B(i_surf).wfm_cor_i2q2_2      = wfm_AC_2(i_surf,:);
			%L1B(i_surf).wfm_cor_i2q2_comb   = wfm_AC_combined(i_surf,:);
            L1B(i_surf).phase_difference    = phase_difference(i_surf,:);
            L1B(i_surf).coherence           = wvfm_coh(i_surf,:);
        end
        parfor_progress;
        
    end %end different surfaces
    
    parfor_progress(0);
    %close pools
    if str2double(version_matlab(end-5:end-2))>2013
        poolobj = gcp('nocreate');
        delete(poolobj);
    else
        matlabpool('close');
    end
else
    for i_surf = 1: N_surfaces
        disp(strcat('Surf: ',num2str(i_surf)));
        
        % Select indexes where to compute backprojection according to
        % integration time
        indexes_integration = logical(ones(1, N_total_pulses));
        if cnf.FFt.T_integration ~= -1
            indexes_integration = (time_pulse>=(time_surf(i_surf)-cnf.FFt.T_integration/2) & time_pulse<=(time_surf(i_surf)+cnf.FFt.T_integration/2));
        end
        
        % Compute backprojection on channel rx1:
        [wfm_AC(i_surf,:), wfm_RC_phase_corr_masked, N_averages_SL(i_surf)] = ...
            backprojection_kernel_freq_domain(...
            wfm_cal_gain_corrected,...
            T0_pulse,...
            win_delay_pulse,...
            x_sar_sat_pulse,y_sar_sat_pulse,z_sar_sat_pulse,...
            lat_sar_sat_pulse,lon_sar_sat_pulse,alt_sar_sat_pulse,...
            pitch_sar_pulse,roll_sar_pulse,yaw_sar_pulse,...
            x_vel_sat_sar_pulse,y_vel_sat_sar_pulse,z_vel_sat_sar_pulse,...
            x_surf(i_surf),y_surf(i_surf),z_surf(i_surf),win_delay_surf(i_surf),...
            N_total_pulses,N_samples_rg,...
            indexes_integration, ...
            cnf,cst,chd,...
            'ProcessID_pulse',Process_ID_pulse,...
            'i_sample_start_chd',chd.RMC_start_sample);
        
        % Eventually compute backprojection on channel rx2:
        if(strcmp(cnf.processing_mode,'SIN'))
            [wfm_AC_2(i_surf,:), wfm_RC_phase_corr_masked_2, N_averages_SL_2(i_surf)] = ...
                backprojection_kernel_freq_domain(...
                wfm_cal_gain_corrected_2,...
                T0_pulse,...
                win_delay_pulse, ...
                x_sar_sat_pulse,y_sar_sat_pulse,z_sar_sat_pulse,...
                lat_sar_sat_pulse,lon_sar_sat_pulse,alt_sar_sat_pulse,...
                pitch_sar_pulse,roll_sar_pulse,yaw_sar_pulse,...
                x_vel_sat_sar_pulse,y_vel_sat_sar_pulse,z_vel_sat_sar_pulse,...
                x_surf(i_surf),y_surf(i_surf),z_surf(i_surf),win_delay_surf(i_surf),...
                N_total_pulses,N_samples_rg,...
                indexes_integration, ...
                cnf,cst,chd,...
                'ProcessID_pulse',Process_ID_pulse,...
                'i_sample_start_chd',chd.RMC_start_sample);
            
            % Evaluate cross product
            cross_product(i_surf,:) = sum(wfm_RC_phase_corr_masked.*conj(wfm_RC_phase_corr_masked_2), 1)/numel(wfm_RC_phase_corr_masked(:,1));
            
			% Compute the combined waveform
            wfm_AC_combined(i_surf,:) =  nanmean(cross_product(i_surf,:).^2,1); 
			
            % Evaluate phase difference
            phase_difference(i_surf,:) = angle(cross_product(i_surf,:)); %recoverying phase between [-pi,pi] (wraped)
            
            % Average waveform power per surface
            wfm_AC_avPower_1(i_surf,:) = sum(wfm_RC_phase_corr_masked.*conj(wfm_RC_phase_corr_masked), 1)/numel(wfm_RC_phase_corr_masked(:,1));
            wfm_AC_avPower_2(i_surf,:) = sum(wfm_RC_phase_corr_masked_2.*conj(wfm_RC_phase_corr_masked_2), 1)/numel(wfm_RC_phase_corr_masked_2(:,1));
            
            % Compute SARin power
            wfm_avPower_SARin(i_surf,:) = sum(0.5*(wfm_RC_phase_corr_masked.*conj(wfm_RC_phase_corr_masked) + ...
                wfm_RC_phase_corr_masked_2.*conj(wfm_RC_phase_corr_masked_2)))/ ...
                numel(wfm_RC_phase_corr_masked(:,1));
            
            % Evaluate coherence
            wvfm_coh(i_surf,:) = sqrt(cross_product(i_surf,:).*conj(cross_product(i_surf,:))./(wfm_AC_avPower_1(i_surf,:).*wfm_AC_avPower_2(i_surf,:)));
            
            
        end
        
        % Pack data to write file:
        %L1B(i_surf).beam_ang_surf
        %L1B(i_surf).burst_index
        
        L1B(i_surf).T0_sar_surf    = T0_pulse(1);
        L1B(i_surf).time_surf      = time_surf(i_surf);
        
        L1B(i_surf).x_sat          = x_sat_surf(i_surf);
        L1B(i_surf).y_sat          = y_sat_surf(i_surf);
        L1B(i_surf).z_sat          = z_sat_surf(i_surf);
        
        L1B(i_surf).lat_sat        = lat_surf(i_surf);
        L1B(i_surf).lon_sat        = lon_surf(i_surf);
        L1B(i_surf).alt_sat        = alt_sat_surf(i_surf);
        L1B(i_surf).alt_rate_sat   = alt_rate_sat_surf(i_surf);
        
        L1B(i_surf).pitch_surf     = pitch_sat_surf(i_surf);
        L1B(i_surf).roll_surf      = roll_sat_surf(i_surf);
        L1B(i_surf).yaw_surf       = yaw_sat_surf(i_surf);
        
        L1B(i_surf).x_vel_sat      = x_vel_sat_surf(i_surf);
        L1B(i_surf).y_vel_sat      = y_vel_sat_surf(i_surf);
        L1B(i_surf).z_vel_sat      = z_vel_sat_surf(i_surf);
        
        L1B(i_surf).win_delay_surf = win_delay_surf(i_surf);
        
        L1B(i_surf).wfm_scaling_factor     = 1; % Set to one provionally
        
        L1B(i_surf).wfm_cor_i2q2        = wfm_AC(i_surf,:);
        if(strcmp(cnf.processing_mode,'SIN')) 
            L1B(i_surf).wfm_cor_i2q2_2      = wfm_AC_2(i_surf,:);
%			L1B(i_surf).wfm_cor_i2q2_comb   = wfm_AC_combined(i_surf,:);
            L1B(i_surf).phase_difference    = phase_difference(i_surf,:);
            L1B(i_surf).coherence           = wvfm_coh(i_surf,:);
        end
        
    end %end different surfaces
    
end

%scaled waveforms considering all azimuth values?
% max_wfm_AC = max(max(abs(wfm_AC)));
% wfm_sc = max_wfm_AC / (2^16-2);
% wfm_AC_sc=(round(abs(wfm_AC)./ wfm_sc));
% 
% max_wfm_AC_2 = max(max(abs(wfm_AC_2)));
% wfm_sc_2 = max_wfm_AC_2 / (2^16-2);
% wfm_AC_sc_2=(round(abs(wfm_AC_2)./ wfm_sc_2));


% % Plot results in case of SARIN
% if(strcmp(cnf.processing_mode,'SIN'))
%     
%     % Plot coherence between channels
%     figure;
%     imagesc(wvfm_coh);
%     colorbar
%     title('Coherence')
%     
%     % Plot phase difference
%     figure;
%     imagesc(phase_difference);
%     title('Phase difference')
%     colorbar
%     
% end

%     % Sigma0 Scaling Factor
%     disp(['Processing sigma0 scaling factor...']);    
%     [sigma0_scaling_factor_surf] = sigma0_scal_fact_FF(...
%             N_surfaces, N_samples_az, ...
%             x_surf, y_surf, z_surf, ...
%             x_sar_sat_pulse, y_sar_sat_pulse, z_sar_sat_pulse, ...
%             x_vel_sat_sar_pulse, y_vel_sat_sar_pulse, z_vel_sat_sar_pulse, ...
%             time_sar_ku_pulse, time_surf, ...
%             pri_sar_pre_dat_pulse, 1, 1);

 %% 11. MULTILOOKING (INCOHERENT AVERAGING)
 if ~cnf.trp_flag
     ProcessID=ones(1,N_surfaces);     
     disp(['Multilooking...']);   
     if(strcmp(cnf.processing_mode,'SIN'))
            [pow_ML_FF, ~, pow_ML_FF_2, ~, pow_ML_FF_comb,~, phase_difference_ML, wvfm_coh_ML,... 
             lat_surf_ML, lon_surf_ML, win_delay_surf_ML, time_surf_ML, ...
             x_vel_sat_surf_ML,y_vel_sat_surf_ML,z_vel_sat_surf_ML,...
             alt_rate_sat_surf_ML,alt_sat_surf_ML,...
             pitch_sat_surf_ML,roll_sat_surf_ML,yaw_sat_surf_ML,...
             num_aver_echoes_ML, ...
             T0_surf_ML,pri_surf_ML,N_looks,looks_flag]...
             ...
             = multilooking_FF ...
             (wfm_AC, wfm_AC_2, wfm_AC_combined,phase_difference, wvfm_coh, ...
             lat_surf, lon_surf, win_delay_surf, time_surf,...
             x_vel_sat_surf, y_vel_sat_surf, z_vel_sat_surf,...
             alt_rate_sat_surf, alt_sat_surf, ...
             pitch_sat_surf, roll_sat_surf, yaw_sat_surf,...
             N_averages_SL, ...
             T0_pulse, pri_pulse,...
             cnf.FFt.out_spatial_res,...
             cnf.zp_fact_range,...
             cst.semi_major_axis,cst.semi_minor_axis,cst.c,cst.pi,...
             ProcessID,chd.i_sample_start,cnf.rmc_margin,cnf.processing_mode, cnf.FFt.ref_alignment_ML);
     
         for i_surf=1:size(pow_ML_FF,1)
             % Insert variables in struct
             L1B(i_surf).pow_ML_FF = pow_ML_FF(i_surf,:);
             L1B(i_surf).pow_ML_FF_2 = pow_ML_FF_2(i_surf,:);
			 L1B(i_surf).pow_ML_FF_comb = pow_ML_FF_comb(i_surf,:);
             L1B(i_surf).phase_difference_ML = phase_difference_ML(i_surf,:);
             L1B(i_surf).wvfm_coh_ML = wvfm_coh_ML(i_surf,:);
             L1B(i_surf).lat_surf_ML = lat_surf_ML(i_surf);
             L1B(i_surf).lon_surf_ML = lon_surf_ML(i_surf);
             L1B(i_surf).win_delay_surf_ML = win_delay_surf_ML(i_surf);
             L1B(i_surf).time_surf_ML = time_surf_ML(i_surf);
             L1B(i_surf).x_vel_sat_surf_ML = x_vel_sat_surf_ML(i_surf);
             L1B(i_surf).y_vel_sat_surf_ML = y_vel_sat_surf_ML(i_surf);
             L1B(i_surf).z_vel_sat_surf_ML = z_vel_sat_surf_ML(i_surf);
             L1B(i_surf).alt_rate_sat_surf_ML = alt_rate_sat_surf_ML(i_surf);
             L1B(i_surf).alt_sat_surf_ML = alt_sat_surf_ML(i_surf);
             L1B(i_surf).T0_surf_ML = T0_surf_ML(i_surf);
             L1B(i_surf).pri_surf_ML = pri_surf_ML(i_surf);
             L1B(i_surf).N_looks = N_looks(i_surf);
             L1B(i_surf).pitch_sat_surf_ML = pitch_sat_surf_ML(i_surf);
             L1B(i_surf).roll_sat_surf_ML = roll_sat_surf_ML(i_surf);
             L1B(i_surf).yaw_sat_surf_ML = yaw_sat_surf_ML(i_surf);
             L1B(i_surf).wfm_scaling_factor_ML     = 1; % Set to one provionally
         end
     
     elseif(strcmp(cnf.processing_mode,'SAR'))
            [pow_ML_FF, ~, ~, ~, ~, ~,~,~,... 
             lat_surf_ML, lon_surf_ML, win_delay_surf_ML, time_surf_ML, ...
             x_vel_sat_surf_ML,y_vel_sat_surf_ML,z_vel_sat_surf_ML,...
             alt_rate_sat_surf_ML,alt_sat_surf_ML,...
             pitch_sat_surf_ML,roll_sat_surf_ML,yaw_sat_surf_ML,...
             num_aver_echoes_ML, ...
             T0_surf_ML,pri_surf_ML,N_looks,looks_flag]...
             ...
             = multilooking_FF ...
             (wfm_AC, [], [], [],[], ...
             lat_surf, lon_surf, win_delay_surf, time_surf,...
             x_vel_sat_surf, y_vel_sat_surf, z_vel_sat_surf,...
             alt_rate_sat_surf, alt_sat_surf, ...
             pitch_sat_surf, roll_sat_surf, yaw_sat_surf,...
             N_averages_SL, ...
             T0_pulse, pri_pulse,...
             cnf.FFt.out_spatial_res,...
             cnf.zp_fact_range,...
             cst.semi_major_axis,cst.semi_minor_axis,cst.c,cst.pi,...
             ProcessID,chd.i_sample_start,cnf.rmc_margin,cnf.processing_mode,cnf.FFt.ref_alignment_ML);
         
         for i_surf=1:size(pow_ML_FF,1)
             % Insert variables in struct
             L1B(i_surf).pow_ML_FF = pow_ML_FF(i_surf,:);
             L1B(i_surf).lat_surf_ML = lat_surf_ML(i_surf);
             L1B(i_surf).lon_surf_ML = lon_surf_ML(i_surf);
             L1B(i_surf).win_delay_surf_ML = win_delay_surf_ML(i_surf);
             L1B(i_surf).time_surf_ML = time_surf_ML(i_surf);
             L1B(i_surf).x_vel_sat_surf_ML = x_vel_sat_surf_ML(i_surf);
             L1B(i_surf).y_vel_sat_surf_ML = y_vel_sat_surf_ML(i_surf);
             L1B(i_surf).z_vel_sat_surf_ML = z_vel_sat_surf_ML(i_surf);
             L1B(i_surf).alt_rate_sat_surf_ML = alt_rate_sat_surf_ML(i_surf);
             L1B(i_surf).alt_sat_surf_ML = alt_sat_surf_ML(i_surf);
             L1B(i_surf).T0_surf_ML = T0_surf_ML(i_surf);
             L1B(i_surf).pri_surf_ML = pri_surf_ML(i_surf);
             L1B(i_surf).N_looks = N_looks(i_surf);
             L1B(i_surf).pitch_sat_surf_ML = pitch_sat_surf_ML(i_surf);
             L1B(i_surf).roll_sat_surf_ML = roll_sat_surf_ML(i_surf);
             L1B(i_surf).yaw_sat_surf_ML = yaw_sat_surf_ML(i_surf);
             L1B(i_surf).wfm_scaling_factor_ML     = 1; % Set to one provionally
         end
         
     end
 end
 
%% Writting routines

if(cnf.writting_flag(4))
    % Create Single-Looked data file
    [filesBulk] = create_NetCDF_FFL1B(filesBulk,N_surfaces, 'FF-SL', cnf, chd, cst,0);
    type_FF = 'FF-SL';
    
    if ~cnf.trp_flag
        % Create Multi-Looked data file
        [filesBulk] = create_NetCDF_FFL1B(filesBulk, size(pow_ML_FF,1), 'FF-ML', cnf, chd, cst,0);
        type_FF = 'FF-ML';
    end

% Write Single-Looked data file
for i_surf = 1:N_surfaces
    %     if(cnf.writting_flag(2))
    %         if i_surf==1
    %             filesBulk.ncid = netcdf.open(filesBulk.filename_L1Bs,'WRITE');
    %         end
    %     end
    
    if(cnf.writting_flag(4)) % Index 4 corresponds to L1B_VHR_time
        type_FF = 'FF-SL';
        if i_surf==1
            filesBulk.ncid = netcdf.open(filesBulk.filename_L1BFF_SL, 'WRITE');
        end
        write_NetCDF_FFL1B(filesBulk,L1A_buffer, L1B(i_surf), i_surf-1, type_FF, cnf, chd, cst);
    end
    
    
    
    % %%
    %         indexes_to_remove = find(burst_indexes_in_memory(1:end) < min(L1B(i_surf).burst_index));
    %
    % %         if(cnf.verify_L1BS_internal)
    % %             plots_L1BS;
    % %         end
    %
    %         if(~isempty(indexes_to_remove))
    %             [L1A_buffer(burst_indexes_in_memory(indexes_to_remove))] = empty_L1A_struct(L1A_buffer(burst_indexes_in_memory(indexes_to_remove)));
    %             burst_indexes_in_memory(indexes_to_remove)=[];
    %
    %         end
    %
    %         [L1B(i_surf)] = empty_L1BS_struct(L1B(i_surf));
    
    %     if(L1BS_buffer(i_surf).surface_type_flag == 4)
    %         %TRP location being processed
    %         disp(['Computing the results for the record ' num2str(i_surf)]);
    %     end
     
end

% Write Multi-Looked data file
if ~cnf.trp_flag
    for i_surf = 1:size(pow_ML_FF,1)
        %     if(cnf.writting_flag(2))
        %         if i_surf==1
        %             filesBulk.ncid = netcdf.open(filesBulk.filename_L1Bs,'WRITE');
        %         end
        %     end
        
        if(cnf.writting_flag(4)) % Index 4 corresponds to L1B_VHR_time
            type_FF = 'FF-ML';
            if i_surf==1
                filesBulk.ncid = netcdf.open(filesBulk.filename_L1BFF_ML, 'WRITE');
            end
            write_NetCDF_FFL1B(filesBulk,L1A_buffer, L1B(i_surf), i_surf-1, type_FF, cnf, chd, cst);
        end
    end
end

end

if(cnf.writting_flag(4))
    netcdf.close(filesBulk.ncid);
end

%JPLZ: Provisionally save the ML data in a .mat. We have to implement
%saving different .nc files (1 for SL, 1 for ML)
% if ~cnf.trp_flag
%     L0_date_find = strfind(filesBulk.filename_L0,'_202');
%     Band=strsplit(filesBulk.filename_L0,'PIC_SIRS___');
%     Band = char(Band(2)); Band=Band(1:2);
%     filesBulk.filename_L1B_FFt = strcat(filesBulk.outputPath,'PICE_1B_FF-ML_variables_',Band,'_',...
%         filesBulk.filename_L0((L0_date_find(1)+1):(L0_date_find(1)+30)),...
%         '_isd','.mat');
%     if cnf.writting_flag(4)
%         if(strcmp(cnf.processing_mode,'SIN'))
%             save(strcat(filesBulk.filename_L1B_FFt,'.mat'),...
%                 'pow_ML_FF','pow_ML_FF_noalign','pow_ML_FF_2','pow_ML_FF_noalign_2','lat_surf_ML','lon_surf_ML','win_delay_surf_ML',...
%                 'time_surf_ML',...
%                 'x_vel_sat_surf_ML','y_vel_sat_surf_ML','z_vel_sat_surf_ML',...
%                 'alt_rate_sat_surf_ML','alt_sat_surf_ML',...
%                 'pitch_sat_surf_ML','roll_sat_surf_ML','yaw_sat_surf_ML',...
%                 'T0_surf_ML','pri_surf_ML','N_looks')
%         else
%             save(strcat(filesBulk.filename_L1B_FFt,'.mat'),...
%                 'pow_ML_FF','pow_ML_FF_noalign','lat_surf_ML','lon_surf_ML','win_delay_surf_ML',...
%                 'time_surf_ML',...
%                 'x_vel_sat_surf_ML','y_vel_sat_surf_ML','z_vel_sat_surf_ML',...
%                 'alt_rate_sat_surf_ML','alt_sat_surf_ML',...
%                 'pitch_sat_surf_ML','roll_sat_surf_ML','yaw_sat_surf_ML',...
%                 'T0_surf_ML','pri_surf_ML','N_looks')
%         end
%         
%     end
% end
% %SAVING RESULTS IN .mat
% L0_date_find = strfind(filesBulk.filename_L0,'_202');
% Band=strsplit(filesBulk.filename_L0,'PIC_SIRS___');
% Band = char(Band(2)); Band=Band(1:2);
% filesBulk.filename_L1B_FFt = strcat(filesBulk.outputPath,'PICE_1B_VHRt_',Band,'_',...
%     filesBulk.filename_L0((L0_date_find(1)+1):(L0_date_find(1)+30)),...
%     '_isd','.mat');
% if cnf.writting_flag(4)
%     if(strcmp(cnf.processing_mode,'SIN'))
%         save(strcat(filesBulk.filename_L1B_FFt,'.mat'),...
%             'wfm_AC','wfm_AC_2','wfm_AC_avPower_1','wfm_AC_avPower_2','wfm_avPower_SARin',...
%             'wfm_AC_sc','wfm_AC_sc_2',...
%             'wfm_RC_phase_corr_masked','wfm_RC_phase_corr_masked_2',...
%             'phase_difference','wvfm_coh',...
%             'wfm_cal_gain_corrected','wfm_cal_gain_corrected_2',...
%             'lat_surf', 'lon_surf', 'alt_surf', 'win_delay_surf',...
%             'time_surf',...
%             'x_vel_sat_surf','y_vel_sat_surf','z_vel_sat_surf',...
%             'alt_rate_sat_surf','alt_sat_surf',...
%             'pitch_sat_surf','roll_sat_surf','yaw_sat_surf',...
%             '-v7.3')
%     else
%         save(strcat(filesBulk.filename_L1B_FFt,'.mat'),...
%             'wfm_AC',...
%             'wfm_RC_phase_corr_masked',...
%             'wfm_cal_gain_corrected',...
%             'lat_surf', 'lon_surf', 'alt_surf', 'win_delay_surf',...
%             'time_surf',...
%             'x_vel_sat_surf','y_vel_sat_surf','z_vel_sat_surf',...
%             'alt_rate_sat_surf','alt_sat_surf',...
%             'pitch_sat_surf','roll_sat_surf','yaw_sat_surf',...
%             '-v7.3')
%     end
% end
% figure; mesh(((wfm_sc) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(abs((wfm_AC'))));
% figure; mesh(((wfm_sc) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(abs(10.*log10(wfm_AC'))));
% waveform_scale_factor_l1b_echo(100) = single(max(L1B(100).wfm_cor_i2q2)/ (2^16-2));

%% TRP analysis
if cnf.trp_flag
    tic
    font_size = 20;
    % Path to save different plots
    plotsPath = strcat(filesBulk.outputPath,'plots/');
    mkdir(plotsPath);
    %[~,filename_L1B_FFt]=fileparts(filesBulk.filename_L1B_FFt);
    [~,filename_L1B_FFt]=fileparts(filesBulk.filename_L1BFF_SL); %JPLZ: removed a t
    
    % 2D Re-sampling of the complex matrix
    wfm_AC_interp=interpolsinc_2D(wfm_AC,cnf.FFt.zp_ac_TRP,cnf.FFt.zp_al_TRP);
    %wfm_AC_interp=interpolsinc_2D(wfm_AC_2,cnf.FFt.zp_ac_TRP,cnf.FFt.zp_al_TRP); % For channel 2 in SARIn case

    % Compute distances along and slant range w.r.t TRP location
    % along
    [arclen,~]=distance(lat_surf,lon_surf,ones(1,length(lat_surf))*chd.lat_trp,ones(1,length(lat_surf))*chd.lon_trp,[cst.semi_major_axis sqrt(1-(cst.semi_minor_axis/cst.semi_major_axis).^2)]);
%     Eccentricity = sqrt(2*cst.flat_coeff-cst.flat_coeff^2);
%     [arclen,~]=distance(lat_surf,lon_surf,ones(1,length(lat_surf))*chd.lat_trp,ones(1,length(lat_surf))*chd.lon_trp,[cst.semi_major_axis Eccentricity]);
    [~,min_dist]=min(arclen);
    distances_over_arclengths=arclen;
    distances_over_arclengths(1:min_dist-1)=-1.0*distances_over_arclengths(1:min_dist-1);
    distances_over_arclengths_interp=interp(distances_over_arclengths,cnf.FFt.zp_al_TRP);
    % across
    slant_range_wrt_TRP_range_interp = ((1:1:chd.N_samples_sar*cnf.FFt.zp_ac_TRP)-((chd.N_samples_sar/2)*cnf.FFt.zp_ac_TRP+1))*cst.c/2*chd.T0_nom/cnf.FFt.zp_ac_TRP+...
        alt_sat_surf(TRP_surf_idx)-(win_delay_surf(TRP_surf_idx)*cst.c/2)-chd.alt_trp;
    
    % 2D image of the PTR
    max_data  = max(20*log10(abs(wfm_AC_interp(:))));
    max_image = 0.0;
    min_image = max_image-60.0;
    
    min_max_along = [-cnf.FFt.surf_dist_trp,cnf.FFt.surf_dist_trp];
    min_max_across = [-cnf.FFt.surf_dist_trp,cnf.FFt.surf_dist_trp];
    toc
    figure;
    imagesc(slant_range_wrt_TRP_range_interp,distances_over_arclengths_interp,20*log10(abs(wfm_AC_interp(1:length(distances_over_arclengths_interp),:)))-max_data);
    colormap('jet'); colorbar; caxis([min_image,max_image]);
    title('Focused Image');
    set(gca,'XLim',min_max_along,'FontSize',font_size);
    set(gca,'YLim',min_max_along,'FontSize',font_size);
    xlabel('Slant range w.r.t TRP position [m]'); ylabel('Along-track distance w.r.t TRP location [m]');
    %ylim(min_max_along); xlim(min_max_across);
    print('-dpng','-r150',strcat(plotsPath,'2D_PTR_interp_',filename_L1B_FFt,'.png'))
    
    tic
    % Along and across track cuts of the 2D PTR
    % along and across track cuts in the maximums
    [~,pos_max] = max(abs(wfm_AC_interp(:)));
    [pos_max_along_focused,pos_max_across_focused] = ind2sub(size(wfm_AC_interp),pos_max);
    cut_IRF_across = abs(wfm_AC_interp(pos_max_along_focused,:));
    cut_IRF_along = abs(wfm_AC_interp(1:length(distances_over_arclengths_interp),pos_max_across_focused));
    
    %Compute retracked range and height error
    R_retracked = (win_delay_surf(TRP_surf_idx)+(pos_max_across_focused-((chd.N_samples_sar/2)*cnf.FFt.zp_ac_TRP+1))*chd.T0_nom/cnf.FFt.zp_ac_TRP)*cst.c/2;
    H_measured = alt_sat_surf(TRP_surf_idx)-R_retracked;
    Height_error = abs(H_measured-chd.alt_trp);
    %Along-track location error (position of peak in along-track);
    Along_track_error = abs(distances_over_arclengths_interp(pos_max_along_focused));
    %Datation error: translate along-track error location to time with
    %the ground velocity of the beam
    %compute mean values of velocities and orbital height
    vs=mean(sqrt(sum(([x_vel_sat_sar_pulse.',y_vel_sat_sar_pulse.',z_vel_sat_sar_pulse.']).^2,2)));
    H_orb = mean(alt_sar_sat_pulse);
    vg = vs/(1+H_orb/cst.earth_radius);
    Datation_error = Along_track_error/vg;
    
    
    
    min_max_along = [-cnf.FFt.surf_dist_trp,cnf.FFt.surf_dist_trp];
    min_max_across = [-cnf.FFt.surf_dist_trp,+cnf.FFt.surf_dist_trp];
    toc
    figure;
    subplot(1,2,1);
    %resolution from 3dB
    [res_rg]=cc_ldb(20*log10(cut_IRF_across),slant_range_wrt_TRP_range_interp);
    % sidelobe levels
    [PSL_rg] = psl(cut_IRF_across.^2);
    plot(slant_range_wrt_TRP_range_interp,20*log10(cut_IRF_across)-max_data)
    set(gca,'XLim',min_max_across,'FontSize',font_size);
    set(gca,'YLim',[min_image,max_image],'FontSize',font_size);
    %figlabels('Shift [samples]','Beams','',['Range corrections' ' PT #' tag_pt(end:end) ],font_size);

    %xlim(min_max_across); ylim([min_image,max_image]);
    title(strcat('\delta_{ac}=',num2str(res_rg),{' [m] '}))%,...
    %                      'PSL_{l}=',num2str(PSL_rg(1)),{' [dB], '},...
    %                      'PSL_{r}=',num2str(PSL_rg(2)),{' [dB] '}));
    xlabel('Slant range w.r.t TRP position [m]'); ylabel('[dB]');
    subplot(1,2,2);
    plot(distances_over_arclengths_interp,20*log10(cut_IRF_along)-max_data)
    %resolution
    [res_az]=cc_ldb(20*log10(cut_IRF_along),distances_over_arclengths_interp);
    % sidelobe levels
    [PSL_az] = psl(cut_IRF_along.^2);
    set(gca,'XLim',min_max_across,'FontSize',font_size);
    set(gca,'YLim',[min_image,max_image],'FontSize',font_size);
    %xlim(min_max_along); ylim([min_image,max_image]);
    title(strcat('\delta_{al}=',num2str(res_az),{' [m] '}))%,...
    %                      'PSL_{l}=',num2str(PSL_az(1)),{' [dB], '},...
    %                      'PSL_{r}=',num2str(PSL_az(2)),{' [dB] '}));
    xlabel('Along-track distance w.r.t TRP location [m]'); ylabel('[dB]');
    print('-dpng','-r150',strcat(plotsPath,'CUTs_PTR_interp_',filename_L1B_FFt,'.png'));
    
    %Display numerical results of TRP analysis
    disp(strcat('Height error [m]: ',num2str(Height_error)));
    disp(strcat('Along-track error [m]: ',num2str(Along_track_error)));
    disp(strcat('Datation error [s]: ',num2str(Datation_error)));
    disp(strcat('Cross-track (Range) resolution [m]: ',num2str(res_rg)));
    disp(strcat('Along-track (azimuth) resolution [m]: ',num2str(res_az)));
    disp(strcat('PSL across-track [dB]: ',num2str(PSL_rg)));
    disp(strcat('PSL along-track [dB]: ',num2str(PSL_az)));
    
    %     save(strrep(filesBulk.filename_L1B_FFt,'.mat','_interp.mat'),...
    %     'wfm_AC_interp',...
    %     'lat_surf', 'lon_surf', 'alt_surf', 'win_delay_surf',...
    %     'time_surf',...
    %     'x_vel_sat_surf','y_vel_sat_surf','z_vel_sat_surf',...
    %     'alt_rate_sat_surf','alt_sat_surf',...
    %     'pitch_sat_surf','roll_sat_surf','yaw_sat_surf',...
    %     'distances_over_arclengths_interp','slant_range_wrt_TRP_range_interp',...
    %     'cut_IRF_across','cut_IRF_along',...
    %     '-v7.3')
    
    
end

end

%%TESTING, along track positioning without interpolating
%     % Along and across track cuts of the 2D PTR
%     % along and across track cuts in the maximums
%     [~,pos_max] = max(abs(wfm_AC(:)));
%     [pos_max_along_focused,pos_max_across_focused] = ind2sub(size(wfm_AC),pos_max);
%     cut_IRF_across = abs(wfm_AC(pos_max_along_focused,:));
%     cut_IRF_along = abs(wfm_AC(1:length(distances_over_arclengths),pos_max_across_focused));

%%PLOTS FOR ANTENNA DISTORTION
% Ku1
% figure;plot(slant_range_wrt_TRP_range_interp,20*log10(cut_IRF_across)-max_data)
% hold on;plot(slant_range_wrt_TRP_range_interp.*0.325, 20*log10(sinc(slant_range_wrt_TRP_range_interp)))
% Why doesn't align when multiplying x-axis by res_rg? ~0.03 bias

% Ka
% figure;plot(slant_range_wrt_TRP_range_interp,20*log10(cut_IRF_across)-max_data)
% hold on;plot(slant_range_wrt_TRP_range_interp.*res_rg, 20*log10(sinc(slant_range_wrt_TRP_range_interp)))
% Why doesn't align when multiplying x-axis by res_rg? ~0.03 bias

% Ka
% figure;plot(slant_range_wrt_TRP_range_interp,20*log10(cut_IRF_across)-max_data)
% hold on;plot(slant_range_wrt_TRP_range_interp.*0.33, 20*log10(sinc(slant_range_wrt_TRP_range_interp)))
% Why doesn't align when multiplying x-axis by res_rg? ~0.03 bias