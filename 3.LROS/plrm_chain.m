% plrm_chain_clean
function L1B = plrm_chain(L1A)
    global N_samples N_ku_pulses_burst_chd zp_fact_range_cnf
    global N_bursts_cycle_chd pi_cst
    global chirp_slope_ku_chd wv_length_ku c_cst
    global power_tx_ant_ku_chd antenna_gain_ku_chd pulse_length_chd
    global earth_radius_cst
    global ref_burst fromSARin
    
    if size(L1A,2) == N_bursts_cycle_chd
        wfm_L1A_fft_sq_av1                  = zeros(N_bursts_cycle_chd,N_samples*zp_fact_range_cnf);
        if fromSARin
            wfm_L1A_fft_sq_av2              = zeros(N_bursts_cycle_chd,N_samples*zp_fact_range_cnf);
        end
        sigma0_scaling_factor_lrm_ku_burst  = zeros(1,N_bursts_cycle_chd);
        altitude_rate                       = zeros(1,N_bursts_cycle_chd);
        roll                                = zeros(1,N_bursts_cycle_chd);
        pitch                               = zeros(1,N_bursts_cycle_chd);
        yaw                                 = zeros(1,N_bursts_cycle_chd);
        samples                             = 0:(N_samples-1);
        pulses_and_samples                  = repmat(samples,N_ku_pulses_burst_chd,1);
        doppler_seconds                     = zeros(1,N_bursts_cycle_chd);
        wd_shift1                           = zeros(1,N_bursts_cycle_chd);
        max_value_burst                     = zeros(1,N_bursts_cycle_chd);
        min_value_burst                     = zeros(1,N_bursts_cycle_chd);
        
        for i_burst = 1:N_bursts_cycle_chd
            %% Compute Doppler correction
            doppler_seconds(i_burst)        = - 2 * L1A(i_burst).alt_rate_sar_sat / (wv_length_ku * chirp_slope_ku_chd); %in seconds
            
            L1A(i_burst).win_delay_sar_ku   = L1A(i_burst).win_delay_sar_ku + doppler_seconds(i_burst);
        end %Two different loops so the ref_burst window delay always has the doppler correction applied to it, even if we change the 'ref_burst'
        
        for i_burst = 1:N_bursts_cycle_chd
            if fromSARin
                wfm_L1A1 = squeeze(L1A(i_burst).wfm_cal_gain_corrected);
                wfm_L1A2 = squeeze(L1A(i_burst).wfm_cal_gain_corrected_2);
            else
                %% Compute window delay shift
                wd_shift1(i_burst)      = ((L1A(i_burst).win_delay_sar_ku - L1A(ref_burst).win_delay_sar_ku) - (L1A(i_burst).alt_sar_sat - L1A(ref_burst).alt_sar_sat)*2/c_cst) ./ L1A(i_burst).T0_sar;
                
                wfm_L1A_aux = squeeze(L1A(i_burst).wfm_cal_gain_corrected);
                if wfm_L1A_aux == 0
                    i_burst = i_burst;
                end
                %% Apply window delay shift
                wfm_L1A1 = wfm_L1A_aux .* exp(2i*pi_cst/N_samples*wd_shift1(i_burst).*pulses_and_samples);
            end
            
            %% FFT range and power waveforms
            wfm_L1A_fft1 = abs(fftshift(fft(wfm_L1A1.',N_samples*zp_fact_range_cnf),1).').^2/(N_samples*zp_fact_range_cnf);
            if fromSARin
                wfm_L1A_fft2 = abs(fftshift(fft(wfm_L1A2.',N_samples*zp_fact_range_cnf),1).').^2/(N_samples*zp_fact_range_cnf);
            end
            
            
            %% Intraburst AVERAGE
            wfm_L1A_fft_sq_av1(i_burst,:)       = mean(wfm_L1A_fft1);
            if fromSARin
                wfm_L1A_fft_sq_av2(i_burst,:)       = mean(wfm_L1A_fft2);
            end
            
            %% max values
            max_value_burst(i_burst)    = max(max(wfm_L1A_fft1.'));
            min_value_burst(i_burst)    = min(max(wfm_L1A_fft1.'));
            max_value(i_burst)          = max(wfm_L1A_fft_sq_av1(i_burst,:));
            
            %% Waveform scaling factor
            sigma0_scaling_factor_lrm_ku_burst(i_burst) = ...
                + 30*log10(4)...
                + 30*log10(L1A(i_burst).win_delay_sar_ku*0.5*c_cst)...
                + 20*log10(pi_cst)...
                + 10*log10(pulse_length_chd * chirp_slope_ku_chd)...
                + 10*log10(1 + (L1A(i_burst).win_delay_sar_ku*0.5*c_cst)/earth_radius_cst)...
                - 10*log10(power_tx_ant_ku_chd)...
                - 2*antenna_gain_ku_chd...
                - 20*log10(wv_length_ku) ...
                - 10*log10(c_cst)...
                + 7.25; %bias between pLRM and SAR
            
            
            %% Parameters to average
            altitude_rate(i_burst)   = L1A(i_burst).alt_rate_sar_sat;
            roll(i_burst)            = L1A(i_burst).roll_sar;
            pitch(i_burst)           = L1A(i_burst).pitch_sar;
            yaw(i_burst)             = L1A(i_burst).yaw_sar;
            
            dry_tropo_correction(i_burst)               = L1A(i_burst).dry_tropo_correction_bursts;
            wet_tropo_correction(i_burst)               = L1A(i_burst).wet_tropo_correction_bursts;
            inverse_baro_correction(i_burst)            = L1A(i_burst).inverse_baro_correction_bursts;
            Dynamic_atmospheric_correction(i_burst)     = L1A(i_burst).Dynamic_atmospheric_correction_bursts;
            GIM_iono_correction(i_burst)                = L1A(i_burst).GIM_iono_correction_bursts;
            model_iono_correction(i_burst)              = L1A(i_burst).model_iono_correction_bursts;
            ocean_equilibrium_tide(i_burst)             = L1A(i_burst).ocean_equilibrium_tide_bursts;
            long_period_tide_height(i_burst)            = L1A(i_burst).long_period_tide_height_bursts;
            ocean_loading_tide(i_burst)                 = L1A(i_burst).ocean_loading_tide_bursts;
            solid_earth_tide(i_burst)                   = L1A(i_burst).solid_earth_tide_bursts;
            geocentric_polar_tide(i_burst)              = L1A(i_burst).geocentric_polar_tide_bursts;
            
        end
        
        %% Prepare Writing
        if fromSARin
            ref_burst = 1;
        end
        %---------- A. Time variables -----------------------------------------
        L1B.time_l1b_plrm                   = L1A(ref_burst).days + L1A(ref_burst).seconds + L1A(ref_burst).microseconds; %to check units of microseconds
        
        %---------- B. Orbit and attitude variables ---------------------------
        L1B.x                               = L1A(ref_burst).x_sar_sat;
        L1B.y                               = L1A(ref_burst).y_sar_sat;
        L1B.z                               = L1A(ref_burst).z_sar_sat;
        
        L1B.x_vel                           = L1A(ref_burst).x_vel_sat_sar;
        L1B.y_vel                           = L1A(ref_burst).y_vel_sat_sar;
        L1B.z_vel                           = L1A(ref_burst).z_vel_sat_sar;
        
        L1B.latitude                        = L1A(ref_burst).lat_sar_sat;
        L1B.longitude                       = L1A(ref_burst).lon_sar_sat;
        if (L1B.longitude > 180)
            L1B.longitude                   = L1B.longitude - 360;
        end
        L1B.altitude                        = L1A(ref_burst).alt_sar_sat;
        L1B.altitude_rate                   = mean(altitude_rate);
        
        L1B.roll                            = mean(roll);
        L1B.pitch                           = mean(pitch);
        L1B.yaw                             = mean(yaw);
        
        L1B.doppler                         = mean(doppler_seconds);
        L1B.win_delay_sar                   = L1A(ref_burst).win_delay_sar_ku; %Doppler has already been compensated (subracted from the WD) before
        L1B.T0_sar                          = L1A(ref_burst).T0_sar;
        
        
        %% Interburst AVERAGE
        if fromSARin == 0
            L1B.wfm_L1B_RC                  = mean(wfm_L1A_fft_sq_av1);
        else
            L1B.wfm_L1B_RC                  = mean([wfm_L1A_fft_sq_av1;wfm_L1A_fft_sq_av2]);
        end
        
        L1B.sigma0_scaling_factor_RC        = mean(sigma0_scaling_factor_lrm_ku_burst);
        
        
        %----------N. Geophysical Corrections variables -------------------------
        L1B.dry_tropo_correction            = mean(dry_tropo_correction(i_burst));
        L1B.wet_tropo_correction            = mean(wet_tropo_correction(i_burst));
        L1B.inverse_baro_correction         = mean(inverse_baro_correction(i_burst));
        L1B.Dynamic_atmospheric_correction  = mean(Dynamic_atmospheric_correction(i_burst));
        L1B.GIM_iono_correction             = mean(GIM_iono_correction(i_burst));
        L1B.model_iono_correction           = mean(model_iono_correction(i_burst));
        L1B.ocean_equilibrium_tide          = mean(ocean_equilibrium_tide(i_burst));
        L1B.long_period_tide_height         = mean(long_period_tide_height(i_burst));
        L1B.ocean_loading_tide              = mean(ocean_loading_tide(i_burst));
        L1B.solid_earth_tide                = mean(solid_earth_tide(i_burst));
        L1B.geocentric_polar_tide           = mean(geocentric_polar_tide(i_burst));
    
    
    end
