%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT
% --------------------------------------------------------
% Sentinel 3 Next Generation Topography
% ----------------------------------------------------------
% Author:    Pol Villalvilla-Ornat  / isardSAT
%
% v1.0 21/08/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READING
wfm_AC = zeros(N_bursts,256);
beam_ang = linspace(89.9*pi/180,90.5*pi/180,N_bursts);

for i_burst=1:N_bursts
    if i_burst==1
        fidL1A = netcdf.open(filesBulk.filename_L0,'NC_NOWRITE'); %open file
    end
    netCDF_L1A = readanyNETCDF_record(fidL1A,i_burst,'nb');
    L1A = netCDF_L1A.data;
    if i_burst==N_bursts
        netcdf.close(fidL1A); %close the netcdf file
    end
    %% adaptL1A data
    altimeter_clock = 1/chd.T0_nom;
    position_vector = L1A.com_position_vector.';
    altitude_rate = double(L1A.com_altitude_rate) * netCDF_L1A.attributes.com_altitude_rate.scale_factor;
    velocity_vector = double(L1A.com_velocity_vector) * netCDF_L1A.attributes.com_velocity_vector.scale_factor;
    off_nadir_pitch_angle = double(L1A.off_nadir_pitch_angle) * netCDF_L1A.attributes.off_nadir_pitch_angle.scale_factor;
    off_nadir_roll_angle = double(L1A.off_nadir_roll_angle) * netCDF_L1A.attributes.off_nadir_roll_angle.scale_factor;
    off_nadir_yaw_angle = double(L1A.off_nadir_yaw_angle) * netCDF_L1A.attributes.off_nadir_yaw_angle.scale_factor;
    power_scaling_to_antenna = double(L1A.power_scaling_to_antenna) * netCDF_L1A.attributes.power_scaling_to_antenna.scale_factor;
    pri = double(L1A.pulse_repetition_interval) * netCDF_L1A.attributes.pulse_repetition_interval.scale_factor;
    tracker_range_calibrated = double(L1A.tracker_range_calibrated) * netCDF_L1A.attributes.tracker_range_calibrated.scale_factor + netCDF_L1A.attributes.tracker_range_calibrated.add_offset;
    
    
    %%
    wfm_iq = L1A.rx1_complex_waveforms_i_samples.' + 1i*L1A.rx1_complex_waveforms_q_samples.';
    N_samples_az = size(wfm_iq,1);
    cnf.zp_fact_azimuth = 1;
    cnf.zp_fact_range_uf_cnf = 1;
    
    
    %% Preliminary Datation
    
    
    %% Window Delay
    
    
    %% Final Datation
    
    
    %% Instrument Gain
    
    
    %% Waveforms Correction
    
    
    %% Azimuth procesing
    wfm_beams_focused = fft(fftshift(wfm_iq,1),N_samples_az*cnf.zp_fact_azimuth,1)/sqrt(N_samples_az*cnf.zp_fact_azimuth);
    
    
    
    %% Geo Corrections
    apply_doppler = 0;
    apply_slant_range = 0;
    
    
    %%% Doppler Correction
    
    norm_vel_sat = norm(velocity_vector);
    % Range correction computation
    doppler_range = (- cst.c / chd.wv_length .* norm_vel_sat .* cos(beam_ang(i_burst))) .* chd.pulse_length / chd.bw;
    % Doppler correction computation
    T0_sar_surf = chd.T0_nom;
    doppler_corr = 2 / cst.c .* doppler_range ./ T0_sar_surf;
    %         L1BS.doppler_corr(i_beam) = 0; %testing
    doppler_corr(isnan(doppler_corr))=0;
    
    %%% Slant range correction
    lla = ecef2lla(position_vector,cst.flat_coeff,cst.semi_major_axis);
    lat_sat = lla(1,1);
    lon_sat = lla(1,2);
    altitude = lla(1,3);
    
    lat_surf = lat_sat;
    lon_surf = lon_sat;
    alt_surf = altitude - tracker_range_calibrated;
    
    position_trp = [chd.x_trp, chd.y_trp, chd.z_trp];
    
    range_sat_surf = norm(position_vector - position_trp);
    % Range delay computation
    slant_range_corr_time = (tracker_range_calibrated * 2 / cst.c) - (range_sat_surf * 2 / cst.c);
    slant_range_corr = slant_range_corr_time / T0_sar_surf;
    
    shift = -doppler_corr*apply_doppler + slant_range_corr*apply_slant_range;
    % the Doppler correction Sign is changed after teleconf with Michele as they are generating the sinusoidal with a negative carrier frequency
    
    
    
    shift_coarse = round(shift);
    shift_fine = shift - shift_coarse;
    
    % Apply corrections
    i_samples               = 0:(chd.N_samples_sar-1);
    beam_geo_corr = wfm_beams_focused .* exp(2i*cst.pi./chd.N_samples_sar.*shift*i_samples);
    
    
    
    
    %% Range compression
    
    non_centered_spectra=fftshift(beam_geo_corr(:,:),2);
    beams_zp_fft = fft([non_centered_spectra(:,1:chd.N_samples_sar/2),...
        zeros(N_samples_az,(cnf.zp_fact_range_uf_cnf-1)*chd.N_samples_sar),...
        non_centered_spectra(:,chd.N_samples_sar/2+1:chd.N_samples_sar)],...
        chd.N_samples_sar * cnf.zp_fact_range_uf_cnf,2);
    
    % figure; imagesc(20*log10(abs(beams_zp_fft)))
    % colormap('jet'); colorbar;
    
    %% Azimuth procesing
    % wfm_beams_focused = fft(fftshift(beams_zp_fft,1),N_samples_az*cnf.zp_fact_azimuth,1)/sqrt(N_samples_az*cnf.zp_fact_azimuth);
    
    
    %% Multilooking
    N_av_SL = numel(beams_zp_fft(:,1)) - sum(isnan(beams_zp_fft(:,1)));
    wfm_AC(i_burst,:) = nansum(beams_zp_fft)/N_av_SL';
    
    
    
    %%
    % Power waveform
    % beams_rng_cmpr = abs(wfm_beams_focused).^2;
    %
    % figure; imagesc(20*log10(abs(beams_rng_cmpr)))
    % colormap('jet'); colorbar;
    
    % %%
    %
    % % zp_al_TRP = 16;
    % % zp_ac_TRP = 16;
    % % wfm_AC_interp=interpolsinc_2D(wfm_beams_focused,zp_ac_TRP,zp_al_TRP);
    % max_data  = max(20*log10(abs(wfm_beams_focused(:))));
    % max_image = 0.0;
    % min_image = max_image-60.0;
    %
    % figure; imagesc(20*log10(abs(wfm_beams_focused))-max_data)
    % colormap('jet'); colorbar; caxis([min_image,max_image]);
    %
    %
    % [~,pos_max] = max(abs(wfm_beams_focused(:)));
    %         [pos_max_along_focused,pos_max_across_focused] = ind2sub(size(wfm_beams_focused),pos_max);
    %         cut_IRF_across = abs(wfm_beams_focused(pos_max_along_focused,:));
    %         cut_IRF_along = abs(wfm_beams_focused(:,pos_max_across_focused));
    %
    % figure; plot(cut_IRF_across)
    % figure; plot(cut_IRF_along)
    %
    % a
    
    %% Sigma-0 Scaling Factor
    
end
[~,pos_max] = max(max(wfm_AC));

step_size = tracker_range_calibrated / pos_max;

% Create an array of range values
range_array = (0:step_size:(tracker_range_calibrated-step_size));

% If needed, fill the rest of the array with values beyond tracker_range
range_array = [range_array, (tracker_range_calibrated+step_size):step_size:(step_size*(chd.N_samples_sar*cnf.zp_fact_range_uf_cnf))];

figure; imagesc(range_array, 1:N_bursts,20*log10(abs(wfm_AC)))
figure; plot(range_array,(sum(abs(wfm_AC),1)).^2)

%% Estimate coordinates to plot
L1Aall = readanyNETCDF_V3(filesBulk.filename_L0);
position_vector = L1Aall.data.com_position_vector;
for kk=1:N_bursts
    lla = ecef2lla(position_vector(:,kk).',cst.flat_coeff,cst.semi_major_axis);
    lat_sat(kk) = lla(1,1);
    lon_sat(kk) = lla(1,2);
    altitude(kk) = lla(1,3);
end

position_TRP = [chd.x_trp, chd.y_trp, chd.z_trp];
lla = ecef2lla(position_TRP,cst.flat_coeff,cst.semi_major_axis);
lat_sat_trp = lla(1,1);
lon_sat_trp = lla(1,2);
altitude_trp = lla(1,3);

figure; plot(lat_sat, lon_sat)
hold on
plot(lat_sat_trp, lon_sat_trp,'o')
xlabel('latitude')
ylabel('longitude')
grid
legend('Satellite track','PT')


figure; scatter3(lat_sat, lon_sat,altitude_trp*ones(1,1969))
hold on; scatter3(lat_sat_trp, lon_sat_trp,altitude_trp)
xlabel('latitude')
ylabel('longitude')
zlabel('altitude')

legend('Satellite track','PT')

%%
for kk=1:N_bursts
    ecef = lla2ecef([lat_sat(kk),lon_sat(kk),altitude_trp],cst.flat_coeff,cst.semi_major_axis);
    x_g(kk) = ecef(1,1);
    y_g(kk) = ecef(1,2);
    z_g(kk) = ecef(1,3);
end

across_track_d = min(sqrt( (x_g-position_TRP(1)).^2 + (y_g-position_TRP(2)).^2 + (z_g-position_TRP(3)).^2) )




