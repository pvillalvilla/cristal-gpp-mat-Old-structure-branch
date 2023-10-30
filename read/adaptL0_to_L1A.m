%v1.0 First Version
%v1.1 Addn (Power, CoG to aned Calibratiotenna)
%v1.2 Added new cnf variables move_CoG2ant_alongtrack_cnf window_delay_source_cnf
function [L1A] = adaptL0_to_L1A(netCDF_L0,filename_L0, L1A, i_burst, cnf, chd, cst)
%v1.3 assigned ins_tracking_mode to netCDF_L0.data.instrument_mode_id


% global diff_wd_delay instr_delay_g
% L0         = [];
eccentricity = sqrt(cst.flat_coeff*(2-cst.flat_coeff));

L1A = setfield(L1A,'ins_id',netCDF_L0.data.instrument_mode_id);
L1A = setfield(L1A,'ins_tracking_mode',netCDF_L0.data.instrument_mode_id); %0=SAR_CB, 1=SAR_OB, 2=SARIN_CB 3=SARIN_OB
L1A = setfield(L1A,'inst_id_sar_isp',netCDF_L0.data.isp_type); % 0=LRM, 1=RAW, 2=RMC
% Check with ins_id when is defined if the processing mode is aligned with the data
L1A = setfield(L1A,'ProcessID',cnf.processing_mode); % read from Configuration File
L1A = setfield(L1A,'pri_sar',double(netCDF_L0.data.pulse_repetition_interval).*netCDF_L0.attributes.pulse_repetition_interval.scale_factor); %[seconds]
L1A = setfield(L1A,'bri_sar',double(netCDF_L0.data.burst_repetition_interval).*netCDF_L0.attributes.burst_repetition_interval.scale_factor); %[seconds]
L1A = setfield(L1A,'ambiguity_order_sar_isp',double(netCDF_L0.data.ambiguity_rank)); %[counts]
L1A = setfield(L1A,'burst',netCDF_L0.data.burst_counter);

pitch = double(netCDF_L0.data.off_nadir_pitch_angle).*netCDF_L0.attributes.off_nadir_pitch_angle.scale_factor ;
yaw = double(netCDF_L0.data.off_nadir_yaw_angle).*netCDF_L0.attributes.off_nadir_yaw_angle.scale_factor;
x_cog_ant= -1.632;  % [meters] TO BE Provided by Michele and read from CHD!!!!!!!!
x_cog_ant= -2.049-0.04997;  % [meters] TO BE Provided by Michele and read from CHD!!!!!!!!
x_cog_ant= 0.00;  % [meters] TO BE Provided by Michele and read from CHD!!!!!!!!
z_cog_ant= 0;  % [meters] TO BE Provided by Michele and read from CHD!!!!!!!!

if(cnf.move_CoG2ant_alongtrack)

    try
        com_pos = ncread(filename_L0,'com_position_vector',[1 i_burst],[3 2]); 
        x_pos_pair = com_pos(1,:);
        y_pos_pair = com_pos(2,:);
        z_pos_pair = com_pos(3,:);

        x_cog_ant_aux   = x_cog_ant * cos(pitch) + z_cog_ant * sin(pitch);
        x_cog_ant_corr  = x_cog_ant_aux * cos(yaw);
        
        x_vel   = double(netCDF_L0.data.com_velocity_vector(1)).*netCDF_L0.attributes.com_velocity_vector.scale_factor;
        y_vel   = double(netCDF_L0.data.com_velocity_vector(2)).*netCDF_L0.attributes.com_velocity_vector.scale_factor;
        z_vel   = double(netCDF_L0.data.com_velocity_vector(3)).*netCDF_L0.attributes.com_velocity_vector.scale_factor;
        vel     = sqrt((x_vel)^2+(y_vel)^2+(z_vel)^2); % m/s
        
        t_cog_ant_corr = x_cog_ant_corr / vel;
        corrected_time = double(netCDF_L0.data.data_record_time)+t_cog_ant_corr; % 'seconds since 2000-1-1 0:0:0'
    catch %last burst gives error, fill with the previous burst values
        com_pos = ncread(filename_L0,'com_position_vector',[1 i_burst-1],[3 2]); 
        x_pos_pair = com_pos(1,:);
        y_pos_pair = com_pos(2,:);
        z_pos_pair = com_pos(3,:);
        
        x_cog_ant_aux   = x_cog_ant * cos(pitch) + z_cog_ant * sin(pitch);
        x_cog_ant_corr  = x_cog_ant_aux * cos(yaw);
    end

    %lla = ecef2lla([x_pos_pair,y_pos_pair,z_pos_pair],cst.flat_coeff,cst.semi_major_axis); % error when computing the last one
    % test, take 1st value only (antena1)
    lla(1,:) = ecef2lla([x_pos_pair(1),y_pos_pair(1),z_pos_pair(1)],cst.flat_coeff,cst.semi_major_axis); % error when computing the last one
    lla(2,:) = ecef2lla([x_pos_pair(2),y_pos_pair(2),z_pos_pair(2)],cst.flat_coeff,cst.semi_major_axis); % error when computing the last one

    
    [arclen_v,az_v]     = distance('gc',lla(1,1),lla(1,2),lla(2,1),lla(2,2),[cst.semi_major_axis eccentricity]);
    [lats,lons]         = reckon(lla(1,1),lla(1,2),x_cog_ant_corr,az_v,[cst.semi_major_axis eccentricity]);
    lat_sat             = lats;
    lon_sat             = lons;
    alt_sat             = lla(1,3);
else
    
    lla = ecef2lla([double(netCDF_L0.data.com_position_vector(1)),double(netCDF_L0.data.com_position_vector(2)),double(netCDF_L0.data.com_position_vector(3))],cst.flat_coeff,cst.semi_major_axis);
    lat_sat = lla(1,1);
    lon_sat = lla(1,2);
    alt_sat = lla(1,3);
    
end


% lat_sat = lla(1,1);
% lon_sat = lla(1,2);

alt_sat2 = lla(1,3); 
L1A = setfield(L1A,'lat_sar_sat',lat_sat);
L1A = setfield(L1A,'lon_sar_sat',lon_sat);
% alt_sat = double(netCDF_L0.data.alt_l1a_echo)*double(netCDF_L0.attributes.alt_l1a_echo.scale_factor)+netCDF_L0.attributes.alt_l1a_echo.add_offset;
L1A = setfield(L1A,'alt_sar_sat',alt_sat);
%the difference using the alt_sat instead of alt_sat2 is -4.5009e-05 meters for the record cheched

alt_rate=double(netCDF_L0.data.com_altitude_rate)*double(netCDF_L0.attributes.com_altitude_rate.scale_factor);
L1A = setfield(L1A,'alt_rate_sar_sat',alt_rate);

L1A = setfield(L1A,'x_vel_sat_sar',double(netCDF_L0.data.com_velocity_vector(1)).*netCDF_L0.attributes.com_velocity_vector.scale_factor);
L1A = setfield(L1A,'y_vel_sat_sar',double(netCDF_L0.data.com_velocity_vector(2)).*netCDF_L0.attributes.com_velocity_vector.scale_factor);
L1A = setfield(L1A,'z_vel_sat_sar',double(netCDF_L0.data.com_velocity_vector(3)).*netCDF_L0.attributes.com_velocity_vector.scale_factor);
L1A = setfield(L1A,'roll_sar',    double(netCDF_L0.data.off_nadir_roll_angle).*netCDF_L0.attributes.off_nadir_roll_angle.scale_factor);
L1A = setfield(L1A,'pitch_sar',   double(netCDF_L0.data.off_nadir_pitch_angle).*netCDF_L0.attributes.off_nadir_pitch_angle.scale_factor);
L1A = setfield(L1A,'yaw_sar',     double(netCDF_L0.data.off_nadir_yaw_angle).*netCDF_L0.attributes.off_nadir_yaw_angle.scale_factor);
L1A = setfield(L1A,'confi_block_degraded',0); %double(netCDF_L1A.data.nav_bul_status_l1a_echo)
% L0 = setfield(L0,'mea_conf_data_fbr',0);


%'Distance between the altimeter reference point and the surface height 
% associated to a range gate used as reference inside the tracking window 
%(reference tracking point), corrected for USO frequency drift and internal path correction
if cnf.window_delay_source == 0 %L0
    HPR_b=floor(double(netCDF_L0.data.cor2_applied_l1a_echo)/chd.N_bursts_cycle_sar);
    h = zeros(1,chd.N_pulses_burst);
    for i_pulse=1:chd.N_pulses_burst
        h(i_pulse) = double(netCDF_L0.data.h0_applied_l1a_echo) + (double(netCDF_L0.data.burst_count_cycle_l1a_echo)-1) * floor(HPR_b/2^4);
    end

    CAI_tmp = floor(h(1) / 2^8); %[truncate result, algorithm b)]
    FAI_tmp = h(1) - CAI_tmp* 2^8;
    if (0 <= FAI_tmp)&& (FAI_tmp<= 127)
        FAI = FAI_tmp;
        CAI = CAI_tmp;
        rx_delay =  CAI *2^8*3.125/64*10^-9;

    elseif (128 <= FAI_tmp)&& (FAI_tmp<= 255)
        FAI = FAI_tmp - 256;
        CAI = CAI_tmp + 1;
        rx_delay =  CAI *2^8*3.125/64*10^-9;
    else
        CAI = 0;
        disp(['ERROR in CAI computation for record ' num2str(netCDF_L0.data.burst_count_cycle_l1a_echo) ' FAI= ' num2str(FAI_tmp)]);
        rx_delay = double(netCDF_L0.data.h0_applied_l1a_echo)* 3.125/64*10^-9 + double(netCDF_L0.data.cor2_applied_l1a_echo)*3.125/1024*10^-9 * (double(netCDF_L0.data.burst_count_cycle_l1a_echo)-1)/chd.N_bursts_cycle_sar;  %[s] scaling factor 3.125/64*10^-9 from attributes units
    end
    instr_delay = double(netCDF_L0.data.int_path_cor_ku_l1a_echo)   .* netCDF_L0.attributes.int_path_cor_ku_l1a_echo.scale_factor   + ...   %[m]
                  double(netCDF_L0.data.uso_cor_l1a_echo)           .* netCDF_L0.attributes.uso_cor_l1a_echo.scale_factor +...    %[m]
                  double(netCDF_L0.data.cog_cor_l1a_echo)           .* netCDF_L0.attributes.cog_cor_l1a_echo.scale_factor .* (cos(pitch));                  %[m]

    %           z_cog_ant_corr(i_burst,i_pulse) = x_cog_ant * sin(pitch_pre_dat(i_burst)) + z_cog_ant * cos(pitch_pre_dat(i_burst));      
    %           delta_pitch(i_burst,i_pulse) = - 2 * z_cog_ant_corr(i_burst,i_pulse) / cst.c;

    range_delay = rx_delay * cst.c/2 + instr_delay; %[m]
    % diff_wd_delay = [diff_wd_delay range_delay_2-range_delay];
    % instr_delay_g = [instr_delay_g ; double(netCDF_L1A.data.int_path_cor_ku_l1a_echo)   .* netCDF_L1A.attributes.int_path_cor_ku_l1a_echo.scale_factor   ...   %[m]
    %               double(netCDF_L1A.data.uso_cor_l1a_echo)           .* netCDF_L1A.attributes.uso_cor_l1a_echo.scale_factor           ...    %[m]
    %               double(netCDF_L1A.data.cog_cor_l1a_echo)           .* netCDF_L1A.attributes.cog_cor_l1a_echo.scale_factor];

else %L1A
    range_delay = double(netCDF_L0.data.tracker_range_calibrated)*double(netCDF_L0.attributes.tracker_range_calibrated.scale_factor)+netCDF_L0.attributes.tracker_range_calibrated.add_offset;
%     instr_delay = double(netCDF_L0.data.cog_cor_l1a_echo)           .* netCDF_L0.attributes.cog_cor_l1a_echo.scale_factor .* (cos(pitch));
%     range_delay = range_delay + instr_delay;
end
    L1A = setfield(L1A,'win_delay',range_delay/cst.c*2);

L1A = setfield(L1A,'h0_comp_sar_isp',double(netCDF_L0.data.h0_initial_altitude_instruction));%[seconds]
L1A = setfield(L1A,'cor2_comp_sar_isp',double(netCDF_L0.data.cor2_altitude_rate_estimation));%[seconds]
% L1A = setfield(L1A,'cai_sar_isp',double(netCDF_L0.data.coarse_altitude_instruction));%[seconds]
% L1A = setfield(L1A,'fai_sar_isp',double(netCDF_L0.data.fine_altitude_instruction));
% L0 = setfield(L0,'ATT1_science',-1.0*double(netCDF_L0.data.agc_ku_l1a_echo)*double(netCDF_L0.attributes.agc_ku_l1a_echo.scale_factor)+netCDF_L0.attributes.agc_ku_l1a_echo.add_offset);
% L0 = setfield(L0,'ATT2_science',0);
% L0 = setfield(L0,'att_isp',62-L0.ATT1_science); %as done in CR2
% L0 = setfield(L0,'tot_fixed_gain_1',0);
% L0 = setfield(L0,'tot_fixed_gain_2',0);
% L0 = setfield(L0,'transmit_power',0);
L1A = setfield(L1A,'doppler_range_correction',0);
% L0 = setfield(L0,'instrument_range_correction_tx_rx',double(netCDF_L0.data.int_path_cor_ku_l1a_echo)*double(netCDF_L0.attributes.int_path_cor_ku_l1a_echo.scale_factor)+netCDF_L0.attributes.int_path_cor_ku_l1a_echo.add_offset);
% L0 = setfield(L0,'instrument_range_correction_rx',0);
% L1A = setfield(L1A,'instrument_sigma0_correction_tx_rx',double(netCDF_L0.data.power_scaling_to_antenna)*double(netCDF_L0.attributes.power_scaling_to_antenna.scale_factor)); %SIGMA0 CAL
% L0 = setfield(L0,'instrument_sigma0_correction_rx',double(netCDF_L0.data.scale_factor_ku_l1a_echo)*double(netCDF_L0.attributes.scale_factor_ku_l1a_echo.scale_factor)+netCDF_L0.attributes.scale_factor_ku_l1a_echo.add_offset); %[dB] %SCALE FACTOR SIGMA0 , GOES TO SCALING FACTOR ALGORITHM
L1A = setfield(L1A,'internal_phase_correction',0);
L1A = setfield(L1A,'external_phase_correction',0);
L1A = setfield(L1A,'noise_power',0);
L1A = setfield(L1A,'phase_slope_correction',0);
L1A = setfield(L1A,'gain_corr_instr_sar_1',0);


L1A = setfield(L1A,'dry_tropo_correction_bursts',0);
L1A = setfield(L1A,'wet_tropo_correction_bursts',0);
L1A = setfield(L1A,'inverse_baro_correction_bursts',0);
L1A = setfield(L1A,'Dynamic_atmospheric_correction_bursts',0);
L1A = setfield(L1A,'GIM_iono_correction_bursts',0);
L1A = setfield(L1A,'model_iono_correction_bursts',0);
L1A = setfield(L1A,'ocean_equilibrium_tide_bursts',0);
L1A = setfield(L1A,'long_period_tide_height_bursts',0);
L1A = setfield(L1A,'ocean_loading_tide_bursts',0);
L1A = setfield(L1A,'solid_earth_tide_bursts',0);
L1A = setfield(L1A,'geocentric_polar_tide_bursts',0);
L1A = setfield(L1A,'surface_type_flag_bursts',0);

L1A = setfield(L1A,'nimp_sar_isp',netCDF_L0.data.ambiguity_rank); %[counts]);



i_pulse = double(netCDF_L0.data.rx1_complex_waveforms_i_samples);
if(sum(isnan(i_pulse)))
  disp('NaNs found')
end
q_pulse = double(netCDF_L0.data.rx1_complex_waveforms_q_samples);
if(sum(isnan(q_pulse)))
  disp('NaNs found')
end
wfm_cal_gain_CORRECTED = (i_pulse+1i.*q_pulse).';
L1A = setfield(L1A,'wfm_cal_gain_corrected',wfm_cal_gain_CORRECTED);

if(strcmp(cnf.processing_mode,'SIN'))
    i_pulse_2 = double(netCDF_L0.data.rx2_complex_waveforms_i_samples);
    q_pulse_2 = double(netCDF_L0.data.rx2_complex_waveforms_q_samples);
    wfm_cal_gain_CORRECTED_2 = (i_pulse_2+1i.*q_pulse_2).';
    L1A = setfield(L1A,'wfm_cal_gain_corrected_2',wfm_cal_gain_CORRECTED_2);
end


if cnf.height_rate_application
    %Assuming FAI is not applied on-board
    [L1A] = height_rate_alignment(L1A, cnf, chd, cst);  
end


lat_sar_surf = lat_sat;
lon_sar_surf = lon_sat;
alt_sar_surf = alt_sat - range_delay;


% geod2cart(SURF)
p = lla2ecef([lat_sar_surf,lon_sar_surf,alt_sar_surf],cst.flat_coeff,cst.semi_major_axis);
x_sar_surf = p(1).';
y_sar_surf = p(2).';
z_sar_surf = p(3).';


L1A = setfield(L1A,'time_rx_1st',double(netCDF_L0.data.data_record_time)); % time of the first RX pulse of the burst
L1A = setfield(L1A,'x_sar_sat',double(netCDF_L0.data.com_position_vector(1)));
L1A = setfield(L1A,'y_sar_sat',double(netCDF_L0.data.com_position_vector(2)));
L1A = setfield(L1A,'z_sar_sat',double(netCDF_L0.data.com_position_vector(3)));
L1A = setfield(L1A,'lat_sar_surf',lat_sar_surf);
L1A = setfield(L1A,'lon_sar_surf',lon_sar_surf);
L1A = setfield(L1A,'alt_sar_surf',alt_sar_surf);
L1A = setfield(L1A,'x_sar_surf',x_sar_surf);
L1A = setfield(L1A,'y_sar_surf',y_sar_surf);
L1A = setfield(L1A,'z_sar_surf',z_sar_surf);
[~,doppler_ang_sar_sat] = compute_height_rate(1, L1A.x_vel_sat_sar, L1A.y_vel_sat_sar, L1A.z_vel_sat_sar, L1A.x_sar_sat, L1A.y_sar_sat, L1A.z_sar_sat,L1A.x_sar_surf,L1A.y_sar_surf,L1A.z_sar_surf, cst);


L1A = setfield(L1A,'doppler_ang_sar_sat',doppler_ang_sar_sat);

L1A = setfield(L1A,'T0_sar',chd.T0_nom);


%beam angles
L1A = setfield(L1A,'N_beams',0);
L1A = setfield(L1A,'surf_loc_index',0);
L1A = setfield(L1A,'beam_ang',0);
L1A = setfield(L1A,'beam_ang_nadir_index',0);
L1A = setfield(L1A,'beam_ang_index',0);
L1A = setfield(L1A,'start_beam',0);
L1A = setfield(L1A,'end_beam',0);
L1A = setfield(L1A,'beam_index',0);

%azimuth processing
L1A = setfield(L1A,'beams_focused_shifted',zeros(chd.N_pulses_burst,chd.N_samples_sar));
if(strcmp(cnf.processing_mode,'SIN'))
    L1A = setfield(L1A,'beams_focused_shifted_2',zeros(chd.N_pulses_burst,chd.N_samples_sar));
end
p = lla2ecef([L1A.lat_sar_sat.',L1A.lon_sar_sat.',L1A.alt_sar_sat.'],cst.flat_coeff,cst.semi_major_axis);
L1A.x_sar_sat = p(:,1).';
L1A.y_sar_sat = p(:,2).';
L1A.z_sar_sat = p(:,3).';

L1A.lat_sar_surf = L1A.lat_sar_sat;
L1A.lon_sar_surf = L1A.lon_sar_sat;
L1A.alt_sar_surf = L1A.alt_sar_sat - L1A.win_delay * cst.c/2;

% geod2cart(SURF)
p = lla2ecef([L1A.lat_sar_surf,L1A.lon_sar_surf,L1A.alt_sar_surf],cst.flat_coeff,cst.semi_major_axis);
L1A.x_sar_surf = p(1).';
L1A.y_sar_surf = p(2).';
L1A.z_sar_surf = p(3).';



end

