%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% 
% This code implements the CODING & PACKING 
% algorithm for Level-1BS HR CRISTAL product 
%
% ---------------------------------------------------------
% Objective: Pack variables and write the NETCDF
% 
% INPUTs : Workspace
% OUTPUTs: NETCDF file with structure as defined in CRIS-DS-ISR-GS-0007-Annex-A_v5a_PFS_L1
%
% ----------------------------------------------------------
% Author:    Juan Pedro López-Zaragoza / isardSAT
% Reviewer:  
% Last rev.: 
% 
% Versions
% 1.0 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils

function write_NetCDF_L1Bs_HR(files,L1A_buffer,L1BS, L1B,i_surf_stacked,cnf,chd,cst)
i_surf_stacked= i_surf_stacked-1; % Because netcdf file variables start at 0, while MATLAB at 1

%{
i_surf = i_surf_stacked-1;
global cst.c chd.bw_ku_chd  cnf.zp_fact_range cst.sec_in_day cst.pi
global mission cnf.mode chd.N_samples_sar N_max_beams_stack_chd cnf.processing_mode
global cnf.compute_L1_stadistics
global cnf.include_wfms_aligned
%added EM: 04.10.2016
global ACDC_application_cnf

%}
%global netcdf_type
% ncid = netcdf.open(files.filename_netCDF,'WRITE');

% ku_rec_dimension = netcdf.defDim(ncid,'time_l1b_echo',N_surfs_loc_estimated);
% nl_dimension = netcdf.defDim(ncid,'Nl',N_max_beams_stack_chd);
% ns_dimension = netcdf.defDim(ncid,'Ns',chd.N_samples_sar*cnf.zp_fact_range);
% space_3D_dimension = netcdf.defDim(ncid,'space_3D',3);

% t6 = tic;

% date_creation = datestr(now, '_yyyymmdd_');
% switch mission
%     case 'CR2'        
%         switch cnf.mode
%             case 'SAR'                
%                 files.filename_netCDF = strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA____',...
%                                 files.sph.product_info.product_id(20:20+30),...
%                                 date_creation,...
%                                 'isd','.nc');
%         end
%     case {'S3_','S3A','S3B'} 
%     case {'S6_'} 
% end


% dimensions_key = 'Dimensions';
% format_key = 'Format';
% data_type_key = 'DataType';
% fill_value_key='FillValue';
% 
% netcdf_v4_format = 'netcdf4';
% 
% ku_rec_dimension = 'time_l1b_echo');
% nl_dimension = 'Nl';
 nl_dimension_size = chd.N_max_beams_stack;
% ns_dimension = 'Ns';
ns_dimension_size=chd.N_samples_sar*cnf.zp_fact_range;
% space_3D_dimension = 'space_3D';
 space_3D_dimension_size = 3;
np_dimension_size= chd.N_pulses_burst;


% int8_type = 'int8';
% uint8_type = 'uint8';
% int16_type = 'int16';
% uint16_type = 'uint16';
% int32_type = 'int32';
% uint32_type = 'uint32';
% uint64_type = 'uint64';
% float_type = 'single';
% double_type='double';


%% CODING L1B

%----------A. Time & counters variables ----------------------------------------------

l1b_record_counter=int32(L1BS.surf_counter);
time = double(L1BS.time_surf);
% (TO DO!) Put time in TAI, we are using UTC for the time being
time_tai=double(L1BS.time_surf);
% (TO DO!) Copy value from L1A when it is written
%tm_source_sequence_counter=int16(L1A.tm_source_sequence_counter) ;

%----------B. Orbit and attitude variables ------------------------------

altitude=int32((L1BS.alt_sat-700000) * 1e4);
altitude_rate=int32(L1BS.alt_rate_sat * 1e4);
latitude=int32(L1BS.lat_surf * 1e6); % we save the surface lat/lon. The difference wrt to the sat ones is minimal
longitude=int32(L1BS.lon_surf * 1e6);
% (TO DO!) Misspointing angles still to implement when we get the attitude file
% off_nadir_roll_angle_ant1
% off_nadir_pitch_angle_ant1
% off_nadir_yaw_angle_ant1
% off_nadir_roll_angle_ant2
% off_nadir_pitch_angle_ant2
% off_nadir_yaw_angle_ant2
% off_nadir_roll_angle
% off_nadir_pitch_angle
% off_nadir_yaw_angle
%(TO DO!) Put orbit type once we get the orbit file
%orbit_type_flag
position_vector=double([L1BS.x_sat; L1BS.y_sat; L1BS.z_sat]);
velocity_vector=int32(1e4.*[L1BS.x_vel_sat' L1BS.y_vel_sat' L1BS.z_vel_sat']);

%----------C. Configuration and quality variables --------------------------------------

% (TO DO!) Implement these flags
% mcd_flags
% processing_configuration_flags
% hr_mcd_flags
% hr_processing_configuration_flags
% iris_instrument_configuration_flags
% iris_mode_flag
% range_oversampling_factor

%----------------D. Altimeter range variables ---------------------------

% (TO DO!) Add all the corrections range related variables
% range_cor_com_ant1
% range_cor_com_ant2
% range_cor_com
% range_cor_doppler
% range_cor_internal_delay_cal_rx1
% range_cor_internal_delay_cal_rx2
% range_cor_internal_delay_att_rx1
% range_cor_internal_delay_att_rx2
% range_cor_uso
% range_cor_reference_sample
%tracker_range_calibrated=int32(L1BS.win_delay_surf*cst.c/2-700000)*1e4;
tracker_range_calibrated=(L1BS.win_delay_surf*cst.c/2-700000)*1e4;

%tracker_range_calibrated_gnss

%----------E. Altimeter power variables ------------------------------

% (TO DO!) Implement all the altimeter power related variables
% altimeter_power_drift_rx1
% altimeter_power_drift_rx2
% attenuation_calibrated_rx1
% attenuation_calibrated_rx2
% power_scaling_to_antenna_rx1
% power_scaling_to_antenna_rx2
% variable_digital_gain
% cal1_power_rx1
% cal1_power_rx2

%----------F. Altimeter engineering variables --------------------------
[~,beam_index_nadir]=min(abs(-pi/2+L1BS.beam_ang_surf'));
burst_index_nadir= L1BS.burst_index(beam_index_nadir);

T0_sar_surf_nadir = L1BS.T0_sar_surf(beam_index_nadir);
altimeter_clock = int32(1./T0_sar_surf_nadir - 6e8)* 1e6;
% (TO DO!) Implement mean altitude instruction
% hn_mean=;
pri_sar_surf            = L1A_buffer(burst_index_nadir).pri_sar;
pulse_repetition_interval = int32(pri_sar_surf  * 1e12);
tm_h0=L1A_buffer(burst_index_nadir).h0_comp_sar_isp;
tm_cor2=L1A_buffer(burst_index_nadir).cor2_comp_sar_isp;

%----------G. Altimeter characterization variables --------------------------

%(TO DO!) Implement altimeter characterization variables
% residual_fixed_digital_gain_lrm
% residual_fixed_digital_gain_raw
% residual_fixed_digital_gain_rmc
% range_cor_external_group_delay_rx1
% range_cor_external_group_delay_rx2
% g_scaling_rx1
% g_scaling_rx2
% antenna_gain_ant1
% antenna_gain_ant2
% antenna_gain_ant

%------------------H. Flag variables --------------------------

%(TO DO!) Implement flag variables
% manoeuvre_flag
% surface_classification_flag

%------------------J. Waveform related variables --------------------------

% (TO DO!) Implement all the variables below once we have the right L1A chain and products
% burst_phase_array_cor_rx1
% burst_phase_array_cor_rx2
% burst_power_array_cor_rx1
% burst_power_array_cor_rx2
% instr_ext_phase_cor
% instr_int_phase_cor
% phase_slope_cor
% interf_base_vector
% pnr_estimation
% 
% ptr_main_lobe_width
% sig0_scaling_factor
% snr_instr_estimation

%----------------N. Look related variables --------------------------
% (TO DO!) Update using the l1a variable tm_source_seq_counter
look_counter= int32(L1BS.burst_index);
% (TO DO!) We have to write correctly the time stamp for each beam in the stack using the new L1A format/functions. For now we use the initial time of each burst
j=1;
for i=look_counter(1):look_counter(end)
   look_time(j) = L1A_buffer(i).time; % JPLZ: it does not work assigning the variaqbles inside the struct if we don't put it inside a loop
   j=j+1;
end
look_time=double(look_time');

look_i_samples_rx1_V = real(L1BS.beam_geo_corr);
look_q_samples_rx1_V = imag(L1BS.beam_geo_corr);
look_i_samples_rx2_V = real(L1BS.beam_geo_corr_2);
look_q_samples_rx2_V = imag(L1BS.beam_geo_corr_2);

    for i_surf=1:size(look_i_samples_rx1_V,1)
        i_scale_factor(i_surf) = single(max(abs(look_i_samples_rx1_V(i_surf,:))) / (2^7-1) * 1e6);
        q_scale_factor(i_surf) = single(max(abs(look_q_samples_rx1_V(i_surf,:))) / (2^7-1) * 1e6);
        iq_scale_factor_rx1(i_surf)= max([i_scale_factor(i_surf) q_scale_factor(i_surf)])';
        look_i_samples_rx1(i_surf,:) = int8(round(look_i_samples_rx1_V(i_surf,:) ./ iq_scale_factor_rx1(i_surf) .* 1e6));
        look_q_samples_rx1(i_surf,:) = int8(round(look_q_samples_rx1_V(i_surf,:) ./ iq_scale_factor_rx1(i_surf) .* 1e6));

        if strcmp(cnf.processing_mode,'SIN')
            i_scale_factor_2(i_surf) = single(max(abs(look_i_samples_rx1_V(i_surf,:))) / (2^7-1) * 1e6);
            q_scale_factor_2(i_surf) = single(max(abs(look_q_samples_rx1_V(i_surf,:))) / (2^7-1) * 1e6);
            iq_scale_factor_rx2(i_surf)= max([i_scale_factor_2(i_surf) q_scale_factor_2(i_surf)])';
            look_i_samples_rx2(i_surf,:) = int8(round(look_i_samples_rx2_V(i_surf,:) ./ iq_scale_factor_rx2(i_surf) .* 1e6));
            look_q_samples_rx2(i_surf,:) = int8(round(look_q_samples_rx2_V(i_surf,:) ./ iq_scale_factor_rx2(i_surf) .* 1e6));          
        end
    end
    
%     look_i_samples_rx1=look_i_samples_rx1.';
%     look_q_samples_rx1=look_q_samples_rx1.';
%     look_i_samples_rx2=look_i_samples_rx2.';
%     look_q_samples_rx2=look_q_samples_rx2.';

%     look_[iq]_samples_V(time, looks, samples) = look_[iq]_samples(time,looks, samples) * iq_scale_factor(time,looks)
%     
%     for n = 1:size(L1B.iq_scale_factor_rx1_ku, 1)
%               i_samples(n, :) = L1B.i_samples_rx1_ku(:,n,point_target_idx_ku).'.*L1B.iq_scale_factor_rx1_ku(n, point_target_idx_ku);
%               q_samples(n, :) = L1B.q_samples_rx1_ku(:,n,point_target_idx_ku).'.*L1B.iq_scale_factor_rx1_ku(n, point_target_idx_ku);
%           end
%         wfm_AC_rx1_ku = abs(i_samples + 1i*q_samples).^2;
%         
%         
% 
%               i_samples_V_nc = double(i_samples_nc(:,1:555,62)) .* iq_scale_factor_nc(1:555,62).';
%               q_samples_V_nc = double(q_samples_nc(:,1:555,62)) .* iq_scale_factor_nc(1:555,62).';
%        
%          wfm_AC_rx1_nc = abs(i_samples_V_nc + 1i*q_samples_V_nc).^2; 
    
%----------------O. Look characterization variables --------------------------

look_index = int8(L1BS.beam_index);
look_angle = int16(L1BS.look_ang_surf);
doppler_angle = int16(L1BS.doppler_ang_surf );
pointing_angle = int16(L1BS.pointing_ang_surf);
slant_range_correction_applied = int32(L1BS.slant_range_corr*1e3);
doppler_correction_applied = int32(L1BS.doppler_corr*1e3);
stack_mask = uint8(L1BS.stack_mask);

%----------R. Rate & collocation variables --------------------------

% (TO DO!) Implement these variables
%posting_rate
%ku_ka_collocation_flag

%----------T. Thermistor temperature variables --------------------------

% (TO DO!) Implement these variables
% thr0_txrf_ku
% thr1_txrf_ka
% thr2_rxrf_ku1
% thr3_rxrf_ku2
% thr4_rxrf_ka
% thr5_lo
% thr6_sspa_ku
% thr7_sspa_ka
% thr8_txnum
% thr9_rxnum_ku
% thr10_rxnum_ka
% thr11_dcdc_isps_1
% thr12_dcdc_isps_2
% thr13_formatter
% thr14_sequencer


% %Compute nadir burst for a given surface
% [~,beam_index_nadir]=min(abs(-pi/2+L1BS.beam_ang_surf'));
% 
% 
% burst_index_nadir= L1BS.burst_index(beam_index_nadir);
% 
% %source_seq_count_sar_isp_surf = L1A_buffer(burst_index_nadir).source_seq_count_sar_isp;
% ins_id=L1A_buffer(burst_index_nadir).ins_id;
% %ins_loop_stat=L1A_buffer(burst_index_nadir).ins_loop_stat;
% h0_comp_sar_isp=L1A_buffer(burst_index_nadir).h0_comp_sar_isp;
% cor2_comp_sar_isp=L1A_buffer(burst_index_nadir).cor2_comp_sar_isp;
% %ATT1_science=L1A_buffer(burst_index_nadir).ATT1_science;    
% %ATT2_science=L1A_buffer(burst_index_nadir).ATT2_science;     
% surface_type_flag = L1A_buffer(burst_index_nadir).surface_type_flag_bursts;
% %USO_correction = L1A_buffer(burst_index_nadir).USO_correction;
% %instrument_range_correction_tx_rx= L1A_buffer(burst_index_nadir).instrument_range_correction_tx_rx;
% %{
% if strcmp(cnf.mode,'SIN') && strcmp(mission,'CR2') && strcmp(cnf.processing_mode,'SIN')
%  instrument_range_correction=L1A_buffer(burst_index_nadir).instrument_range_correction_rx;
% end
% %}
% T0_sar_surf_nadir = L1BS.T0_sar_surf(beam_index_nadir);
% %      ATT1_science_corr=L1A_buffer(burst_index_nadir).ATT1_science_corr;  
% %      ATT2_science_corr=L1A_buffer(burst_index_nadir).ATT2_science_corr;  
% %      ATT1_delta_corr=L1A_buffer(burst_index_nadir).ATT1_delta_corr;
% %      ATT2_delta_corr=L1A_buffer(burst_index_nadir).ATT2_delta_corr;  
% 
% pri_sar_surf            = L1A_buffer(burst_index_nadir).pri_sar;
% 
% %geophysical corrections
% %{
% dry_tropo_correction=L1A_buffer(burst_index_nadir).dry_tropo_correction_bursts;
% wet_tropo_correction=L1A_buffer(burst_index_nadir).wet_tropo_correction_bursts;
% inverse_baro_correction=L1A_buffer(burst_index_nadir).inverse_baro_correction_bursts;
% Dynamic_atmospheric_correction=L1A_buffer(burst_index_nadir).Dynamic_atmospheric_correction_bursts;
% GIM_iono_correction=L1A_buffer(burst_index_nadir).GIM_iono_correction_bursts;
% cnf.model_iono_correction=L1A_buffer(burst_index_nadir).cnf.model_iono_correction_bursts;
% ocean_equilibrium_tide=L1A_buffer(burst_index_nadir).ocean_equilibrium_tide_bursts;
% long_period_tide_height=L1A_buffer(burst_index_nadir).long_period_tide_height_bursts;
% ocean_loading_tide=L1A_buffer(burst_index_nadir).ocean_loading_tide_bursts;
% solid_earth_tide=L1A_buffer(burst_index_nadir).solid_earth_tide_bursts;
% geocentric_polar_tide=L1A_buffer(burst_index_nadir).geocentric_polar_tide_bursts;
% %}
% %----------A. Time variables ----------------------------------------------
% % leap seconds in 2010 34s, 
% % +1 the 1st of July 2012 TAI_2012 = (12*365+4+181)*3600*24 + 35
% % +1 the 1st of July 2015 TAI_2015 = (15*365+4+181)*3600*24 + 36
% 
% TAI_2012 = (12*365+4+181)*3600*24 + 35;
% TAI_2015 = (15*365+4+181)*3600*24 + 36;
% time_l1b_echo = double(L1BS.time_surf);
% 
% % add leap seconds to the TAI time. Only valid for 
% if(time_l1b_echo < TAI_2012)
%     time_l1b_echo = time_l1b_echo - 34;
% elseif(time_l1b_echo > TAI_2015)
%     time_l1b_echo = time_l1b_echo - 36;    
% else
%     time_l1b_echo = time_l1b_echo - 35;
% end
%   
% UTC_day_l1b_echo = int16(floor(L1BS.time_surf./ cst.sec_in_day));
% UTC_sec_l1b_echo = double((time_l1b_echo - double(UTC_day_l1b_echo) * cst.sec_in_day));
% %isp_coarse_time_l1b_echo=uint32();
% %isp_fine_time_l1b_echo=int32();
% %sral_fine_time_l1b_echo=uint32();
% 
% 
% %tm_source_sequence_counter_ku = uint16(source_seq_count_sar_isp_surf);
% %l1b_record_counter_ku = uint16(0:(length(win_delay_surf)-1));
% 
% %----------B. Orbit and attitude variables ------------------------------
% lat_l1b_echo = int32(L1BS.lat_sat * 1e6);
% lon_l1b_echo = int32(L1BS.lon_sat * 1e6); 
% alt_l1b_echo = int32((L1BS.alt_sat-700000) * 1e4); 
% orb_alt_rate_l1b_echo = int16(L1BS.alt_rate_sat * 1e2);
% 
% satellite_mispointing_l1b_sar_echo_ku = int32([L1BS.pitch_surf * 180/cst.pi * 1e7; L1BS.roll_surf * 180/cst.pi * 1e7; L1BS.yaw_surf * 180/cst.pi * 1e7]);
% 
% %----------C. Flag time variables --------------------------------------
% 
% %----------D. Position/Velocity variables ------------------------------
% x_pos_l1b_echo=double(L1BS.x_sat');
% y_pos_l1b_echo=double(L1BS.y_sat');
% z_pos_l1b_echo=double(L1BS.z_sat');
% 
% x_vel_l1b_echo=double(L1BS.x_vel_sat');
% y_vel_l1b_echo=double(L1BS.y_vel_sat');
% z_vel_l1b_echo=double(L1BS.z_vel_sat');
% 
% %----------E. Navigation Bulletin --------------------------------------
% %seq_count_l1b_echo=uint16(source_seq_count_sar_isp_surf);
% 
% %----------F. Operating instrument & tracking --------------------------
% oper_instr_l1b_echo=int8(ins_id);%
% 
% %SAR_cnf.mode_l1b_echo=int8(ins_loop_stat);
% 
% 
% %---------- G. H0, COR2 and AGC ----------------------------------------
% % switch netcdf_type
% %     case 'netcdf4'
% %         h0_applied_l1b_echo = uint32(h0_comp_sar_isp/(3.125/64*1e-9));
% %     case 'netcdf3'
% %         h0_applied_l1b_echo = int64(h0_comp_sar_isp/(3.125/64*1e-9));
% % end
% h0_applied_l1b_echo = int64(h0_comp_sar_isp);        
% cor2_applied_l1b_echo=int16(cor2_comp_sar_isp/(3.125/1024*1e-9));
% %agccode_ku_l1b_echo=int8(-1.*(ATT1_science+ATT2_science));
% 
% 
% %-----------H. Surface Type flag----------------------------------------
% surf_type_l1b_echo =int8(surface_type_flag);
% 
% 
% %---------- I. Altimeter range and Corrections -------------------------
% range_ku_l1b_echo=int32((L1BS.win_delay_surf*cst.c/2-700000)*1e4);
% if cnf.include_wfms_aligned
%     %EM:14.04.2016
%     range_ku_surf_aligned_l1b_echo=int32((L1BS.win_delay_surf_aligned*cst.c/2-700000)*1e4);
% end
% %{
% uso_cor_l1b_echo=int32(USO_correction*cst.c/2*1e4);
% int_path_cor_ku_l1b_echo=int32(instrument_range_correction_tx_rx*1e4);
% if strcmp(cnf.mode,'SIN') && strcmp(mission,'CR2')&& strcmp(cnf.processing_mode,'SIN')
%     int_path_2_cor_ku_l1b_echo=int32(instrument_range_correction*1e4);
% end
% %}
% range_rate_l1b_echo=int32(T0_sar_surf_nadir*cst.c/2*1e3);
% win_delay_l1b_echo = double(L1BS.win_delay_surf*1e3);
% 
% %---------- J. AGC and Sigma0 scalings --------------------------------
% scale_factor_ku_l1b_echo=int32(L1B.wfm_scaling_factor*1e2);% sigma zero scaling factor
% 
% %---------- K. Stack characterization--------------------------------------
% % switch netcdf_type
% %     case 'netcdf4'
% %         nb_stack_l1b_echo=uint16(L1BS.N_beams_stack');
% %     case 'netcdf3'
% %         nb_stack_l1b_echo=int32(L1BS.N_beams_stack');
% % end
% nb_stack_l1b_echo=int32(L1B.N_beams_contributing);
% nb_stack_start_stop_l1b_echo=int32(L1B.N_beams_start_stop');
% %     aux=squeeze(L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,:));
% %     finite_indx=isfinite(aux);
% %     max_stack_l1b_echo = uint32(max(aux(finite_indx))*1e2);
% %     clear aux;
% if L1B.start_beam~=0
%     look_angle_start_l1b_echo = int16(L1BS.look_ang_surf(L1B.start_beam)'*1e6);
%     look_angle_stop_l1b_echo = int16(L1BS.look_ang_surf(L1B.stop_beam)*1e6);
%     doppler_angle_start_l1b_echo = int16((L1BS.doppler_ang_surf(L1B.start_beam))'*1e6);
%     doppler_angle_stop_l1b_echo  = int16((L1BS.doppler_ang_surf(L1B.stop_beam))'*1e6);
%     pointing_angle_start_l1b_echo = int16(L1BS.pointing_ang_surf(L1B.start_beam)*1e6)';
%     pointing_angle_stop_l1b_echo  = int16(L1BS.pointing_ang_surf(L1B.stop_beam)'*1e6);
% else
%     look_angle_start_l1b_echo = int16(0);
%     look_angle_stop_l1b_echo = int16(0);
%     doppler_angle_start_l1b_echo = int16(0);
%     doppler_angle_stop_l1b_echo  = int16(0);
%     pointing_angle_start_l1b_echo = int16(0);
%     pointing_angle_stop_l1b_echo  = int16(0);
% end
% 
% % switch netcdf_type
% %     case 'netcdf4'
% %         stdev_stack_l1b_echo = uint32(L1B.stack_std' * 1e6);
% %     case 'netcdf3'
% %         stdev_stack_l1b_echo = int64(L1B.stack_std' * 1e6);
% % end
% if cnf.compute_L1_stadistics
%     stdev_stack_l1b_echo = int64(L1B.stack_std' * 1e6);
%     skew_stack_l1b_echo = int32(L1B.stack_skewness' * 1e6);
%     kurt_stack_l1b_echo = int32(L1B.stack_kurtosis' * 1e6);
%     gaussian_fitting_centre_look_l1b_echo = int16(L1B.stack_look_ang_centre' * 1e6);
%     gaussian_fitting_centre_pointing_l1b_echo = int16(L1B.stack_pointing_ang_centre' * 1e6);    
% end
% 
% %------------- L. Altimeter engineering variables -------------------------
% altimeter_clock_l1b_echo = int32(1./T0_sar_surf_nadir - chd.bw)* 1e9;
% pri_l1b_echo = int64(pri_sar_surf  * 1e12);
% % switch netcdf_type
% %     case 'netcdf4'
% %         pri_lrm_l1b_echo = uint32(pri_sar_isp_surf .* T0_sar_surf_nadir * 1e19);
% %     case 'netcdf3'
% %         pri_lrm_l1b_echo = int64(pri_sar_isp_surf .* T0_sar_surf_nadir * 1e19);
% % end
% 
% 
% %------------- M. Waveform related variables -----------------------------
% waveform_scale_factor_l1b_echo = single(max(L1B.wfm_cor_i2q2)/ (2^16-2));
% i2q2_meas_ku_l1b_echo = int32(round(L1B.wfm_cor_i2q2./ waveform_scale_factor_l1b_echo));
% 
% if strcmp(cnf.processing_mode,'SIN')
%  %  waveform_scale_factor_l1b_echo_2 = single(max(L1B.wfm_cor_i2q2_2)/ (2^16-2));
% %   i2q2_meas_ku_l1b_echo_2 = int32(round(L1B.wfm_cor_i2q2_2./ waveform_scale_factor_l1b_echo_2));
%    phase_diff_meas_ku_l1b_echo = double(L1B.phase_difference);
%    coherence_meas_ku_l1b_echo = double(L1B.coherence);
% end
% 
% if cnf.include_wfms_aligned
%     i2q2_meas_ku_wdcorr_l1b_echo = int32(round(L1B.wfm_cor_i2q2_wdcorr./ waveform_scale_factor_l1b_echo));
% end
% switch netcdf_type
%     case 'netcdf4'
%         i2q2_meas_ku_l1b_echo = uint16(round(L1B.wfm_cor_i2q2.*10^0.3 ./ waveform_scale_factor_l1b_echo));
%         i2q2_meas_ku_wdcorr_l1b_echo = uint16(round(L1B.wfm_cor_i2q2_wdcorr.*10^0.3 ./ waveform_scale_factor_l1b_echo));
%     case 'netcdf3'
%         i2q2_meas_ku_l1b_echo = int32(round(L1B.wfm_cor_i2q2.*10^0.3 ./ waveform_scale_factor_l1b_echo));
%         i2q2_meas_ku_wdcorr_l1b_echo = int32(round(L1B.wfm_cor_i2q2_wdcorr.*10^0.3 ./ waveform_scale_factor_l1b_echo));
% end

% EM: 22.09.2016
%{
if strcmp(cnf.mode,'SIN') && strcmp(mission,'CR2')&& strcmp(cnf.processing_mode,'SIN')
    waveform_scale_factor_l1b_echo_2 = single(max(L1B.wfm_cor_i2q2_2)/ (2^16-2));
    i2q2_meas_ku_l1b_echo_2 = int32(round(L1B.wfm_cor_i2q2_2./ waveform_scale_factor_l1b_echo_2));
    
    coherence_meas_ku_l1b_echo=int16(L1B.coherence*1e3);
    
    phase_difference_meas_ku_l1b_echo=int32(L1B.phase_difference*1e6);
          
    
    if cnf.include_wfms_aligned
        
        i2q2_meas_ku_wdcorr_l1b_echo_2 = int32(round(L1B.wfm_cor_i2q2_wdcorr_2./ waveform_scale_factor_l1b_echo_2));
        
        coherence_meas_ku_wdcorr_l1b_echo= int16(L1B.coherence_wdcorr*1e3);
        
        phase_difference_meas_ku_wdcorr_l1b_echo=int32(L1B.phase_difference_wdcorr*1e6);  

    end
           
end
%}

%------------- N. Geophysical Corrections variables ---------------------
%{
dry_tropo_correction=int32(dry_tropo_correction.*1e3);
wet_tropo_correction=int32(wet_tropo_correction.*1e3);
inverse_baro_correction=int32(inverse_baro_correction.*1e3);
Dynamic_atmospheric_correction=int32(Dynamic_atmospheric_correction.*1e3);
GIM_iono_correction=int32(GIM_iono_correction.*1e3);
cnf.model_iono_correction=int32(cnf.model_iono_correction.*1e3);
ocean_equilibrium_tide=int32(ocean_equilibrium_tide.*1e3);
long_period_tide_height=int32(long_period_tide_height.*1e3);
ocean_loading_tide=int32(ocean_loading_tide.*1e3);
solid_earth_tide=int32(solid_earth_tide.*1e3);
geocentric_polar_tide=int32(geocentric_polar_tide.*1e3);
%}

%added EM 04.10.2016
% %------------- P. ACDC variables ----------------------------------------
%{
if ACDC_application_cnf
    %-------------- Waveform & Scaling ------------------------------------   
    %ACDC waveform
    waveform_scale_factor_ACDC_echo = single(max(L1B.ACDC.waveform)/ (2^16-2));
    i2q2_meas_ku_ACDC_echo = int32(round(L1B.ACDC.waveform./ waveform_scale_factor_ACDC_echo));
    %ACDC fitted waveform
    waveform_scale_factor_fit_ACDC_echo = single(max(L1B.ACDC.ml_wav_fitted_ACDC)/ (2^16-2));
    i2q2_fit_ku_ACDC_echo = int32(round(L1B.ACDC.ml_wav_fitted_ACDC./ waveform_scale_factor_fit_ACDC_echo));
    
    %L1B fitted waveform
    waveform_scale_factor_fit_echo = single(max(L1B.ACDC.ml_wav_conv)/ (2^16-2));
    i2q2_fit_ku_echo = int32(round(L1B.ACDC.ml_wav_conv./ waveform_scale_factor_fit_echo));

    %-------------- Range index -------------------------------------------
    range_index_ACDC_echo=int32(L1B.ACDC.range_index.*1e4);
    
    %-------------- Geophysical retrievals --------------------------------
    %******************* retracked range **********************************   
    retracked_range_ACDC_20_ku=int32((L1B.ACDC.tracker_range-700000).*1e4);
    %******************* epoch ********************************************
    epoch_ACDC_20_ku=int32(L1B.ACDC.retracking_cor.*2/cst.c.*1e15);    
    %******************** SWH  ********************************************
    swh_ACDC_20_ku=int16(L1B.ACDC.Hs.*1e3);
    swh_ACDC_PRE_20_ku=int16(L1B.ACDC.Hs_conv.*1e3);
    %******************** SSH  ********************************************
    ssh_ACDC_20_ku=int32(L1B.ACDC.SSH.*1e4);
    ssh_ACDC_PRE_20_ku=int32(L1B.ACDC.SSH_conv.*1e4);
    %******************** SIGMA0 ******************************************    
    sig0_ACDC_20_ku=int16(L1B.ACDC.sigma0.*1e2);
    sig0_ACDC_PRE_20_ku=int16(L1B.ACDC.sigma0_conv.*1e2);
    %********************* Pu: fitted peak power **************************    
    Pu_analytical_ACDC_20_ku=int16(L1B.ACDC.Pu.*1e2);
    %********************** Pearson Correlation Coefficient ***************   
    Pearson_corr_ACDC_20_ku=int16(L1B.ACDC.corr_coeff.*1e2);
    Pearson_corr_ACDC_PRE_20_ku=int16(L1B.ACDC.corr_coeff_conv.*1e2);
    %************************ Fitting flag ********************************    
    Flag_fitting_ACDC_20_ku=int8(L1B.ACDC.flag_fitting);
      
end
%}

%%----------0. Define groups ----------------------------------------------
global_ncid = netcdf.inqNcid(files.ncid_L1Bs,'global');
global_ku_ncid = netcdf.inqNcid(global_ncid,'ku');
global_ka_ncid = netcdf.inqNcid(global_ncid,'ka');

data_ncid = netcdf.inqNcid(files.ncid_L1Bs,'data');
data_ku_ncid = netcdf.inqNcid(data_ncid,'ku');
data_ka_ncid = netcdf.inqNcid(data_ncid,'ka');

if strcmp(chd.band,'Ku')
    global_currentband_ncid=global_ku_ncid;
    data_currentband_ncid=data_ku_ncid;
elseif strcmp(chd.band,'Ka')
    global_currentband_ncid=global_ka_ncid;
    data_currentband_ncid=data_ka_ncid;
end

%% PACKING L1Bs
%----------A. Time & counters variables ----------------------------------------------
% Write the dimension variables only for the first record
if i_surf_stacked==0
    var_id = netcdf.inqVarID(files.ncid_L1Bs,'looks');
    netcdf.putVar(files.ncid_L1Bs,var_id,0, nl_dimension_size, 0:nl_dimension_size-1);
    
    var_id = netcdf.inqVarID(files.ncid_L1Bs,'samples_ov');
    netcdf.putVar(files.ncid_L1Bs,var_id,0 ,ns_dimension_size, 0:ns_dimension_size-1);
    
    var_id = netcdf.inqVarID(files.ncid_L1Bs,'space_3d');
    netcdf.putVar(files.ncid_L1Bs,var_id,0, space_3D_dimension_size, 0:space_3D_dimension_size-1);
    
    var_id = netcdf.inqVarID(files.ncid_L1Bs,'pulses');
    netcdf.putVar(files.ncid_L1Bs,var_id,0 ,np_dimension_size, 0:np_dimension_size-1);
end

var_id = netcdf.inqVarID(data_currentband_ncid,'l1b_record_counter');
netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,l1b_record_counter);

var_id = netcdf.inqVarID(data_currentband_ncid,'time');
netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,time);

% var_id = netcdf.inqVarID(data_currentband_ncid,'time_tai');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,time_tai);

% var_id = netcdf.inqVarID(data_currentband_ncid,'tm_source_sequence_counter');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,tm_source_sequence_counter);

%----------B. Orbit and attitude variables ------------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'altitude');
netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,altitude);

var_id = netcdf.inqVarID(data_currentband_ncid,'altitude_rate');
netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,altitude_rate);

var_id = netcdf.inqVarID(data_currentband_ncid,'latitude');
netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,latitude);

var_id = netcdf.inqVarID(data_currentband_ncid,'longitude');
netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,longitude);

% if strcmp(chd.band,'Ku')
%     var_id = netcdf.inqVarID(data_ku_ncid,'off_nadir_roll_angle_ant1');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,off_nadir_roll_angle_ant1);
%     
%     var_id = netcdf.inqVarID(data_ku_ncid,'off_nadir_pitch_angle_ant1');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,off_nadir_pitch_angle_ant1);
%     
%     var_id = netcdf.inqVarID(data_ku_ncid,'off_nadir_yaw_angle_ant1_name');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,off_nadir_yaw_angle_ant1);
% end
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(data_ku_ncid,'off_nadir_roll_angle_ant2');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,off_nadir_roll_angle_ant2);
%     
%     var_id = netcdf.inqVarID(data_ku_ncid,'off_nadir_pitch_angle_ant2');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,off_nadir_pitch_angle_ant2);
%     
%     var_id = netcdf.inqVarID(data_ku_ncid,'off_nadir_yaw_angle_ant2');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,off_nadir_yaw_angle_ant2);
% end
% 
% if strcmp(chd.band,'Ka')
%     var_id = netcdf.inqVarID(data_ka_ncid,'off_nadir_roll_angle');
%     netcdf.putVar(data_ka_ncid,var_id,i_surf_stacked,off_nadir_roll_angle);
%     
%     var_id = netcdf.inqVarID(data_ka_ncid,'off_nadir_pitch_angle');
%     netcdf.putVar(data_ka_ncid,var_id,i_surf_stacked,off_nadir_pitch_angle);
%     
%     var_id = netcdf.inqVarID(data_ka_ncid,'off_nadir_yaw_angle');
%     netcdf.putVar(data_ka_ncid,var_id,i_surf_stacked,off_nadir_yaw_angle);
% end
% 
% var_id = netcdf.inqVarID(global_currentband_ncid,'orbit_type_flag');
% netcdf.putVar(global_currentband_ncid,var_id,i_surf_stacked,orbit_type_flag);

var_id = netcdf.inqVarID(data_currentband_ncid,'position_vector');
netcdf.putVar(data_currentband_ncid,var_id,[0 i_surf_stacked],[space_3D_dimension_size 1], position_vector);

var_id = netcdf.inqVarID(data_currentband_ncid,'velocity_vector');
netcdf.putVar(data_currentband_ncid,var_id,[0 i_surf_stacked],[space_3D_dimension_size 1], velocity_vector);

%----------C. Configuration and quality variables --------------------------------------
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'mcd_flags');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,mcd_flags);
% 
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'processing_configuration_flags');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,processing_configuration_flags);
% 
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'hr_mcd_flags');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,hr_mcd_flags);
% 
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'hr_processing_configuration_flags');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,hr_processing_configuration_flags);
% 
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'iris_instrument_configuration_flags');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,iris_instrument_configuration_flags);
% 
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'iris_mode_flag');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,iris_mode_flag);
% 
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'range_oversampling_factor');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,range_oversampling_factor);
% 
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'telemetry_type_flag');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,telemetry_type_flag);

%----------------D. Altimeter range variables ---------------------------
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'range_cor_com_ant1');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,range_cor_com_ant1);
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(files.ncid_L1Bs,'range_cor_com_ant2');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,range_cor_com_ant2);
% end
% 
% if strcmp(chd.band,'Ka')
%     var_id = netcdf.inqVarID(files.ncid_L1Bs,'range_cor_com');
%     netcdf.putVar(data_ka_ncid,var_id,i_surf_stacked,range_cor_com);
% end
% 
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'range_cor_doppler');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,range_cor_doppler);
% 
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'range_cor_internal_delay_cal_rx1');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,range_cor_internal_delay_cal_rx1);
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(files.ncid_L1Bs,'range_cor_internal_delay_cal_rx2');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,range_cor_internal_delay_cal_rx2);
% end
% 
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'range_cor_internal_delay_att_rx1');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,range_cor_internal_delay_att_rx1);
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(files.ncid_L1Bs,'range_cor_internal_delay_cal_rx2');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,range_cor_internal_delay_cal_rx2);
% end
% 
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'range_cor_uso');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,range_cor_uso);
% 
% var_id = netcdf.inqVarID(files.ncid_L1Bs,'range_cor_reference_sample');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,range_cor_reference_sample);

var_id = netcdf.inqVarID(data_currentband_ncid,'tracker_range_calibrated');
netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,tracker_range_calibrated);

% var_id = netcdf.inqVarID(files.ncid_L1Bs,'tracker_range_calibrated_gnss');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,tracker_range_calibrated_gnss);

%----------E. Altimeter power variables ------------------------------
% var_id = netcdf.inqVarID(data_currentband_ncid,'altimeter_power_drift_rx1');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,altimeter_power_drift_rx1);
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(data_ku_ncid,'altimeter_power_drift_rx2');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,altimeter_power_drift_rx2);
% end
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'attenuation_calibrated_rx1');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,attenuation_calibrated_rx1);
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(data_ku_ncid,'attenuation_calibrated_rx2');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,attenuation_calibrated_rx2);
% end
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'power_scaling_to_antenna_rx1');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,power_scaling_to_antenna_rx1);
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(data_ku_ncid,'power_scaling_to_antenna_rx2');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,power_scaling_to_antenna_rx2);
% end
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'variable_digital_gain');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,variable_digital_gain);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'cal1_power_rx1');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,cal1_power_rx1);
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(data_ku_ncid,'cal1_power_rx2');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,cal1_power_rx2);
% end

%----------F. Altimeter engineering variables --------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'altimeter_clock');
netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,altimeter_clock);

% var_id = netcdf.inqVarID(files.ncid_L1Bs,'hn_mean');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,hn_mean);

var_id = netcdf.inqVarID(data_currentband_ncid,'pulse_repetition_interval');
netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,pulse_repetition_interval);

var_id = netcdf.inqVarID(data_currentband_ncid,'tm_cor2');
netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,tm_cor2);

var_id = netcdf.inqVarID(data_currentband_ncid,'tm_h0');
netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,tm_h0);

%----------G. Altimeter characterization variables --------------------------
% var_id = netcdf.inqVarID(global_currentband_ncid,'residual_fixed_digital_gain_raw');
% netcdf.putVar(global_currentband_ncid,var_id,i_surf_stacked,residual_fixed_digital_gain_raw);
% 
% var_id = netcdf.inqVarID(global_currentband_ncid,'residual_fixed_digital_gain_rmc');
% netcdf.putVar(global_currentband_ncid,var_id,i_surf_stacked,residual_fixed_digital_gain_rmc);
% 
% var_id = netcdf.inqVarID(global_currentband_ncid,'range_cor_external_group_delay_rx1');
% netcdf.putVar(global_currentband_ncid,var_id,i_surf_stacked,range_cor_external_group_delay_rx1);
% 
% var_id = netcdf.inqVarID(global_currentband_ncid,'range_cor_external_group_delay_rx2');
% netcdf.putVar(global_currentband_ncid,var_id,i_surf_stacked,range_cor_external_group_delay_rx2);
% 
% var_id = netcdf.inqVarID(global_currentband_ncid,'g_scaling_rx1');
% netcdf.putVar(global_currentband_ncid,var_id,i_surf_stacked,g_scaling_rx1);
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(global_ku_ncid,'g_scaling_rx2');
%     netcdf.putVar(global_ku_ncid,var_id,i_surf_stacked,g_scaling_rx2);
% end
% 
% var_id = netcdf.inqVarID(global_currentband_ncid,'antenna_gain_ant1');
% netcdf.putVar(global_currentband_ncid,var_id,i_surf_stacked,antenna_gain_ant1);
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(global_ku_ncid,'antenna_gain_ant2');
%     netcdf.putVar(global_ku_ncid,var_id,i_surf_stacked,antenna_gain_ant2);
% end
% 
% if strcmp(chd.band,'Ka')
%     var_id = netcdf.inqVarID(global_ka_ncid,'antenna_gain_ant');
%     netcdf.putVar(global_ka_ncid,var_id,i_surf_stacked,antenna_gain_ant);
% end

%------------------H. Flag variables --------------------------
% var_id = netcdf.inqVarID(data_currentband_ncid,'manoeuvre_flag');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,manoeuvre_flag);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'surface_classification_flag');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,surface_classification_flag);

%------------------J. Waveform related variables --------------------------
% var_id = netcdf.inqVarID(global_currentband_ncid,'burst_phase_array_cor_rx1');
% netcdf.putVar(global_currentband_ncid,var_id,i_surf_stacked,burst_phase_array_cor_rx1);
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(global_ku_ncid,'burst_phase_array_cor_rx2');
% netcdf.putVar(global_ku_ncid,var_id,i_surf_stacked,burst_phase_array_cor_rx2);
% end
% 
% var_id = netcdf.inqVarID(global_currentband_ncid,'burst_power_array_cor_rx1');
% netcdf.putVar(global_currentband_ncid,var_id,i_surf_stacked,burst_power_array_cor_rx1);
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(global_ku_ncid,'burst_power_array_cor_rx2');
%     netcdf.putVar(global_ku_ncid,var_id,i_surf_stacked,burst_power_array_cor_rx2);
% end
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(global_ku_ncid,'instr_ext_phase_cor');
%     netcdf.putVar(global_ku_ncid,var_id,i_surf_stacked,instr_ext_phase_cor);
% end
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(global_ku_ncid,'instr_int_phase_cor');
%     netcdf.putVar(global_ku_ncid,var_id,i_surf_stacked,instr_ext_phase_cor);
% end
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(data_ku_ncid,'phase_slope_cor');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,phase_slope_cor);
% end
% 
% if strcmp(cnf.processing_mode,'SIN')
%     var_id = netcdf.inqVarID(data_ku_ncid,'interf_base_vector');
%     netcdf.putVar(data_ku_ncid,var_id,i_surf_stacked,interf_base_vector);
% end
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'pnr_estimation_name');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,pnr_estimation_name);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'ptr_main_lobe_width_name');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,ptr_main_lobe_width_name); 
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'sig0_scaling_factor');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,sig0_scaling_factor); 
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'snr_instr_estimation');
% netcdf.putVar(data_currentband_ncid,var_id,i_surf_stacked,snr_instr_estimation);

%----------------N. Look related variables --------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'look_counter');
netcdf.putVar(data_currentband_ncid,var_id,[0 i_surf_stacked],[size(look_counter,1) 1],look_counter);

var_id = netcdf.inqVarID(data_currentband_ncid,'look_time');
netcdf.putVar(data_currentband_ncid,var_id,[0 i_surf_stacked],[size(look_time,1) 1],look_time);

var_id = netcdf.inqVarID(data_currentband_ncid,'look_i_samples_rx1');
netcdf.putVar(data_currentband_ncid,var_id,[0 0 i_surf_stacked],[ns_dimension_size size(look_i_samples_rx1,1) 1],look_i_samples_rx1.');

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'look_i_samples_rx2');
    netcdf.putVar(data_ku_ncid,var_id,[0 0 i_surf_stacked],[ns_dimension_size size(look_i_samples_rx2,1) 1],look_i_samples_rx2.');
end

var_id = netcdf.inqVarID(data_currentband_ncid,'look_q_samples_rx1');
netcdf.putVar(data_currentband_ncid,var_id,[0 0 i_surf_stacked],[ns_dimension_size size(look_q_samples_rx1,1) 1],look_q_samples_rx1.');

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'look_q_samples_rx2');
    netcdf.putVar(data_ku_ncid,var_id,[0 0 i_surf_stacked],[ns_dimension_size size(look_q_samples_rx2,1) 1],look_q_samples_rx2.');
end

var_id = netcdf.inqVarID(data_currentband_ncid,'iq_scale_factor_rx1');
netcdf.putVar(data_currentband_ncid,var_id,[0 i_surf_stacked],[size(iq_scale_factor_rx1,2) 1],iq_scale_factor_rx1);

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(data_ku_ncid,'iq_scale_factor_rx2');
    netcdf.putVar(data_ku_ncid ,var_id,[0 i_surf_stacked],[size(iq_scale_factor_rx1,2) 1],iq_scale_factor_rx2);
end

%----------------O. Look characterization variables --------------------------
var_id = netcdf.inqVarID(data_currentband_ncid,'look_index');
netcdf.putVar(data_currentband_ncid ,var_id,[0 i_surf_stacked],[size(look_index,1) 1],look_index);

var_id = netcdf.inqVarID(data_currentband_ncid,'look_angle');
netcdf.putVar(data_currentband_ncid ,var_id,[0 i_surf_stacked],[size(look_angle,2) 1],look_angle);

var_id = netcdf.inqVarID(data_currentband_ncid,'doppler_angle');
netcdf.putVar(data_currentband_ncid ,var_id,[0 i_surf_stacked],[size(doppler_angle,2) 1],doppler_angle);

var_id = netcdf.inqVarID(data_currentband_ncid,'pointing_angle');
netcdf.putVar(data_currentband_ncid ,var_id,[0 i_surf_stacked],[size(pointing_angle,2) 1],pointing_angle);

var_id = netcdf.inqVarID(data_currentband_ncid,'slant_range_correction_applied');
netcdf.putVar(data_currentband_ncid ,var_id,[0 i_surf_stacked],[size(slant_range_correction_applied,2) 1],slant_range_correction_applied);

var_id = netcdf.inqVarID(data_currentband_ncid,'doppler_correction_applied');
netcdf.putVar(data_currentband_ncid ,var_id,[0 i_surf_stacked],[size(doppler_correction_applied,2) 1],doppler_correction_applied);

% (TO DO!) Check how to write the mask inside this variable. If I'm not mistaken according to PFS we have to write in dimensions [looks X time], while in the GPP we have [Range samples X looks X time]
% var_id = netcdf.inqVarID(data_currentband_ncid,'stack_mask');
% netcdf.putVar(data_currentband_ncid ,var_id,[0 i_surf_stacked],[size(stack_mask,1) 1],stack_mask);

%----------R. Rate & collocation variables --------------------------
% var_id = netcdf.inqVarID(global_currentband_ncid,'posting_rate');
% netcdf.putVar(global_currentband_ncid,var_id,i_surf_stacked,posting_rate);
% 
% var_id = netcdf.inqVarID(global_currentband_ncid,'ku_ka_collocation_flag');
% netcdf.putVar(global_currentband_ncid ,var_id,i_surf_stacked,ku_ka_collocation_flag);

%----------T. Thermistor temperature variables --------------------------
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr0_txrf_ku');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr0_txrf_ku);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr1_txrf_ka');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr1_txrf_ka);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr2_rxrf_ku1');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr2_rxrf_ku1);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr3_rxrf_ku2');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr3_rxrf_ku2);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr4_rxrf_ka');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr4_rxrf_ka);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr5_lo');
% netcdf.putVar(data_currentband_ncidd ,var_id,i_surf_stacked,thr5_lo);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr6_sspa_ku');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr6_sspa_ku);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr7_sspa_ka');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr7_sspa_ka);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr8_txnum');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr8_txnum);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr9_rxnum_ku');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr9_rxnum_ku);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr10_rxnum_ka');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr10_rxnum_ka);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr11_dcdc_isps_1');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr11_dcdc_isps_1);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr12_dcdc_isps_2');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr12_dcdc_isps_2);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr13_formatter');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr13_formatter);
% 
% var_id = netcdf.inqVarID(data_currentband_ncid,'thr14_sequencer');
% netcdf.putVar(data_currentband_ncid ,var_id,i_surf_stacked,thr14_sequencer);


% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);
% netcdf.close(files.ncid_L1Bs);