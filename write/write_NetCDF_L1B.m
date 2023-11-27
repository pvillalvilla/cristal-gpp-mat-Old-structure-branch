%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the CODING & PACKING 
% algorithm for Level-1B product using Sentinel-3 like format
% Ref: Product Data Format Specification - SRAL/MWR Level 1 & 2 Instrument
% Products issue 2.0
%
% ---------------------------------------------------------
% Objective: Pack variables and write the NETCDF
% 
% INPUTs : Workspace
% OUTPUTs: TM Structure as defined on isardSAT_JasonCS_DPM
%
% ----------------------------------------------------------
% Author:    Eduard Makhoul/ isardSAT
%            Gorka Moyano  / isardSAT
%            Roger Escola  / isardSAT
%            Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT ()
% 
% Versions
% 1.0 
% 1.1 Updated time conversion for data between 2010 and 2016 (CR2)
% 2.0 Transformed to a function. Writting one record
% 2.1 changed chd.N_samples_sar_sar_chd by chd.N_samples_sar
% 2.2 (04.10.2016, EM) Adding ACDC output related variables
% 2.3 nb_stack_l1b_echo is L1B.N_beams_contributing instead of L1BS.N_beams_stack
% 2.4 added cnf flag cnf.processing_mode
% 2.5 changed int32 for doubles and int64 for int32
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils

function write_NetCDF_L1B(files,L1A_buffer,L1BS, L1B,i_surf_stacked,cnf,chd,cst)
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
% nl_dimension_size = N_max_beams_stack_chd;
% ns_dimension = 'Ns';
% space_3D_dimension = 'space_3D';
% space_3D_dimension_size = 3;



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

%Compute nadir burst for a given surface
[~,beam_index_nadir]=min(abs(-pi/2+L1BS.beam_ang_surf'));


burst_index_nadir= L1BS.burst_index(beam_index_nadir);

%source_seq_count_sar_isp_surf = L1A_buffer(burst_index_nadir).source_seq_count_sar_isp;
ins_id=L1A_buffer(burst_index_nadir).ins_id;
%ins_loop_stat=L1A_buffer(burst_index_nadir).ins_loop_stat;
h0_comp_sar_isp=L1A_buffer(burst_index_nadir).h0_comp_sar_isp;
cor2_comp_sar_isp=L1A_buffer(burst_index_nadir).cor2_comp_sar_isp;
%ATT1_science=L1A_buffer(burst_index_nadir).ATT1_science;    
%ATT2_science=L1A_buffer(burst_index_nadir).ATT2_science;     
surface_type_flag = L1A_buffer(burst_index_nadir).surface_type_flag_bursts;
%USO_correction = L1A_buffer(burst_index_nadir).USO_correction;
%instrument_range_correction_tx_rx= L1A_buffer(burst_index_nadir).instrument_range_correction_tx_rx;
%{
if strcmp(cnf.mode,'SIN') && strcmp(mission,'CR2') && strcmp(cnf.processing_mode,'SIN')
 instrument_range_correction=L1A_buffer(burst_index_nadir).instrument_range_correction_rx;
end
%}
T0_sar_surf_nadir = L1BS.T0_sar_surf(beam_index_nadir);
%      ATT1_science_corr=L1A_buffer(burst_index_nadir).ATT1_science_corr;  
%      ATT2_science_corr=L1A_buffer(burst_index_nadir).ATT2_science_corr;  
%      ATT1_delta_corr=L1A_buffer(burst_index_nadir).ATT1_delta_corr;
%      ATT2_delta_corr=L1A_buffer(burst_index_nadir).ATT2_delta_corr;  

pri_sar_surf            = L1A_buffer(burst_index_nadir).pri_sar;

%geophysical corrections
%{
dry_tropo_correction=L1A_buffer(burst_index_nadir).dry_tropo_correction_bursts;
wet_tropo_correction=L1A_buffer(burst_index_nadir).wet_tropo_correction_bursts;
inverse_baro_correction=L1A_buffer(burst_index_nadir).inverse_baro_correction_bursts;
Dynamic_atmospheric_correction=L1A_buffer(burst_index_nadir).Dynamic_atmospheric_correction_bursts;
GIM_iono_correction=L1A_buffer(burst_index_nadir).GIM_iono_correction_bursts;
cnf.model_iono_correction=L1A_buffer(burst_index_nadir).cnf.model_iono_correction_bursts;
ocean_equilibrium_tide=L1A_buffer(burst_index_nadir).ocean_equilibrium_tide_bursts;
long_period_tide_height=L1A_buffer(burst_index_nadir).long_period_tide_height_bursts;
ocean_loading_tide=L1A_buffer(burst_index_nadir).ocean_loading_tide_bursts;
solid_earth_tide=L1A_buffer(burst_index_nadir).solid_earth_tide_bursts;
geocentric_polar_tide=L1A_buffer(burst_index_nadir).geocentric_polar_tide_bursts;
%}
%----------A. Time variables ----------------------------------------------
% leap seconds in 2010 34s, 
% +1 the 1st of July 2012 TAI_2012 = (12*365+4+181)*3600*24 + 35
% +1 the 1st of July 2015 TAI_2015 = (15*365+4+181)*3600*24 + 36

TAI_2012 = (12*365+4+181)*3600*24 + 35;
TAI_2015 = (15*365+4+181)*3600*24 + 36;
time_l1b_echo = double(L1BS.time_surf);

% add leap seconds to the TAI time. Only valid for 
if(time_l1b_echo < TAI_2012)
    time_l1b_echo = time_l1b_echo - 34;
elseif(time_l1b_echo > TAI_2015)
    time_l1b_echo = time_l1b_echo - 36;    
else
    time_l1b_echo = time_l1b_echo - 35;
end
  
UTC_day_l1b_echo = int16(floor(L1BS.time_surf./ cst.sec_in_day));
UTC_sec_l1b_echo = double((time_l1b_echo - double(UTC_day_l1b_echo) * cst.sec_in_day));
%isp_coarse_time_l1b_echo=uint32();
%isp_fine_time_l1b_echo=int32();
%sral_fine_time_l1b_echo=uint32();


%tm_source_sequence_counter_ku = uint16(source_seq_count_sar_isp_surf);
%l1b_record_counter_ku = uint16(0:(length(win_delay_surf)-1));

%----------B. Orbit and attitude variables ------------------------------
lat_l1b_echo = int32(L1BS.lat_sat * 1e6);
lon_l1b_echo = int32(L1BS.lon_sat * 1e6); 
alt_l1b_echo = int32((L1BS.alt_sat-700000) * 1e4); 
orb_alt_rate_l1b_echo = int16(L1BS.alt_rate_sat * 1e2);

satellite_mispointing_l1b_sar_echo_ku = int32([L1BS.pitch_surf * 180/cst.pi * 1e7; L1BS.roll_surf * 180/cst.pi * 1e7; L1BS.yaw_surf * 180/cst.pi * 1e7]);

%----------C. Flag time variables --------------------------------------

%----------D. Position/Velocity variables ------------------------------
x_pos_l1b_echo=double(L1BS.x_sat');
y_pos_l1b_echo=double(L1BS.y_sat');
z_pos_l1b_echo=double(L1BS.z_sat');

x_vel_l1b_echo=double(L1BS.x_vel_sat');
y_vel_l1b_echo=double(L1BS.y_vel_sat');
z_vel_l1b_echo=double(L1BS.z_vel_sat');

%----------E. Navigation Bulletin --------------------------------------
%seq_count_l1b_echo=uint16(source_seq_count_sar_isp_surf);

%----------F. Operating instrument & tracking --------------------------
oper_instr_l1b_echo=int8(ins_id);%

%SAR_cnf.mode_l1b_echo=int8(ins_loop_stat);


%---------- G. H0, COR2 and AGC ----------------------------------------
% switch netcdf_type
%     case 'netcdf4'
%         h0_applied_l1b_echo = uint32(h0_comp_sar_isp/(3.125/64*1e-9));
%     case 'netcdf3'
%         h0_applied_l1b_echo = int64(h0_comp_sar_isp/(3.125/64*1e-9));
% end
h0_applied_l1b_echo = int64(h0_comp_sar_isp);        
cor2_applied_l1b_echo=int16(cor2_comp_sar_isp/(3.125/1024*1e-9));
%agccode_ku_l1b_echo=int8(-1.*(ATT1_science+ATT2_science));


%-----------H. Surface Type flag----------------------------------------
surf_type_l1b_echo =int8(surface_type_flag);


%---------- I. Altimeter range and Corrections -------------------------
range_ku_l1b_echo=int32((L1BS.win_delay_surf*cst.c/2-700000)*1e4);
if cnf.include_wfms_aligned
    %EM:14.04.2016
    range_ku_surf_aligned_l1b_echo=int32((L1BS.win_delay_surf_aligned*cst.c/2-700000)*1e4);
end
%{
uso_cor_l1b_echo=int32(USO_correction*cst.c/2*1e4);
int_path_cor_ku_l1b_echo=int32(instrument_range_correction_tx_rx*1e4);
if strcmp(cnf.mode,'SIN') && strcmp(mission,'CR2')&& strcmp(cnf.processing_mode,'SIN')
    int_path_2_cor_ku_l1b_echo=int32(instrument_range_correction*1e4);
end
%}
range_rate_l1b_echo=int32(T0_sar_surf_nadir*cst.c/2*1e3);
win_delay_l1b_echo = double(L1BS.win_delay_surf*1e3);

%---------- J. AGC and Sigma0 scalings --------------------------------
scale_factor_ku_l1b_echo=int32(L1B.wfm_scaling_factor*1e2);% sigma zero scaling factor

%---------- K. Stack characterization--------------------------------------
% switch netcdf_type
%     case 'netcdf4'
%         nb_stack_l1b_echo=uint16(L1BS.N_beams_stack');
%     case 'netcdf3'
%         nb_stack_l1b_echo=int32(L1BS.N_beams_stack');
% end
nb_stack_l1b_echo=int32(L1B.N_beams_contributing);
nb_stack_start_stop_l1b_echo=int32(L1B.N_beams_start_stop');
%     aux=squeeze(L1BS.beams_rng_cmpr(1:L1BS.N_beams_stack,:));
%     finite_indx=isfinite(aux);
%     max_stack_l1b_echo = uint32(max(aux(finite_indx))*1e2);
%     clear aux;
if L1B.start_beam~=0
    look_angle_start_l1b_echo = int16(L1BS.look_ang_surf(L1B.start_beam)'*1e6);
    look_angle_stop_l1b_echo = int16(L1BS.look_ang_surf(L1B.stop_beam)*1e6);
    doppler_angle_start_l1b_echo = int16((L1BS.doppler_ang_surf(L1B.start_beam))'*1e6);
    doppler_angle_stop_l1b_echo  = int16((L1BS.doppler_ang_surf(L1B.stop_beam))'*1e6);
    pointing_angle_start_l1b_echo = int16(L1BS.pointing_ang_surf(L1B.start_beam)*1e6)';
    pointing_angle_stop_l1b_echo  = int16(L1BS.pointing_ang_surf(L1B.stop_beam)'*1e6);
else
    look_angle_start_l1b_echo = int16(0);
    look_angle_stop_l1b_echo = int16(0);
    doppler_angle_start_l1b_echo = int16(0);
    doppler_angle_stop_l1b_echo  = int16(0);
    pointing_angle_start_l1b_echo = int16(0);
    pointing_angle_stop_l1b_echo  = int16(0);
end

% switch netcdf_type
%     case 'netcdf4'
%         stdev_stack_l1b_echo = uint32(L1B.stack_std' * 1e6);
%     case 'netcdf3'
%         stdev_stack_l1b_echo = int64(L1B.stack_std' * 1e6);
% end
if cnf.compute_L1_stadistics
    stdev_stack_l1b_echo = int64(L1B.stack_std' * 1e6);
    skew_stack_l1b_echo = int32(L1B.stack_skewness' * 1e6);
    kurt_stack_l1b_echo = int32(L1B.stack_kurtosis' * 1e6);
    gaussian_fitting_centre_look_l1b_echo = int16(L1B.stack_look_ang_centre' * 1e6);
    gaussian_fitting_centre_pointing_l1b_echo = int16(L1B.stack_pointing_ang_centre' * 1e6);    
end

%------------- L. Altimeter engineering variables -------------------------
altimeter_clock_l1b_echo = int32(1./T0_sar_surf_nadir - chd.bw)* 1e9;
pri_l1b_echo = int64(pri_sar_surf  * 1e12);
% switch netcdf_type
%     case 'netcdf4'
%         pri_lrm_l1b_echo = uint32(pri_sar_isp_surf .* T0_sar_surf_nadir * 1e19);
%     case 'netcdf3'
%         pri_lrm_l1b_echo = int64(pri_sar_isp_surf .* T0_sar_surf_nadir * 1e19);
% end


%------------- M. Waveform related variables -----------------------------
waveform_scale_factor_l1b_echo = single(max(L1B.wfm_cor_i2q2)/ (2^16-2));
i2q2_meas_ku_l1b_echo = int32(round(L1B.wfm_cor_i2q2./ waveform_scale_factor_l1b_echo));

if strcmp(cnf.processing_mode,'SIN')
    waveform_scale_factor_l1b_echo_2 = single(max(L1B.wfm_cor_i2q2_2)/ (2^16-2));
    i2q2_meas_ku_l1b_echo_2 = int32(round(L1B.wfm_cor_i2q2_2./ waveform_scale_factor_l1b_echo_2));
    phase_diff_meas_ku_l1b_echo = double(L1B.phase_difference);
    coherence_meas_ku_l1b_echo = double(L1B.coherence);
end

if cnf.include_wfms_aligned
    i2q2_meas_ku_wdcorr_l1b_echo = int32(round(L1B.wfm_cor_i2q2_wdcorr./ waveform_scale_factor_l1b_echo));
end
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

%% PACKING L1B
%----------A. Time variables ----------------------------------------------
var_id = netcdf.inqVarID(files.ncid,'time_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,time_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'UTC_day_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,UTC_day_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'UTC_sec_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,UTC_sec_l1b_echo);

%----------B. Orbit and attitude variables ------------------------------
var_id = netcdf.inqVarID(files.ncid,'lat_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,lat_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'lon_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,lon_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'alt_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,alt_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'orb_alt_rate_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,orb_alt_rate_l1b_echo);
%{
var_id = netcdf.inqVarID(files.ncid,'satellite_mispointing_l1b_sar_echo_ku');
netcdf.putVar(files.ncid,var_id,[0 i_surf],[3 1],satellite_mispointing_l1b_sar_echo_ku);
%}

%----------C. Flag time variables --------------------------------------

%----------D. Position/Velocity variables ------------------------------
var_id = netcdf.inqVarID(files.ncid,'x_pos_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,x_pos_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'y_pos_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,y_pos_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'z_pos_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,z_pos_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'x_vel_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,x_vel_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'y_vel_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,y_vel_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'z_vel_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,z_vel_l1b_echo);

%----------E. Navigation Bulletin ------------------------------

% var_id = netcdf.inqVarID(files.ncid,'seq_count_l1b_echo');
% netcdf.putVar(files.ncid,var_id,i_surf_stacked,seq_count_l1b_echo);

%----------F. Operating instrument & tracking --------------------------
% var_id = netcdf.inqVarID(files.ncid,'oper_instr_l1b_echo');
% netcdf.putVar(files.ncid,var_id,i_surf_stacked,oper_instr_l1b_echo);

% var_id = netcdf.inqVarID(files.ncid,'SAR_cnf.mode_l1b_echo');
% netcdf.putVar(files.ncid,var_id,i_surf_stacked,SAR_cnf.mode_l1b_echo);

%---------- G. H0, COR2 and AGC ----------------------------------------
var_id = netcdf.inqVarID(files.ncid,'h0_applied_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,h0_applied_l1b_echo);
%ncwrite(files.filename_netCDF,'h0_applied_l1b_echo',h0_applied_l1b_echo,i_surf_stacked);

var_id = netcdf.inqVarID(files.ncid,'cor2_applied_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,cor2_applied_l1b_echo);

% var_id = netcdf.inqVarID(files.ncid,'agccode_ku_l1b_echo');
% netcdf.putVar(files.ncid,var_id,i_surf_stacked,agccode_ku_l1b_echo);

%---------- H. Surface type -----------------------------------------------
var_id = netcdf.inqVarID(files.ncid,'surf_type_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,surf_type_l1b_echo);

%---------- I. Altimeter range and Corrections ----------------------------
var_id = netcdf.inqVarID(files.ncid,'range_ku_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,range_ku_l1b_echo);
if cnf.include_wfms_aligned
    %EM:14.04.2016
    var_id = netcdf.inqVarID(files.ncid,'range_ku_surf_aligned_l1b_echo');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,range_ku_surf_aligned_l1b_echo);
end

% var_id = netcdf.inqVarID(files.ncid,'uso_cor_l1b_echo');
% netcdf.putVar(files.ncid,var_id,i_surf_stacked,uso_cor_l1b_echo);

% var_id = netcdf.inqVarID(files.ncid,'int_path_cor_ku_l1b_echo');
% netcdf.putVar(files.ncid,var_id,i_surf_stacked,int_path_cor_ku_l1b_echo);
%{
if strcmp(cnf.mode,'SIN') && strcmp(mission,'CR2')&& strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(files.ncid,'int_path_2_cor_ku_l1b_echo');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,int_path_2_cor_ku_l1b_echo);
    
end
%}
var_id = netcdf.inqVarID(files.ncid,'range_rate_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,range_rate_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'win_delay_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,win_delay_l1b_echo);

%---------- J. AGC and Sigma0 scalings --------------------------------
var_id = netcdf.inqVarID(files.ncid,'scale_factor_ku_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,scale_factor_ku_l1b_echo);

%---------- K. Stack characterization--------------------------------------
var_id = netcdf.inqVarID(files.ncid,'nb_stack_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,nb_stack_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'nb_stack_start_stop_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,nb_stack_start_stop_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'look_angle_start_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,look_angle_start_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'look_angle_stop_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,look_angle_stop_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'doppler_angle_start_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,doppler_angle_start_l1b_echo);
 
var_id = netcdf.inqVarID(files.ncid,'doppler_angle_stop_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,doppler_angle_stop_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'pointing_angle_start_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,pointing_angle_start_l1b_echo);
 
var_id = netcdf.inqVarID(files.ncid,'pointing_angle_stop_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,pointing_angle_stop_l1b_echo);

if cnf.compute_L1_stadistics
    var_id = netcdf.inqVarID(files.ncid,'skew_stack_l1b_echo');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,skew_stack_l1b_echo);
    
    var_id = netcdf.inqVarID(files.ncid,'kurt_stack_l1b_echo');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,kurt_stack_l1b_echo);
    
    var_id = netcdf.inqVarID(files.ncid,'stdev_stack_l1b_echo');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,stdev_stack_l1b_echo);
%   ncwrite(files.filename_netCDF,'stdev_stack_l1b_echo',stdev_stack_l1b_echo,i_surf_stacked);
    
    var_id = netcdf.inqVarID(files.ncid,'gaussian_fitting_centre_look_l1b_echo');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,gaussian_fitting_centre_look_l1b_echo);
    
    var_id = netcdf.inqVarID(files.ncid,'gaussian_fitting_centre_pointing_l1b_echo');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,gaussian_fitting_centre_pointing_l1b_echo);
    
    var_id = netcdf.inqVarID(files.ncid,'beam_form_l1b_echo');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,int32(100.0*1e2));
    %ncwrite(files.filename_netCDF,'beam_form_l1b_echo',uint16(100.0*1e2),i_surf_stacked);
    
end

%------------- L. Altimeter engineering variables -------------------------
var_id = netcdf.inqVarID(files.ncid,'altimeter_clock_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,altimeter_clock_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'pri_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,pri_l1b_echo);
%ncwrite(files.filename_netCDF,'pri_lrm_l1b_echo',pri_lrm_l1b_echo,i_surf_stacked);

%------------- M. Waveform related variables -----------------------------
var_id = netcdf.inqVarID(files.ncid,'i2q2_meas_ku_l1b_echo');
netcdf.putVar(files.ncid,var_id,[0 i_surf_stacked],[chd.N_samples_sar*cnf.zp_fact_range 1], i2q2_meas_ku_l1b_echo.');
%ncwrite(files.filename_netCDF,'i2q2_meas_ku_l1b_echo',i2q2_meas_ku_l1b_echo.',[1 i_surf_stacked]);

if strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(files.ncid,'i2q2_meas_ku_l1b_echo_2');
    netcdf.putVar(files.ncid,var_id,[0 i_surf_stacked],[chd.N_samples_sar*cnf.zp_fact_range 1], i2q2_meas_ku_l1b_echo_2.');   
    var_id = netcdf.inqVarID(files.ncid,'waveform_scale_factor_l1b_echo_2');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,waveform_scale_factor_l1b_echo_2);
    var_id = netcdf.inqVarID(files.ncid,'phase_diff_meas_ku_l1b_echo');
    netcdf.putVar(files.ncid,var_id,[0 i_surf_stacked],[chd.N_samples_sar*cnf.zp_fact_range 1], phase_diff_meas_ku_l1b_echo.');
    var_id = netcdf.inqVarID(files.ncid,'coherence_meas_ku_l1b_echo');
    netcdf.putVar(files.ncid,var_id,[0 i_surf_stacked],[chd.N_samples_sar*cnf.zp_fact_range 1], coherence_meas_ku_l1b_echo.');
end


if cnf.include_wfms_aligned
    % %EM:14.04.2016
    var_id = netcdf.inqVarID(files.ncid,'i2q2_meas_ku_wdcorr_l1b_echo');
    netcdf.putVar(files.ncid,var_id,[0 i_surf_stacked],[chd.N_samples_sar*cnf.zp_fact_range 1], i2q2_meas_ku_wdcorr_l1b_echo.');
    %ncwrite(files.filename_netCDF,'i2q2_meas_ku_wdcorr_l1b_echo',i2q2_meas_ku_wdcorr_l1b_echo.',[1 i_surf_stacked]);
end

var_id = netcdf.inqVarID(files.ncid,'stack_mask_range_bin_l1b_echo');
netcdf.putVar(files.ncid,var_id,[0 i_surf_stacked],[chd.N_max_beams_stack 1], int16(ceil(L1B.stack_mask_vector/cnf.zp_fact_range)).');


var_id = netcdf.inqVarID(files.ncid,'waveform_scale_factor_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,waveform_scale_factor_l1b_echo);

var_id = netcdf.inqVarID(files.ncid,'scale_factor_ku_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,scale_factor_ku_l1b_echo);

% EM: 22.09.2016
%{
if strcmp(cnf.mode,'SIN') && strcmp(mission,'CR2') && strcmp(cnf.processing_mode,'SIN')
    var_id = netcdf.inqVarID(files.ncid,'i2q2_meas_ku_l1b_echo_2');
    netcdf.putVar(files.ncid,var_id,[0 i_surf],[chd.N_samples_sar*cnf.zp_fact_range 1], i2q2_meas_ku_l1b_echo_2.');
    
    var_id = netcdf.inqVarID(files.ncid,'coherence_meas_ku_l1b_echo');
    netcdf.putVar(files.ncid,var_id,[0 i_surf],[chd.N_samples_sar*cnf.zp_fact_range 1], coherence_meas_ku_l1b_echo.');
    
    var_id = netcdf.inqVarID(files.ncid,'phase_difference_meas_ku_l1b_echo');
    netcdf.putVar(files.ncid,var_id,[0 i_surf],[chd.N_samples_sar*cnf.zp_fact_range 1], phase_difference_meas_ku_l1b_echo.');               
    
    if cnf.include_wfms_aligned
        
        var_id = netcdf.inqVarID(files.ncid,'i2q2_meas_ku_wdcorr_l1b_echo_2');
        netcdf.putVar(files.ncid,var_id,[0 i_surf],[chd.N_samples_sar*cnf.zp_fact_range 1], i2q2_meas_ku_wdcorr_l1b_echo_2.');
        
        var_id = netcdf.inqVarID(files.ncid,'coherence_meas_ku_wdcorr_l1b_echo');
        netcdf.putVar(files.ncid,var_id,[0 i_surf],[chd.N_samples_sar*cnf.zp_fact_range 1], coherence_meas_ku_wdcorr_l1b_echo.');        
        
        var_id = netcdf.inqVarID(files.ncid,'phase_difference_meas_ku_wdcorr_l1b_echo');
        netcdf.putVar(files.ncid,var_id,[0 i_surf],[chd.N_samples_sar*cnf.zp_fact_range 1], phase_difference_meas_ku_wdcorr_l1b_echo.');                

    end
           
end
%}


%------------- N. Geophysical Corrections variables -----------------------
%{
var_id = netcdf.inqVarID(files.ncid,'dry_tropo_correction_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,dry_tropo_correction);

var_id = netcdf.inqVarID(files.ncid,'wet_tropo_correction_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,wet_tropo_correction);

var_id = netcdf.inqVarID(files.ncid,'inverse_baro_correction_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,inverse_baro_correction);

var_id = netcdf.inqVarID(files.ncid,'Dynamic_atmospheric_correction_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,Dynamic_atmospheric_correction);

var_id = netcdf.inqVarID(files.ncid,'GIM_iono_correction_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,GIM_iono_correction);

var_id = netcdf.inqVarID(files.ncid,'cnf.model_iono_correction_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,cnf.model_iono_correction);

var_id = netcdf.inqVarID(files.ncid,'ocean_equilibrium_tide_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,ocean_equilibrium_tide);

var_id = netcdf.inqVarID(files.ncid,'long_period_tide_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,long_period_tide_height);

var_id = netcdf.inqVarID(files.ncid,'ocean_loading_tide_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,ocean_loading_tide);

var_id = netcdf.inqVarID(files.ncid,'solid_earth_tide_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,solid_earth_tide);

var_id = netcdf.inqVarID(files.ncid,'geocentric_polar_tide_l1b_echo');
netcdf.putVar(files.ncid,var_id,i_surf_stacked,geocentric_polar_tide);
%}
% %------------- O. Processing parameters used ------------------------------

% %------------- P. ACDC variables ----------------------------------------
%added EM 04.10.2016
%{
if ACDC_application_cnf
    %-------------- Waveform & Scaling ------------------------------------
    %ACDC waveform
    var_id = netcdf.inqVarID(files.ncid,'i2q2_meas_ku_ACDC_echo');
    netcdf.putVar(files.ncid,var_id,[0 i_surf],[chd.N_samples_sar*cnf.zp_fact_range 1], i2q2_meas_ku_ACDC_echo.');
    %ncwrite(files.filename_netCDF,'i2q2_meas_ku_ACDC_echo',i2q2_meas_ku_ACDC_echo.',[1 i_surf_stacked]);    
    var_id = netcdf.inqVarID(files.ncid,'waveform_scale_factor_ACDC_echo');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,waveform_scale_factor_ACDC_echo);
    
    %ACDC fitted waveform    
    var_id = netcdf.inqVarID(files.ncid,'i2q2_fit_ku_ACDC_echo');
    netcdf.putVar(files.ncid,var_id,[0 i_surf],[chd.N_samples_sar*cnf.zp_fact_range 1], i2q2_fit_ku_ACDC_echo.');
    %ncwrite(files.filename_netCDF,'i2q2_meas_ku_ACDC_echo',i2q2_meas_ku_ACDC_echo.',[1 i_surf_stacked]);    
    var_id = netcdf.inqVarID(files.ncid,'waveform_scale_factor_fit_ACDC_echo');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,waveform_scale_factor_fit_ACDC_echo);
    
    
    %L1B fitted waveform        
    var_id = netcdf.inqVarID(files.ncid,'i2q2_fit_ku_echo');
    netcdf.putVar(files.ncid,var_id,[0 i_surf],[chd.N_samples_sar*cnf.zp_fact_range 1], i2q2_fit_ku_echo.');
    %ncwrite(files.filename_netCDF,'i2q2_meas_ku_ACDC_echo',i2q2_meas_ku_ACDC_echo.',[1 i_surf_stacked]);    
    var_id = netcdf.inqVarID(files.ncid,'waveform_scale_factor_fit_echo');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,waveform_scale_factor_fit_echo);
    
    
    %-------------- Range index -------------------------------------------
    var_id = netcdf.inqVarID(files.ncid,'range_index_ACDC_echo');
    netcdf.putVar(files.ncid,var_id,[0 i_surf],[chd.N_samples_sar*cnf.zp_fact_range 1], range_index_ACDC_echo.');
    %ncwrite(files.filename_netCDF,'range_index_ACDC_echo',range_index_ACDC_echo.',[1 i_surf_stacked]);
    
    %-------------- Geophysical retrievals --------------------------------
    %******************* retracked range **********************************
    var_id = netcdf.inqVarID(files.ncid,'retracked_range_ACDC_20_ku');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,retracked_range_ACDC_20_ku);
    %******************* epoch ********************************************
    var_id = netcdf.inqVarID(files.ncid,'epoch_ACDC_20_ku');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,epoch_ACDC_20_ku);
    %******************** SWH  ********************************************
    var_id = netcdf.inqVarID(files.ncid,'swh_ACDC_20_ku');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,swh_ACDC_20_ku);
    var_id = netcdf.inqVarID(files.ncid,'swh_ACDC_PRE_20_ku');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,swh_ACDC_PRE_20_ku);
    %******************** SSH  ********************************************
    var_id = netcdf.inqVarID(files.ncid,'ssh_ACDC_20_ku');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,ssh_ACDC_20_ku);
    var_id = netcdf.inqVarID(files.ncid,'ssh_ACDC_PRE_20_ku');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,ssh_ACDC_PRE_20_ku);    
    %******************** SIGMA0 ******************************************
    var_id = netcdf.inqVarID(files.ncid,'sig0_ACDC_20_ku');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,sig0_ACDC_20_ku);
    
    var_id = netcdf.inqVarID(files.ncid,'sig0_ACDC_PRE_20_ku');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,sig0_ACDC_PRE_20_ku);
    %********************* Pu: fitted peak power **************************
    var_id = netcdf.inqVarID(files.ncid,'Pu_analytical_ACDC_20_ku');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,Pu_analytical_ACDC_20_ku);
    %********************** Pearson Correlation Coefficient ***************
    var_id = netcdf.inqVarID(files.ncid,'Pearson_corr_ACDC_20_ku');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,Pearson_corr_ACDC_20_ku);
    
    var_id = netcdf.inqVarID(files.ncid,'Pearson_corr_ACDC_PRE_20_ku');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,Pearson_corr_ACDC_PRE_20_ku);
    %************************ Fitting flag ********************************    
    var_id = netcdf.inqVarID(files.ncid,'Flag_fitting_ACDC_20_ku');
    netcdf.putVar(files.ncid,var_id,i_surf_stacked,Flag_fitting_ACDC_20_ku);
    
end
%}

% %----------  Global Attributes definition -----------------------------------
% %---- attributes inherited from Sentinel-3 product description-------------
% ncwriteatt(files.filename_netCDF,'/','creation_time',date_creation);
% ncwriteatt(files.filename_netCDF,'/','Conventions',netcdf_v4_format);
% ncwriteatt(files.filename_netCDF,'/','mission_name',mission);
% ncwriteatt(files.filename_netCDF,'/','altimeter_sensor_name',altimeter_sensor_name);
% ncwriteatt(files.filename_netCDF,'/','gnss_sensor_name',gnss_sensor_name);
% ncwriteatt(files.filename_netCDF,'/','doris_sensor_name',doris_sensor_name);
% ncwriteatt(files.filename_netCDF,'/','doris_sensor_name',acq_station_name);
% ncwriteatt(files.filename_netCDF,'/','doris_sensor_name',acq_station_name);
% ncwriteatt(files.filename_netCDF,'/','first_meas_time',first_meas_time);
% ncwriteatt(files.filename_netCDF,'/','last_meas_time',last_meas_time);
% ncwriteatt(files.filename_netCDF,'/','xref_altimeter_level0',xref_altimeter_level0);
% ncwriteatt(files.filename_netCDF,'/','xref_altimeter_orbit',xref_altimeter_orbit);
% ncwriteatt(files.filename_netCDF,'/','xref_doris_USO',xref_doris_USO);
% ncwriteatt(files.filename_netCDF,'/','xref_altimeter_ltm_sar_cal1',xref_altimeter_ltm_sar_cal1);
% ncwriteatt(files.filename_netCDF,'/','xref_altimeter_ltm_ku_cal2',xref_altimeter_ltm_ku_cal2);
% ncwriteatt(files.filename_netCDF,'/','xref_altimeter_ltm_c_cal2',xref_altimeter_ltm_c_cal2);
% ncwriteatt(files.filename_netCDF,'/','xref_altimeter_characterisation',xref_altimeter_characterisation);
% ncwriteatt(files.filename_netCDF,'/','semi_major_ellipsoid_axis',semi_major_ellipsoid_axis);
% ncwriteatt(files.filename_netCDF,'/','ellipsoid_flattening',ellipsoid_flattening);
% %--------------- add the attributes related to intermediate product--------
% ncwriteatt(files.filename_netCDF,'/','orbit_phase_code',orbit_phase_code);
% ncwriteatt(files.filename_netCDF,'/','orbit_cycle_num',orbit_cycle_num);
% ncwriteatt(files.filename_netCDF,'/','orbit_REL_Orbit',orbit_REL_Orbit);
% ncwriteatt(files.filename_netCDF,'/','orbit_ABS_Orbit_Start',orbit_ABS_Orbit_Start);
% ncwriteatt(files.filename_netCDF,'/','orbit_Rel_Time_ASC_Node_Start',orbit_Rel_Time_ASC_Node_Start);
% ncwriteatt(files.filename_netCDF,'/','orbit_ABS_Orbit_Stop',orbit_ABS_Orbit_Stop);
% ncwriteatt(files.filename_netCDF,'/','orbit_Rel_Time_ASC_Node_Stop',orbit_Rel_Time_ASC_Node_Stop);
% ncwriteatt(files.filename_netCDF,'/','orbit_Equator_Cross_Time',orbit_Equator_Cross_Time);
% ncwriteatt(files.filename_netCDF,'/','orbit_Equator_Cross_Long',orbit_Equator_Cross_Long);
% ncwriteatt(files.filename_netCDF,'/','orbit_Ascending_Flag',orbit_Ascending_Flag);

% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);
% netcdf.close(files.ncid);