
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
% 2.1 changed N_samples_sar_chd by N_samples
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function prepare_NetCDF_L1B_pLRM(files,L1A_buffer,L1B,i_RC)
i_rc = i_RC-1;
global c_cst bw_ku_chd  zp_fact_range_cnf sec_in_day_cst pi_cst
global mission mode N_samples N_max_beams_stack_chd
global compute_L1_stadistics_cnf
global ref_burst

%global netcdf_type
% ncid = netcdf.open(files.filename_netCDF,'WRITE');

% ku_rec_dimension = netcdf.defDim(ncid,'time_l1b_echo_plrm',N_surfs_loc_estimated);
% nl_dimension = netcdf.defDim(ncid,'Nl',N_max_beams_stack_chd);
% ns_dimension = netcdf.defDim(ncid,'Ns',N_samples*zp_fact_range_cnf);
% space_3D_dimension = netcdf.defDim(ncid,'space_3D',3);

% t6 = tic;

% date_creation = datestr(now, '_yyyymmdd_');
% switch mission
%     case 'CR2'        
%         switch mode
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
% ku_rec_dimension = 'time_l1b_echo_plrm');
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

source_seq_count_sar_isp_surf = L1A_buffer.source_seq_count_sar_isp;
ins_id=L1A_buffer.ins_id;
ins_loop_stat=L1A_buffer.ins_loop_stat;
h0_comp_sar_isp=L1A_buffer.h0_comp_sar_isp;
cor2_comp_sar_isp=L1A_buffer.cor2_comp_sar_isp;
ATT1_science=L1A_buffer.ATT1_science;    
ATT2_science=L1A_buffer.ATT2_science;     
surf_type_flag = L1A_buffer.surface_type_flag_bursts;
USO_correction = L1A_buffer.USO_correction;
instrument_range_correction_tx_rx= L1A_buffer.instrument_range_correction_tx_rx;

if strcmp(mode,'SIN') && strcmp(mission,'CR2')
 instrument_range_correction=L1A_buffer.instrument_range_correction_rx;
end
T0_sar_surf_nadir = L1A_buffer.T0_sar;
%      ATT1_science_corr=L1A_buffer.ATT1_science_corr;  
%      ATT2_science_corr=L1A_buffer.ATT2_science_corr;  
%      ATT1_delta_corr=L1A_buffer.ATT1_delta_corr;
%      ATT2_delta_corr=L1A_buffer.ATT2_delta_corr;  

pri_sar_isp_surf                = L1A_buffer.pri_sar_isp;

%geophysical corrections
dry_tropo_correction            = L1A_buffer.dry_tropo_correction_bursts;
wet_tropo_correction            = L1A_buffer.wet_tropo_correction_bursts;
inverse_baro_correction         = L1A_buffer.inverse_baro_correction_bursts;
Dynamic_atmospheric_correction  = L1A_buffer.Dynamic_atmospheric_correction_bursts;
GIM_iono_correction             = L1A_buffer.GIM_iono_correction_bursts;
model_iono_correction           = L1A_buffer.model_iono_correction_bursts;
ocean_equilibrium_tide          = L1A_buffer.ocean_equilibrium_tide_bursts;
long_period_tide_height         = L1A_buffer.long_period_tide_height_bursts;
ocean_loading_tide              = L1A_buffer.ocean_loading_tide_bursts;
solid_earth_tide                = L1A_buffer.solid_earth_tide_bursts;
geocentric_polar_tide           = L1A_buffer.geocentric_polar_tide_bursts;

%----------A. Time variables ----------------------------------------------
% leap seconds in 2010 34s, 
% +1 the 1st of July 2012 TAI_2012 = (12*365+4+181)*3600*24 + 35
% +1 the 1st of July 2015 TAI_2015 = (15*365+4+181)*3600*24 + 36

TAI_2012 = (12*365+4+181)*3600*24 + 35;
TAI_2015 = (15*365+4+181)*3600*24 + 36;
time_l1b_echo_plrm = double(L1B.time_l1b_plrm);

% add leap seconds to the TAI time. Only valid for 
if(time_l1b_echo_plrm < TAI_2012)
    time_l1b_echo_plrm = time_l1b_echo_plrm - 34;
elseif(time_l1b_echo_plrm > TAI_2015)
    time_l1b_echo_plrm = time_l1b_echo_plrm - 36;    
else
    time_l1b_echo_plrm = time_l1b_echo_plrm - 35;
end
  
UTC_day_l1b_echo_plrm = int16(floor(L1B.time_l1b_plrm./ sec_in_day_cst));
UTC_sec_l1b_echo_plrm = double((time_l1b_echo_plrm - double(UTC_day_l1b_echo_plrm) * sec_in_day_cst));
%isp_coarse_time_l1b_echo_plrm=uint32();
%isp_fine_time_l1b_echo_plrm=int32();
%sral_fine_time_l1b_echo_plrm=uint32();


%tm_source_sequence_counter_ku = uint16(source_seq_count_sar_isp_surf);
%l1b_record_counter_ku = uint16(0:(length(win_delay_sar)-1));

%----------B. Orbit and attitude variables ------------------------------
lat_l1b_echo_plrm = int32(L1B.latitude * 1e6);
lon_l1b_echo_plrm = int32(L1B.longitude * 1e6); 
alt_l1b_echo_plrm = int32((L1B.altitude-700000) * 1e4); 
orb_alt_rate_l1b_echo_plrm = int16(L1B.altitude_rate * 1e2);

satellite_mispointing_l1b_echo_plrm = int32([L1B.pitch * 180/pi_cst * 1e7; L1B.roll * 180/pi_cst * 1e7; L1B.yaw * 180/pi_cst * 1e7]);

%----------C. Flag time variables --------------------------------------

%----------D. Position/Velocity variables ------------------------------
x_pos_l1b_echo_plrm=double(L1B.x');
y_pos_l1b_echo_plrm=double(L1B.y');
z_pos_l1b_echo_plrm=double(L1B.z');

x_vel_l1b_echo_plrm=double(L1B.x_vel');
y_vel_l1b_echo_plrm=double(L1B.y_vel');
z_vel_l1b_echo_plrm=double(L1B.z_vel');

%----------E. Navigation Bulletin --------------------------------------
seq_count_l1b_echo_plrm=uint16(source_seq_count_sar_isp_surf);

%----------F. Operating instrument & tracking --------------------------
oper_instr_l1b_echo_plrm=int8(ins_id);%

SAR_mode_l1b_echo_plrm=int8(ins_loop_stat);


%---------- G. H0, COR2 and AGC ----------------------------------------
% switch netcdf_type
%     case 'netcdf4'
%         h0_applied_l1b_echo_plrm = uint32(h0_comp_sar_isp/(3.125/64*1e-9));
%     case 'netcdf3'
%         h0_applied_l1b_echo_plrm = int64(h0_comp_sar_isp/(3.125/64*1e-9));
% end
h0_applied_l1b_echo_plrm = int64(h0_comp_sar_isp/(3.125/64*1e-9));        
cor2_applied_l1b_echo_plrm=int16(cor2_comp_sar_isp/(3.125/1024*1e-9));
agccode_ku_l1b_echo_plrm=int8(-1.*(ATT1_science+ATT2_science));


%-----------H. Surface Type flag----------------------------------------
surf_type_l1b_echo_plrm =int8(surf_type_flag);


%---------- I. Altimeter range and Corrections -------------------------
range_ku_l1b_echo_plrm=int32((L1B.win_delay_sar*c_cst/2-700000)*1e4);

uso_cor_l1b_echo_plrm=int32(USO_correction*c_cst/2*1e4);
int_path_cor_ku_l1b_echo_plrm=int32(instrument_range_correction_tx_rx*1e4);
if strcmp(mode,'SIN') && strcmp(mission,'CR2')
    int_path_2_cor_ku_l1b_echo_plrm=int32(instrument_range_correction*1e4);
end
range_rate_l1b_echo_plrm=int32(T0_sar_surf_nadir*c_cst/2*1e3);

%---------- J. AGC and Sigma0 scalings --------------------------------
scale_factor_ku_l1b_echo_plrm=int32(L1B.sigma0_scaling_factor_RC*1e2);% sigma zero scaling factor


%------------- L. Altimeter engineering variables -------------------------
altimeter_clock_l1b_echo_plrm = int32(1./T0_sar_surf_nadir - bw_ku_chd)* 1e9;
pri_lrm_l1b_echo_plrm = int64(pri_sar_isp_surf .* T0_sar_surf_nadir * 1e19);
% switch netcdf_type
%     case 'netcdf4'
%         pri_lrm_l1b_echo_plrm = uint32(pri_sar_isp_surf .* T0_sar_surf_nadir * 1e19);
%     case 'netcdf3'
%         pri_lrm_l1b_echo_plrm = int64(pri_sar_isp_surf .* T0_sar_surf_nadir * 1e19);
% end


%------------- M. Waveform related variables -----------------------------
waveform_scale_factor_l1b_echo_plrm = single(max(L1B.wfm_L1B_RC)/ (2^16-2));
i2q2_meas_ku_l1b_echo_plrm = int32(round(L1B.wfm_L1B_RC./ waveform_scale_factor_l1b_echo_plrm));

% switch netcdf_type
%     case 'netcdf4'
%         i2q2_meas_ku_l1b_echo_plrm = uint16(round(L1B.wfm_cor_i2q2_sar_ku.*10^0.3 ./ waveform_scale_factor_l1b_echo_plrm));
%         i2q2_meas_ku_wdcorr_l1b_echo_plrm = uint16(round(L1B.wfm_cor_i2q2_sar_ku_wdcorr.*10^0.3 ./ waveform_scale_factor_l1b_echo_plrm));
%     case 'netcdf3'
%         i2q2_meas_ku_l1b_echo_plrm = int32(round(L1B.wfm_cor_i2q2_sar_ku.*10^0.3 ./ waveform_scale_factor_l1b_echo_plrm));
%         i2q2_meas_ku_wdcorr_l1b_echo_plrm = int32(round(L1B.wfm_cor_i2q2_sar_ku_wdcorr.*10^0.3 ./ waveform_scale_factor_l1b_echo_plrm));
% end



%------------- N. Geophysical Corrections variables ---------------------
dry_tropo_correction=int32(dry_tropo_correction.*1e3);
wet_tropo_correction=int32(wet_tropo_correction.*1e3);
inverse_baro_correction=int32(inverse_baro_correction.*1e3);
Dynamic_atmospheric_correction=int32(Dynamic_atmospheric_correction.*1e3);
GIM_iono_correction=int32(GIM_iono_correction.*1e3);
model_iono_correction=int32(model_iono_correction.*1e3);
ocean_equilibrium_tide=int32(ocean_equilibrium_tide.*1e3);
long_period_tide_height=int32(long_period_tide_height.*1e3);
ocean_loading_tide=int32(ocean_loading_tide.*1e3);
solid_earth_tide=int32(solid_earth_tide.*1e3);
geocentric_polar_tide=int32(geocentric_polar_tide.*1e3);




%% PACKING L1B
%----------A. Time variables ----------------------------------------------
var_id = netcdf.inqVarID(files.ncid,'time_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,time_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'UTC_day_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,UTC_day_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'UTC_sec_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,UTC_sec_l1b_echo_plrm);

%----------B. Orbit and attitude variables ------------------------------
var_id = netcdf.inqVarID(files.ncid,'lat_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,lat_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'lon_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,lon_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'alt_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,alt_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'orb_alt_rate_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,orb_alt_rate_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'satellite_mispointing_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,[0 i_rc],[3 1],satellite_mispointing_l1b_echo_plrm);


%----------C. Flag time variables --------------------------------------

%----------D. Position/Velocity variables ------------------------------
var_id = netcdf.inqVarID(files.ncid,'x_pos_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,x_pos_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'y_pos_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,y_pos_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'z_pos_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,z_pos_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'x_vel_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,x_vel_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'y_vel_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,y_vel_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'z_vel_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,z_vel_l1b_echo_plrm);

%----------E. Navigation Bulletin ------------------------------

% var_id = netcdf.inqVarID(files.ncid,'seq_count_l1b_echo_plrm');
% netcdf.putVar(files.ncid,var_id,i_surf,seq_count_l1b_echo_plrm);

%----------F. Operating instrument & tracking --------------------------
var_id = netcdf.inqVarID(files.ncid,'oper_instr_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,oper_instr_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'SAR_mode_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,SAR_mode_l1b_echo_plrm);

%---------- G. H0, COR2 and AGC ----------------------------------------
var_id = netcdf.inqVarID(files.ncid,'h0_applied_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,h0_applied_l1b_echo_plrm);
%ncwrite(files.filename_netCDF,'h0_applied_l1b_echo_plrm',h0_applied_l1b_echo_plrm,i_surf_stacked);

var_id = netcdf.inqVarID(files.ncid,'cor2_applied_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,cor2_applied_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'agccode_ku_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,agccode_ku_l1b_echo_plrm);

%---------- H. Surface type -----------------------------------------------
var_id = netcdf.inqVarID(files.ncid,'surf_type_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,surf_type_l1b_echo_plrm);

%---------- I. Altimeter range and Corrections ----------------------------
var_id = netcdf.inqVarID(files.ncid,'range_ku_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,range_ku_l1b_echo_plrm);


var_id = netcdf.inqVarID(files.ncid,'uso_cor_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,uso_cor_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'int_path_cor_ku_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,int_path_cor_ku_l1b_echo_plrm);

if strcmp(mode,'SIN') && strcmp(mission,'CR2')
    var_id = netcdf.inqVarID(files.ncid,'int_path_2_cor_ku_l1b_echo_plrm');
    netcdf.putVar(files.ncid,var_id,i_rc,int_path_2_cor_ku_l1b_echo_plrm);
    
end

var_id = netcdf.inqVarID(files.ncid,'range_rate_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,range_rate_l1b_echo_plrm);

%---------- J. AGC and Sigma0 scalings --------------------------------
var_id = netcdf.inqVarID(files.ncid,'scale_factor_ku_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,scale_factor_ku_l1b_echo_plrm);

%------------- L. Altimeter engineering variables -------------------------
var_id = netcdf.inqVarID(files.ncid,'altimeter_clock_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,altimeter_clock_l1b_echo_plrm);

var_id = netcdf.inqVarID(files.ncid,'pri_lrm_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,pri_lrm_l1b_echo_plrm);
%ncwrite(files.filename_netCDF,'pri_lrm_l1b_echo_plrm',pri_lrm_l1b_echo_plrm,i_surf_stacked);

%------------- M. Waveform related variables -----------------------------
var_id = netcdf.inqVarID(files.ncid,'i2q2_meas_ku_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,[0 i_rc],[N_samples*zp_fact_range_cnf 1], i2q2_meas_ku_l1b_echo_plrm.');
%ncwrite(files.filename_netCDF,'i2q2_meas_ku_l1b_echo_plrm',i2q2_meas_ku_l1b_echo_plrm.',[1 i_surf_stacked]);


var_id = netcdf.inqVarID(files.ncid,'waveform_scale_factor_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,waveform_scale_factor_l1b_echo_plrm);


var_id = netcdf.inqVarID(files.ncid,'zero_padding_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,zp_fact_range_cnf);


%------------- N. Geophysical Corrections variables -----------------------

var_id = netcdf.inqVarID(files.ncid,'dry_tropo_correction_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,dry_tropo_correction);

var_id = netcdf.inqVarID(files.ncid,'wet_tropo_correction_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,wet_tropo_correction);

var_id = netcdf.inqVarID(files.ncid,'inverse_baro_correction_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,inverse_baro_correction);

var_id = netcdf.inqVarID(files.ncid,'Dynamic_atmospheric_correction_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,Dynamic_atmospheric_correction);

var_id = netcdf.inqVarID(files.ncid,'GIM_iono_correction_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,GIM_iono_correction);

var_id = netcdf.inqVarID(files.ncid,'model_iono_correction_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,model_iono_correction);

var_id = netcdf.inqVarID(files.ncid,'ocean_equilibrium_tide_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,ocean_equilibrium_tide);

var_id = netcdf.inqVarID(files.ncid,'long_period_tide_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,long_period_tide_height);

var_id = netcdf.inqVarID(files.ncid,'ocean_loading_tide_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,ocean_loading_tide);

var_id = netcdf.inqVarID(files.ncid,'solid_earth_tide_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,solid_earth_tide);

var_id = netcdf.inqVarID(files.ncid,'geocentric_polar_tide_l1b_echo_plrm');
netcdf.putVar(files.ncid,var_id,i_rc,geocentric_polar_tide);

% %------------- O. Processing parameters used ------------------------------
% 
% 
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