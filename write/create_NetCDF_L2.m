%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the CODING & PACKING 
% algorithm for Level-2 product using Sentinel-3 like format
% Ref: S3-TN-ESA-SR-0433 SRAL L1A-1BS IODD V1.4
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
% Last rev.: Monica Roca   / isardSAT 
% 
% Versions
% 1.0 
% 1.1 Updated time conversion for data between 2010 and 2016 (CR2)
% 2.0 Transformed to a function. Writting one record
% 2.1 zp changed to int16 as int8 did not work for zp > 128
% 2.2 case S3_ filename
% 2.3 case 'SIN' added 
% 2.4 added cnf flag processing_mode_cnf
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function [files] = create_NetCDF_L2(files, N_bursts, cnf, chd, cst, resize)

%global N_samples   semi_major_axis_cst flat_coeff_cst;
global bw_ku_chd 
%zp_fact_range_cnf  sec_in_day_cst pi_cst
%global mission mode N_max_beams_stack_chd
%global compute_L1_stadistics_cnf include_wfms_aligned
%global optional_ext_file_flag file_ext_string
%added by EM: 04.10.2016 
%global ACDC_application_cnf cnf_p_ACDC
%global netcdf_type
%global processing_mode_cnf
% t6 = tic;

date_creation = datestr(now, '_yyyymmddTHHMMSS_');

L1A_date_find = strfind(files.filename_L1A,'_201');
if (~resize)
    files.filename_L2 = strcat(files.outputPath,'SR_2_SRA____',...
        files.filename_L1A((L1A_date_find(1)-8):(L1A_date_find(1)-3)),'L2_',...
        files.filename_L1A((L1A_date_find(1)+1):(L1A_date_find(1)+30)),...
        '_isd','_long','.nc');
else
    files.filename_L2 = strcat(files.outputPath,'SR_2_SRA____',...
        files.filename_L1A((L1A_date_find(1)-8):(L1A_date_find(1)-3)),'L2_',...
        files.filename_L1A((L1A_date_find(1)+1):(L1A_date_find(1)+30)),...
        '_isd','.nc');
end

%{
files.filename_L1B = strcat(files.outputPath,'SR_1_SRA____',...
                                files.filename_L1A(end-37:end-3),...
                                '_isd','.nc');



%}

ncid = netcdf.create(files.filename_L2,'NETCDF4');

% switch netcdf_type
%     case 'netcdf4'
%         ncid = netcdf.create(files.filename_netCDF,'NETCDF4');
%     case 'netcdf3'
%         ncid = netcdf.create(files.filename_netCDF,'CLASSIC_MODEL');        
% end

long_name_att = 'long_name';
std_name_att = 'standard_name';
calendar_name_att='calendar';
comment_att = 'comment';
units_att = 'units';
scale_factor_att = 'scale_factor';
add_offset_att = 'add_offset';

netcdf_v4_format = 'netcdf4';

% ku_rec_dimension = 'time_l1b_echo_sar_ku';
% nl_dimension = 'Nl';
% nl_dimension_size = N_max_beams_stack_chd;
% ns_dimension = 'Ns';
% space_3D_dimension = 'space_3D';
% space_3D_dimension_size = 3;

ku_rec_dimension = netcdf.defDim(ncid,'time_L2_echo',L2.ISD_num_surfaces_filtered);
single_dimension = netcdf.defDim(ncid,'single',1);

day_units = 'day';
seconds_units = 'seconds';

degrees_units = 'degrees';
meters_units = 'meters';
mps_units = 'meters per second';
dB_units = 'dB';

int16_type = 'NC_SHORT';
int32_type = 'NC_INT';
double_type= 'NC_DOUBLE';
string_type= 'NC_STRING';


%% CODING L2

%--------- Global Attribtues of the netCDF ---------------------------
%{
        %altimeter sensor name
        altimeter_sensor_name='SRAL '; 
        
        % UTC date of the first measurement
        first_meas_time=SPH.product_info.product_time_info.produc_start_time;
        % UTC date of the last measurement
        last_meas_time=SPH.product_info.product_time_info.produc_stop_time;
        
        % Value of the first valid latitude
        first_meas_lat=L2.ISD_lat_surf(1);
        % Value of the last valid latitude
        last_meas_lat=L2.ISD_lat_surf(last);
        % Value of the first valid longitude
        first_meas_lon=L2.ISD_lon_surf(1);
        % Value of the last valid longitude
        last_meas_lon=L2.ISD_lon_surf(last);
        
        % Semi-major axis of the reference ellipsoid
        semi_major_ellipsoid_axis=num2str(semi_major_axis_cst,15);
        
        % Flattening coeffcient of the reference ellipsoid
        ellipsoid_flattening=num2str(flat_coeff_cst,15);    
        

        % Attributes related to the Orbital information required by Porto
        orbit_cycle_num=((SPH.orbit_info.cycle_num));
        orbit_REL_Orbit=((SPH.orbit_info.rel_orbit));
        orbit_ABS_Orbit_Start=((SPH.orbit_info.ABS_Orbit_Start));
%}
%% PACKING L2


%---------- Time variables ----------------------------------------------
time_utc_name = 'time_utc';
id_aux = netcdf.defVar(ncid,time_utc_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,std_name_att,'time');
netcdf.putAtt(ncid,id_aux,long_name_att,'UTC Seconds since 2000-01-01 00:00:00.0+00:00 (Ku-band)');
netcdf.putAtt(ncid,id_aux,calendar_name_att,'Gregorian');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,comment_att,'time at surface of the SAR measurement(multi-looked waveform).');

UTC_day_L2_echo_name = 'UTC_day_L2_echo';
id_aux = netcdf.defVar(ncid,UTC_day_L2_echo_name,int16_type, ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Days since 2000-01-01 00:00:00.0+00:00 (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,day_units);
netcdf.putAtt(ncid,id_aux,comment_att,'days elapsed since 2000-01-01. To be used to link with L1 and L2 records (time_l1b provides the number of seconds since 2000-01-01).');

UTC_sec_L2_echo_name = 'UTC_sec_L2_echo';
id_aux = netcdf.defVar(ncid,UTC_sec_L2_echo_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Seconds in the day UTC, with microsecond resolution (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,comment_att,'seconds in the day. To be used to link L1 and L2 records (time_l1b provides the number of seconds since 2000-01-01).');

L2_record_counter_name = 'L2_record_counter';
id_aux = netcdf.defVar(ncid,L2_record_counter_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,std_name_att,'records');
netcdf.putAtt(ncid,id_aux,long_name_att,'L2 record counter');
netcdf.putAtt(ncid,id_aux,comment_att,'L2 record counter');


%---------- Orbit variables ------------------------------
Latitude_sat_name = 'latitude_sat';
id_aux = netcdf.defVar(ncid,Latitude_sat_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,std_name_att,'latitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'latitude of satellite COM (positive N, negative S) (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Latitude of the satellite Center of Mass at time stamp [-90, +90]: Positive at Nord, Negative at South');

Longitude_sat_name = 'Longitude_sat';
id_aux = netcdf.defVar(ncid,Longitude_sat_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,std_name_att,'longitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'longitude of satellite COM (positive E, negative W) (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Longitude of the satellite Center of Mass at time stamp [-180, +180]: Positive at East, Negative at West');

com_altitude_sat_name = 'com_altitude_sat';
id_aux = netcdf.defVar(ncid,com_altitude_sat_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'altitude of satellite');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Altitude of the satellite Center of Mass at time stamp');

com_altitude_rate_name = 'com_altitude_rate';
id_aux = netcdf.defVar(ncid,com_altitude_rate_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'altitude rate of satellite');
netcdf.putAtt(ncid,id_aux,units_att,mps_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'Altitude rate of the satellite');

com_position_vector_name = 'com_position_vector';
id_aux = netcdf.defVar(ncid,com_position_vector_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'x/y/z coordinates of satellite');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'x/y/z coordinates of the satellite Center of Mass position');

com_position_velocity_name = 'com_position_velocity';
id_aux = netcdf.defVar(ncid,com_position_velocity_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'x/y/z coordinates of satellite velocity');
netcdf.putAtt(ncid,id_aux,units_att,mps_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,700000);
netcdf.putAtt(ncid,id_aux,comment_att,'x/y/z coordinates of the satellite velocity vector');

off_nadir_roll_angle_pf_name = 'off_nadir_roll_angle_pf';
id_aux = netcdf.defVar(ncid,off_nadir_roll_angle_pf_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Platform roll angle mispointing');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Platform roll angle mispointing');

off_nadir_pitch_angle_pf_name = 'off_nadir_pitch_angle_pf';
id_aux = netcdf.defVar(ncid,off_nadir_pitch_angle_pf_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Platform pitch angle mispointing');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Platform pitch angle mispointing');

off_nadir_yaw_angle_pf_name = 'off_nadir_yaw_angle_pf';
id_aux = netcdf.defVar(ncid,off_nadir_yaw_angle_pf_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Platform yaw angle mispointing');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Platform yaw angle mispointing');

%---------- Instrument --------------------------------------
tm_source_seq_counter_name = 'tm_source_seq_counter';
id_aux = netcdf.defVar(ncid,tm_source_seq_counter_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,std_name_att,'records');
netcdf.putAtt(ncid,id_aux,long_name_att,'Instrument sequence counter of the Closest Burst');
netcdf.putAtt(ncid,id_aux,comment_att,'Instrument sequence counter of the Closest Burst');

tm_mode_id_isp_name = 'tm_mode_id_isp';
id_aux = netcdf.defVar(ncid,tm_mode_id_isp_name,string_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Instrument Mode');
netcdf.putAtt(ncid,id_aux,comment_att,'Instrument Mode (SAR/SARIn, RAW/RMC)');

%---------- Geophysical corrections --------------------------------------
wet_tropo_name = 'wet_tropo';
id_aux = netcdf.defVar(ncid,wet_tropo_name,double_type_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Wet Troposphere Correction');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Wet Troposphere Correction');

geoid_name = 'geoid';
id_aux = netcdf.defVar(ncid,geoid_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Geoid EGM2008');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Geoid EGM2008');

mss_name = 'mss';
id_aux = netcdf.defVar(ncid,mss_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Mean Sea Surface from model');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Mean Sea Surface from model');
%---------- Measurement Parameters --------------------------------------
tm_mode_id_isp_nearest_name = 'tm_mode_id_isp_nearest';
id_aux = netcdf.defVar(ncid,tm_mode_id_isp_nearest_name,string_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Instrument Mode of the Closest Burst');
netcdf.putAtt(ncid,id_aux,comment_att,'Instrument Mode (SAR/SARIn, RAW/RMC) of the Closest Burst');

Latitude_surf_name = 'latitude_surf';
id_aux = netcdf.defVar(ncid,Latitude_surf_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,std_name_att,'latitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'latitude of the measurement at time stamp');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'latitude of the measurement at time stamp [-90, +90]: Positive at Nord, Negative at South');

Longitude_surf_name = 'Longitude_surf';
id_aux = netcdf.defVar(ncid,Longitude_surf_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,std_name_att,'longitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'longitude of the measurement at time stamp (positive E, negative W) (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Longitude of the measurement at time stamp [-180, +180]: Positive at East, Negative at West');

%---------- Empirical retrievals --------------------------------------
height_ku_xxxx_name = 'height_ku_xxxx';
id_aux = netcdf.defVar(ncid,height_ku_xxxx_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Height of the surface, ku');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Height of the surface at the measurement point w.r.t. the reference ellipsoid for Ku-band');

height_ka_xxxx_name = 'height_ka_xxxx';
id_aux = netcdf.defVar(ncid,height_ka_xxxx_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Height of the surface, Ka-band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Height of the surface at the measurement point w.r.t. the reference ellipsoid for Ka-band');

back_ku_xxxx_name = 'backs_ku_xxxx';
id_aux = netcdf.defVar(ncid,backs_ku_xxxx_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Backscatter of the surface, ku');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Fully corrected backscatter including instrument gain correction and bias for Ku-band');

back_ka_xxxx_name = 'backs_ka_xxxx';
id_aux = netcdf.defVar(ncid,backs_ka_xxxx_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Backscatter of the surface, ka');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Fully corrected backscatter including instrument gain correction and bias for Ka-band');

%---------- Physical Retrievals --------------------------------------
ssh_ku_name = 'ssh_ku';
id_aux = netcdf.defVar(ncid,ssh_ku_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Sea Surface Height from the ocean retracker for Ku-band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Sea Surface Height from the ocean retracker for Ku-band');

swh_ku_name = 'swh_ku';
id_aux = netcdf.defVar(ncid,swh_ku_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Significant Wave Height from ocean retracker for Ku-band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Significant Wave Height from ocean retracker for Ku-band');

wind_ku_name = 'wind_ku';
id_aux = netcdf.defVar(ncid,wind_ku_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Wind Speed over ocean for Ku-band');
netcdf.putAtt(ncid,id_aux,units_att,mps_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Wind Speed over ocean for Ku-band');

ssh_ka_name = 'ssh_ka';
id_aux = netcdf.defVar(ncid,ssh_ka_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Sea Surface Height from the ocean retracker for Ka-band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Sea Surface Height from the ocean retracker for Ka-band');

swh_ka_name = 'swh_ka';
id_aux = netcdf.defVar(ncid,swh_ka_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Significant Wave Height from ocean retracker for Ka-band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Significant Wave Height from ocean retracker for Ka-band');

wind_ka_name = 'wind_ka';
id_aux = netcdf.defVar(ncid,wind_ka_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Wind Speed over ocean for Ka-band');
netcdf.putAtt(ncid,id_aux,units_att,mps_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Wind Speed over ocean for Ka-band');

%---------- Sea Ice Retrievals --------------------------------
freeboard_name = 'freeboard';
id_aux = netcdf.defVar(ncid,freeboard_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'The difference between the height of the surface of sea ice and the water in open leads');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'The difference between the height of the surface of sea ice and the water in open leads');

snow_depth_name = 'snow_depth';
id_aux = netcdf.defVar(ncid,snow_depth_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Depth of the Snow retrieved with Ku and Ka band measurements');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Depth of the Snow retrieved with Ku and Ka band measurements');

ice_thickness_name = 'ice_thickness';
id_aux = netcdf.defVar(ncid,ice_thickness_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Derived thickness of the sea ice');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Derived thickness of the sea ice');

%----------  Swath Processing Retrievals --------------------------------
bin_index_name = 'bin_index';
id_aux = netcdf.defVar(ncid,bin_index_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Sample index from the L1B product of the retrieved measurement');
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Sample index from the L1B product of the retrieved measurement');

height_swath_name = 'height_swath';
id_aux = netcdf.defVar(ncid,height_swath_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Height of the surface');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Height of the surface at the measurement point w.r.t. the reference ellipsoid');

Latitude_swath_name = 'latitude_swath';
id_aux = netcdf.defVar(ncid,Latitude_swath_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,std_name_att,'latitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'latitude of the measurement');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Latitude of the measurement [-90, +90]: Positive at Nord, Negative at South');

Longitude_swath_name = 'Longitude_swath';
id_aux = netcdf.defVar(ncid,Longitude_swath_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,std_name_att,'longitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'longitude of the measurement');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Longitude of the measurement [-180, +180]: Positive at East, Negative at West');

Unwrapped_AoA_name = 'Unwrapped_AoA';
id_aux = netcdf.defVar(ncid,Unwrapped_AoA_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,std_name_att,'unwrapped_AoA');
netcdf.putAtt(ncid,id_aux,long_name_att,'Angle of Arrival retrieved after the unwrapping process');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'Angle of Arrival retrieved after the unwrapping process');

coherence_name = 'coherence';
id_aux = netcdf.defVar(ncid,coherence_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Sample index from the L1B product of the retrieved measurement');
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Sample index from the L1B product of the retrieved measurement');

dem_diff_name = 'dem_diff';
id_aux = netcdf.defVar(ncid,dem_diff_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Difference w.r.t. the reference Digital Elevation Model');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Difference w.r.t. the reference Digital Elevation Model');

dem_diff_mad_name = 'dem_diff_mad';
id_aux = netcdf.defVar(ncid,dem_diff_mad_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Median Absolute Deviation w.r.t. the reference Digital Elevation Model');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,0);
netcdf.putAtt(ncid,id_aux,add_offset_att,0);
netcdf.putAtt(ncid,id_aux,comment_att,'Median Absolute Deviation w.r.t. the reference Digital Elevation Model');

%{
%---------- Sigma0 --------------------------------
sig0_ice_sheet_L2_echo_name = 'sig0_ice_sheet_L2_echo';
id_aux = netcdf.defVar(ncid,sig0_ice_sheet_L2_echo_name,int16_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Backscatter coefficient sigma0');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'LRM/PLRM modes : ice sheet (CFI) retracking, SAR mode : ice sheet margin retracking. Instrumental corrections included : AGC instrumental errors correction (agc_cor_[x1]_[x2]) and internal calibration correction (sig0_cal_[x1]_[x2])');

sig0_ice_sheet_L2_echo_mean_error_name = 'sig0_ice_sheet_L2_echo_mean_error';
id_aux = netcdf.defVar(ncid,sig0_ice_sheet_L2_echo_mean_error_name,int32_type,single_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Mean error of the difference between the original values and a fitting of these');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,['The fitting has been performed with a sliding window of ',num2str(L2.WINDOW_SIG0),'samples.'])

sig0_ice_sheet_L2_echo_rmse_name = 'sig0_ice_sheet_L2_echo_rmse';
id_aux = netcdf.defVar(ncid,sig0_ice_sheet_L2_echo_rmse_name,int32_type,single_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Root Mean Squared Error of the difference between the original values and a fitting of these');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,['The fitting has been performed with a sliding window of ',num2str(L2.WINDOW_SIG0),'samples.'])

%----------  Elevation of echoing points --------------------------------
elevation_ice_sheet_L2_echo_name = 'elevation_ice_sheet_L2_echo';
id_aux = netcdf.defVar(ncid,elevation_ice_sheet_L2_echo_name,int32_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Elevation of echoing points. Instrumental corrections included. Geophysical corrections included too.');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'LRM/PLRM modes : ice sheet (CFI) retracking, SAR mode : ice sheet margin retracking. Instrumental corrections included : AGC instrumental errors correction (agc_cor_[x1]_[x2]) and internal calibration correction (sig0_cal_[x1]_[x2])');

elevation_ice_sheet_L2_echo_mean_error_name = 'elevation_ice_sheet_L2_echo_mean_error';
id_aux = netcdf.defVar(ncid,elevation_ice_sheet_L2_echo_mean_error_name,int32_type,single_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Mean error of the difference between the original values and a fitting of these');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,['The fitting has been performed with a sliding window of ',num2str(L2.WINDOW_ELEV),'samples.'])

elevation_ice_sheet_L2_echo_rmse_name = 'elevation_ice_sheet_L2_echo_rmse';
id_aux = netcdf.defVar(ncid,elevation_ice_sheet_L2_echo_rmse_name,int32_type,single_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Root Mean Squared Error of the difference between the original values and a fitting of these');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,['The fitting has been performed with a sliding window of ',num2str(L2.WINDOW_ELEV),'samples.'])
%}


%----------  Global Attributes definition -----------------------------------
%---- attributes inherited from Sentinel-3 product description-------------
id_aux = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,id_aux,'creation_time',date_creation);
netcdf.putAtt(ncid,id_aux,'Conventions',netcdf_v4_format);
netcdf.putAtt(ncid,id_aux,'altimeter_sensor_name',altimeter_sensor_name);
netcdf.putAtt(ncid,id_aux,'first_meas_time',first_meas_time);
netcdf.putAtt(ncid,id_aux,'last_meas_time',last_meas_time);
netcdf.putAtt(ncid,id_aux,'first_meas_lat',first_meas_lat);
netcdf.putAtt(ncid,id_aux,'last_meas_lat',last_meas_lat);
netcdf.putAtt(ncid,id_aux,'first_meas_lon',first_meas_lon);
netcdf.putAtt(ncid,id_aux,'last_meas_lon',last_meas_lon);
netcdf.putAtt(ncid,id_aux,'semi_major_ellipsoid_axis',semi_major_ellipsoid_axis);
netcdf.putAtt(ncid,id_aux,'ellipsoid_flattening',ellipsoid_flattening);
%--------------- add the attributes related to intermediate product--------
netcdf.putAtt(ncid,id_aux,'orbit_cycle_num',orbit_cycle_num);
netcdf.putAtt(ncid,id_aux,'orbit_REL_Orbit',orbit_REL_Orbit);
netcdf.putAtt(ncid,id_aux,'orbit_ABS_Orbit_Start',orbit_ABS_Orbit_Start);

netcdf.endDef(ncid);
% time = toc(t6);
% minutes_reading = floor(time/60);
% secs_reading = time - minutes_reading*60;
% disp([num2str(minutes_reading),' minutes and ',num2str(secs_reading),' seconds passed writting L1B']);

end
