function [data] = readL1B_CS2_ISD (filename_L1B,cnf_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows for reading altimetry data from L1B ISD products processed by isardSAT 
%
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 14/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       MANDATORY:
%           -filename_L1B    =   L1B filename with the fullpath information
%           -cnf_p           =   Structure with configuration parameters of the
%                           L2 processor
%       OPTIONAL:
%           -filename_L1BS   =   L1B-S filename with the fullpath
%           information
% OUTPUT:
%       data        =   structure of data as defined by our L2 processor
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS:
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: Based on read_alt_data.m
% v1.1: Include global attributes sensor name, acq. station, acq date,... for
% the reference L1B product & orbital parameter info in a similar fashion
% as the global attribtues provided for output L1B product according to
% Sentinel-3 format and within the SEOMs projects

%% ------------------ Global Variables ------------------------------------
global c_cst;
global prf_chd;

%% ---------------- Handling input variables ------------------------------
if(nargin<2 || nargin>(2+1*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('filename_L1BS',{''},@(x)ischar(x));
p.parse(varargin{:});
filename_L1BS=char(p.Results.filename_L1BS);
clear p;

%% ------------------------------------------------------------------------- 
% Loading data L1B
% ------------------------------------------------------------------------- 
[~,name,ext]=fileparts(filename_L1B);
switch ext
    case '.nc'
        %% --------------------- NetCDF ---------------------------------------
        % Based on the L1B product generated for DeDop a la Sentinel-3
        % open the netCDF file
        ncid=netcdf.open(filename_L1B,'NC_NOWRITE');
        
        % Get dimensions
        dimid=netcdf.inqDimID(ncid,'time_l1b_echo_sar_ku');
        [~,num_surfaces]=netcdf.inqDim(ncid,dimid);
        dimid=netcdf.inqDimID(ncid,'echo_sample_ind');
        [~,N_samples]=netcdf.inqDim(ncid,dimid);
        netcdf.close(ncid);
        
        data.N_samples=N_samples;
        data.N_records=num_surfaces;
        
        % -----------------------------------------------------------------
        % GEO: geographical information
        % -----------------------------------------------------------------
        data.GEO.TAI.days         =   ncread(filename_L1B,'UTC_day_l1b_echo_sar_ku').';
        aux=ncread(filename_L1B,'UTC_sec_l1b_echo_sar_ku');
        data.GEO.TAI.secs         =   floor(aux).';
        data.GEO.TAI.microsecs    =   ((aux-floor(aux))*1e6).';
        data.GEO.TAI.total        =   ncread(filename_L1B,'time_l1b_echo_sar_ku').';
        data.GEO.LAT              =   ncread(filename_L1B,'lat_l1b_echo_sar_ku').';
        data.GEO.LON              =   ncread(filename_L1B,'lon_l1b_echo_sar_ku').';
        data.GEO.H_rate           =   ncread(filename_L1B,'orb_alt_rate_l1b_echo_sar_ku').';
        vx                        =   ncread(filename_L1B,'x_vel_l1b_echo_sar_ku').';
        vy                        =   ncread(filename_L1B,'y_vel_l1b_echo_sar_ku').';
        vz                        =   ncread(filename_L1B,'z_vel_l1b_echo_sar_ku').';
        data.GEO.V                =   sqrt(vx.^2+vy.^2+vz.^2);
        clear vx vy vz;
        data.GEO.H                =   ncread(filename_L1B,'alt_l1b_echo_sar_ku').';
        attitude=double(ncread(filename_L1B,'satellite_mispointing_l1b_sar_echo_ku')).';
        data.GEO.pitch            =   (attitude(:,1).').*pi/180;
        data.GEO.roll             =   (attitude(:,2).').*pi/180;
        data.GEO.yaw              =   (attitude(:,3).').*pi/180;
        clear attitude;
        
        % -----------------------------------------------------------------
        % MEA: measurements
        % -----------------------------------------------------------------
        data.MEA.win_delay = double(ncread(filename_L1B,'range_ku_l1b_echo_sar_ku')).'*2.0/c_cst;
        
        % ---------------------------------------------------------------------
        % COR: Geophysical corrections
        % ---------------------------------------------------------------------
        %-------- load individual corrections ---------------------------------
        data.COR.dry_trop                   =   double(ncread(filename_L1B,'dry_tropo_correction_l1b_echo_sar_ku')).';
        data.COR.wet_trop                   =   double(ncread(filename_L1B,'wet_tropo_correction_l1b_echo_sar_ku')).';
        data.COR.inv_bar                    =   double(ncread(filename_L1B,'inverse_baro_correction_l1b_echo_sar_ku')).';
        data.COR.dac                        =   double(ncread(filename_L1B,'Dynamic_atmospheric_correction_l1b_echo_sar_ku')).';
        data.COR.gim_ion                    =   double(ncread(filename_L1B,'GIM_iono_correction_l1b_echo_sar_ku')).';
        data.COR.model_ion                  =   double(ncread(filename_L1B,'model_iono_correction_l1b_echo_sar_ku')).';
        data.COR.ocean_equilibrium_tide     =   double(ncread(filename_L1B,'ocean_equilibrium_tide_l1b_echo_sar_ku')).';
        data.COR.ocean_longperiod_tide      =   double(ncread(filename_L1B,'long_period_tide_l1b_echo_sar_ku')).';
        data.COR.ocean_loading_tide         =   double(ncread(filename_L1B,'ocean_loading_tide_l1b_echo_sar_ku')).';
        data.COR.solidearth_tide            =   double(ncread(filename_L1B,'solid_earth_tide_l1b_echo_sar_ku')).';
        data.COR.geocentric_polar_tide      =   double(ncread(filename_L1B,'geocentric_polar_tide_l1b_echo_sar_ku')).';
        %---------- Combined corrections --------------------------------------
        data.COR.prop_GIM_ion   =   data.COR.dry_trop + data.COR.wet_trop + data.COR.gim_ion; % not clear??
        data.COR.prop_Mod_ion   =   data.COR.dry_trop + data.COR.wet_trop + data.COR.model_ion; % not clear ??
        data.COR.surf_dac       =   data.COR.dac; % SSB shoudl be added according to Cristina's code??
        data.COR.surf_invb      =   data.COR.inv_bar; % SSB shoudl be added according to Cristina's code??
        data.COR.geop           =   data.COR.ocean_equilibrium_tide + data.COR.ocean_longperiod_tide...
            + data.COR.ocean_loading_tide + data.COR.solidearth_tide + data.COR.geocentric_polar_tide;
        
        % -----------------------------------------------------------------
        % SURFACE TYPE INFORMATION
        % -----------------------------------------------------------------
        data.surf_type_flag=ncread(filename_L1B,'surf_type_l1b_echo_sar_ku').'; %{0: open ocean or semi-enclosed seas, 1: enclosed seas or lakes, 2: continental ice, 3: land/transponder}
        
        
        % -----------------------------------------------------------------
        % HRM: High-resolution mode: Waveforms
        % -----------------------------------------------------------------
        % ------------------ Waveforms ------------------------------------
        switch cnf_p.mode
            case {'SAR','HR'}
                % ------------ SAR mode ---------------------------------------
                i2q2_meas=double(ncread(filename_L1B,'i2q2_meas_ku_l1b_echo_sar_ku'));
                scale_factor=double(ncread(filename_L1B,'waveform_scale_factor_l1b_echo_sar_ku'));
                %apply scaling factor to waveforms
                data.HRM.power_wav=i2q2_meas.*repmat(scale_factor.',N_samples,1);
                clear i2q2_meas scale_factor;
                
                % need to clarify the use of these parameters
                data.HRM.Neff       =   ncread(filename_L1B,'nb_stack_start_stop_l1b_echo_sar_ku').'; %effective number of beams that form the stack
                data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
                data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
                %             data.HRM.ThN        =   zeros(1,num_surfaces);
                %             data.HRM.wfm_count  =   1:1:num_surfaces;
                
                
                % ----sigma0 scaling factor for conversion from Power to sigma0
                %units
                data.HRM.s0_sf=ncread(filename_L1B,'scale_factor_ku_l1b_echo_sar_ku').';
            case {'SARIN','SARin','SIN'}
                % ------------ SARin mode -------------------------------------
                %TBD with albert how the L1B product is defined in terms of
                %waveforms for the SARin mode
        end
        
        %------------------------------------------------------------------
        %------------- Stack characterization parameters ------------------
        %------------------------------------------------------------------
        
        %----- Dopppler mask ----------------------------------------------        
        data.HRM.Doppler_mask   =   ncread(filename_L1B,'stack_mask_range_bin_l1b_echo_sar_ku')+1; %due to the way the stack mask vector is saved in the netcdf the mask is not exactly the same as original
        scale_factor = ncreadatt(filename_L1B,'pri_lrm_l1b_echo_sar_ku','scale_factor');
        data.HRM.pri_surf=ones(1,num_surfaces).*1.0/prf_chd;%(double(ncread(filename_L1B,'pri_lrm_l1b_echo_sar_ku')).*scale_factor).'; %in seconds
        %the scaling factor is not implemented automatically for variables
        %saved as int64 or long issue in reading such data with matlab in
        %panoply is perfectly read
        %------------- Geometry-related parameters ------------------------
        if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
            switch cnf_p.looks_index_method
                case {'Doppler'}
                    %This option Does not provide the correct information but it is kept for comparison purposes
                    switch cnf_p.fd_method
                        case 'exact'
                            % exact method would require precise information for
                            % each beam within the stack (not currently developed for L1B products for SEOMs)
                            data.HRM.beam_ang_stack   =  ncread(filename_L1BS,'beam_angle_stack_l1b_echo_sar_ku'); %TBD in radians
                            data.HRM.pri_stack        =  ncread(filename_L1BS,'pri_stack_l1b_echo_sar_ku'); %TBD %eventually for sentinel-6 the PRF could change
                            data.HRM.V_stack          =  ncread(filename_L1BS,'v_norm_stack_l1b_echo_sar_ku'); % TBD norm of velocity for each beam within stack
                        case 'approximate'
                            data.HRM.pointing_ang_start_surf=(ncread(filename_L1B,'pointing_angle_start_l1b_echo_sar_ku').'); % in radians
                            data.HRM.pointing_ang_stop_surf=(ncread(filename_L1B,'pointing_angle_stop_l1b_echo_sar_ku').'); % in radians
                            data.HRM.doppler_angle_start_surf=(ncread(filename_L1B,'doppler_angle_start_l1b_echo_sar_ku').'); % in radians
                            data.HRM.doppler_angle_stop_surf=(ncread(filename_L1B,'doppler_angle_stop_l1b_echo_sar_ku').'); % in radians
                    end
                case {'Look_angle'}
                    switch cnf_p.look_ang_method
                        case 'exact'
                            % exact method would require precise information for
                            % each beam within the stack (not currently developed for L1B products for SEOMs)
                            data.HRM.look_ang_stack   =  ncread(filename_L1BS,'look_angle_stack_l1b_echo_sar_ku'); %TBD in radians
                            data.HRM.pri_stack        =  ncread(filename_L1BS,'pri_stack_l1b_echo_sar_ku'); %TBD %eventually for sentinel-6 the PRF could change
                            data.HRM.V_stack          =  ncread(filename_L1BS,'v_norm_stack_l1b_echo_sar_ku'); % TBD norm of velocity for each beam within stack
                        case 'approximate'
                            data.HRM.look_ang_start_surf=(ncread(filename_L1B,'look_angle_start_l1b_echo_sar_ku').'); % in radians
                            data.HRM.look_ang_stop_surf=(ncread(filename_L1B,'look_angle_stop_l1b_echo_sar_ku').'); % in radians
                    end
                    
            end
        end
        
        %for checking purposes
        data.HRM.pointing_ang_start_surf=(ncread(filename_L1B,'pointing_angle_start_l1b_echo_sar_ku').'); % in radians
        data.HRM.pointing_ang_stop_surf=(ncread(filename_L1B,'pointing_angle_stop_l1b_echo_sar_ku').'); % in radians
        data.HRM.doppler_angle_start_surf=(ncread(filename_L1B,'doppler_angle_start_l1b_echo_sar_ku').'); % in radians
        data.HRM.doppler_angle_stop_surf=(ncread(filename_L1B,'doppler_angle_stop_l1b_echo_sar_ku').'); % in radians
        
        data.HRM.look_ang_start_surf=(ncread(filename_L1B,'look_angle_start_l1b_echo_sar_ku').'); % in radians
        data.HRM.look_ang_stop_surf=(ncread(filename_L1B,'look_angle_stop_l1b_echo_sar_ku').'); % in radians
        
        %------------------------------------------------------------------
        %-------------------GLOBAL ATTRIBUTES -----------------------------
        %------------------------------------------------------------------
        % Data file relatged information
        data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name=ncreadatt(filename_L1B,'/','altimeter_sensor_name');
        data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name=ncreadatt(filename_L1B,'/','gnss_sensor_name');
        data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name=ncreadatt(filename_L1B,'/','doris_sensor_name');      
        data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name=ncreadatt(filename_L1B,'/','acq_station_name');
        data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=ncreadatt(filename_L1B,'/','first_meas_time');
        data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=ncreadatt(filename_L1B,'/','last_meas_time');        
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=ncreadatt(filename_L1B,'/','xref_altimeter_level0');
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit=ncreadatt(filename_L1B,'/','xref_altimeter_orbit');
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO=ncreadatt(filename_L1B,'/','xref_doris_USO');
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1=ncreadatt(filename_L1B,'/','xref_altimeter_ltm_sar_cal1');
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2=ncreadatt(filename_L1B,'/','xref_altimeter_ltm_ku_cal2');
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2=ncreadatt(filename_L1B,'/','xref_altimeter_ltm_c_cal2');
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation=ncreadatt(filename_L1B,'/','xref_altimeter_characterisation');
        data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=ncreadatt(filename_L1B,'/','semi_major_ellipsoid_axis');
        data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=ncreadatt(filename_L1B,'/','ellipsoid_flattening');
        % Orbit related information     
        data.GLOBAL_ATT.ORBIT_INFO.orbit_phase_code=ncreadatt(filename_L1B,'/','orbit_phase_code');
        data.GLOBAL_ATT.ORBIT_INFO.orbit_cycle_num=ncreadatt(filename_L1B,'/','orbit_cycle_num');
        data.GLOBAL_ATT.ORBIT_INFO.orbit_REL_Orbit=ncreadatt(filename_L1B,'/','orbit_REL_Orbit');
        data.GLOBAL_ATT.ORBIT_INFO.orbit_ABS_Orbit_Start=ncreadatt(filename_L1B,'/','orbit_ABS_Orbit_Start');
        data.GLOBAL_ATT.ORBIT_INFO.orbit_Rel_Time_ASC_Node_Start=ncreadatt(filename_L1B,'/','orbit_Rel_Time_ASC_Node_Start');
        data.GLOBAL_ATT.ORBIT_INFO.orbit_ABS_Orbit_Stop=ncreadatt(filename_L1B,'/','orbit_ABS_Orbit_Stop');
        data.GLOBAL_ATT.ORBIT_INFO.orbit_Rel_Time_ASC_Node_Stop=ncreadatt(filename_L1B,'/','orbit_Rel_Time_ASC_Node_Stop');
        data.GLOBAL_ATT.ORBIT_INFO.orbit_Equator_Cross_Time=ncreadatt(filename_L1B,'/','orbit_Equator_Cross_Time');
        data.GLOBAL_ATT.ORBIT_INFO.orbit_Equator_Cross_Long=ncreadatt(filename_L1B,'/','orbit_Equator_Cross_Long');
        data.GLOBAL_ATT.ORBIT_INFO.orbit_Ascending_Flag=ncreadatt(filename_L1B,'/','orbit_Ascending_Flag');
        
        
        
    otherwise
        error(strcat('File extension ',cnf_p.mission,' is not currently contemplated or not valid'));
end

end

