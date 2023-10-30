function data = readL1B_S6_ISD (filename_L1B,cnf_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows for reading altimetry data from L1B Sentinel-6 products
% processed by ISD
%
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 15/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       MANDATORY:
%           -filename_L1B    =   L1B filename with the fullpath information
%           -cnf_p           =   Structure with configuration parameters of the
%                                L2 processor
%       OPTIONAL:
%           -filename_L1BS   =   L1B-S filename with the fullpath name
%       
% 
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
% Sentinel-3 format and within the SEOMs projects. It is considered only
% for the final netCDF L1B products as defined in the PSF issue 1.3 
%(not all fields considered in SEOMS Sentinel-3 are available in the Sentinel-6 product format netcdf)

%% ------------------ Global Variables ------------------------------------
global c_cst sec_in_day_cst;

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
name='HR';
switch ext
    case '.nc'
        %% --------------------- NetCDF ---------------------------------------
        % Based on the L1B product generated for DeDop a la Sentinel-3
        % open the netCDF file
        ncid=netcdf.open(filename_L1B,'NC_NOWRITE');
        
        % Get dimensions
        dimid=netcdf.inqDimID(ncid,'time_20_ku');
        [~,num_surfaces]=netcdf.inqDim(ncid,dimid);
        dimid=netcdf.inqDimID(ncid,'ns');
        [~,N_samples]=netcdf.inqDim(ncid,dimid);
        dimid=netcdf.inqDimID(ncid,'nl');
        [~,N_looks_conf_stack]=netcdf.inqDim(ncid,dimid); %total number of looks configured for a stack
        netcdf.close(ncid);
        
        data.N_samples=N_samples;
        data.N_records=num_surfaces;
        
        % -----------------------------------------------------------------
        % GEO: geographical information
        % -----------------------------------------------------------------

        data.GEO.TAI.total             =   ncread(filename_L1B,'time_tai_20_ku').';
        data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/sec_in_day_cst);
        data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*sec_in_day_cst);
        data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
        data.GEO.LAT                   =   ncread(filename_L1B,'latitude_20_ku').';
        data.GEO.LON                   =   ncread(filename_L1B,'longitude_20_ku').'; %defined between [0,360]
        data.GEO.LON(data.GEO.LON>180) =   data.GEO.LON(data.GEO.LON>180)-360.0; 
        data.GEO.H_rate                =   ncread(filename_L1B,'com_altitude_rate_20_ku').'; % m/s
        velocity=double(ncread(filename_L1B,'com_velocity_vector_20_ku')).'; %vector of vx,vx,vz
        data.GEO.V                     =   sqrt(velocity(:,1).^2+velocity(:,2).^2+velocity(:,3).^2).';
        clear velocity;
        data.GEO.H                     =   ncread(filename_L1B,'com_altitude_20_ku').';
        data.GEO.pitch                 =   (ncread(filename_L1B,'off_nadir_pitch_angle_pf_20_ku').').*pi/180;
        data.GEO.roll                  =   (ncread(filename_L1B,'off_nadir_roll_angle_pf_20_ku').').*pi/180;
        data.GEO.yaw                   =   (ncread(filename_L1B,'off_nadir_yaw_angle_pf_20_ku').').*pi/180;
        
        
        % -----------------------------------------------------------------
        % MEA: measurements
        % -----------------------------------------------------------------
        data.MEA.win_delay = double(ncread(filename_L1B,'altimeter_range_calibrated_20_ku')).'*2.0/c_cst;
        
        % ---------------------------------------------------------------------
        % COR: Geophysical corrections
        % ---------------------------------------------------------------------
%         % Currently is not foreseen to be included in the L1B-product of
%         % Sentinel-6
%         data.COR.surf_type_flag=ncread(filename_L1B,'flag_surface_classification_20_ku ').'; %{0: open ocean, 1: land, 2: continental_water, 3: acquatic vegetation, 4: continental_ice_snow, 5: floating_ice, 6: salted_basin}
%         
        % -----------------------------------------------------------------
        % HRM: High-resolution mode: Waveforms
        % -----------------------------------------------------------------
        % ------------------ Waveforms ------------------------------------
        if ~isempty(strfind(name,'HR'))
            % ------------ SAR mode ---------------------------------------
            i2q2_meas=double(ncread(filename_L1B,'hr_power_waveform_20_ku')).';
            scale_factor=ncread(filename_L1B,'waveform_scale_factor_20_ku');
            %apply scaling factor to waveforms
            data.HRM.power_wav=(i2q2_meas.*repmat(scale_factor,1,N_samples)).';
            clear i2q2_meas scale_factor;
            
            
            data.HRM.Neff       =   ncread(filename_L1B,'num_looks_start_stop_20_ku').'; %effective number of beams that form the stack including possuible looks that are set entirely to zero
            data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
            data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
            %             data.HRM.ThN        =   zeros(1,num_surfaces);
            %             data.HRM.wfm_count  =   1:1:num_surfaces;
            
            
            % ----sigma0 scaling factor for conversion from Power to sigma0
            %units
            data.HRM.s0_sf=ncread(filename_L1B,'sigma0_scaling_factor_20_ku').';
            
            
            %--------------------------------------------------------------
            %------------- Stack characterization parameters --------------
            %--------------------------------------------------------------
            %----- Dopppler mask ------------------------------------------
            data.HRM.Doppler_mask   =   ncread(filename_L1B,'stack_mask_start_stop_20_ku')+1; %due to the way the stack mask vector is saved in the netcdf the mask is not exactly the same as original
            % internally when using the ncread function the scaling factor
            % of ZP is considered
            data.HRM.pri_surf=ncread(filename_L1B,'pulse_repetition_interval_20_ku').'; %in seconds
            
            
            %------------- Geometry-related parameters ------------------------
            if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
                switch cnf_p.looks_index_method
                    case {'Doppler'}
                        %This option Does not provide the correct information but it is kept for comparison purposes
                        switch cnf_p.fd_method
                            case 'exact'
                                % exact method would require precise information for
                                % each beam within the stack (not currently developed for L1B products for SEOMs)
                                %considering the definition of Doppler angle in
                                %DPM which is different from PSF S6 issue 1.3 Beam_ang= 90-Dopp_ang
                                data.HRM.beam_ang_stack   =  pi/2+ncread(filename_L1BS,'doppler_angle_20_ku')-ncread(filename_L1BS,'look_angle_20_ku'); %TBD in radians
                                % it is not foreseen to provide the specific
                                % PRI/velocity for each look coming from a different
                                % burst: use the same pri for that surface
                                data.HRM.pri_stack        =  ones(N_looks_conf_stack,1)*(ncread(filename_L1B,'pulse_repetition_interval_20_ku').'); %TBD %eventually for sentinel-6 the PRF could change
                                data.GEO.V_stack          =  ones(N_looks_conf_stack,1)*data.GEO.V;
                            case 'approximate'
                                data.HRM.pointing_ang_start_surf=(ncread(filename_L1B,'pointing_angle_start_20_ku').'); % in radians
                                data.HRM.pointing_ang_stop_surf=(ncread(filename_L1B,'pointing_angle_stop_20_ku').'); % in radians
                                data.HRM.doppler_angle_start_surf=(ncread(filename_L1B,'doppler_angle_start_20_ku').'); % in radians
                                data.HRM.doppler_angle_stop_surf=(ncread(filename_L1B,'doppler_angle_stop_20_ku').'); % in radians
                        end
                    case {'Look_angle'}
                        switch cnf_p.look_ang_method
                            case 'exact'
                                % exact method would require precise information for
                                % each beam within the stack (not currently developed for L1B products for SEOMs)
                                data.HRM.look_ang_stack   =  ncread(filename_L1B,'look_angle_20_ku'); % in radians
                                % it is not foreseen to provide the specific
                                % PRI/velocity for each look coming from a different
                                % burst: use the same pri for that surface
                                data.HRM.pri_stack        =  ones(N_looks_conf_stack,1)*(ncread(filename_L1B,'pulse_repetition_interval_20_ku').'); %TBD %eventually for sentinel-6 the PRF could change
                                data.GEO.V_stack          =  ones(N_looks_conf_stack,1)*data.GEO.V; % TBD norm of velocity for each beam within stack
                            case 'approximate'
                                data.HRM.look_ang_start_surf=(ncread(filename_L1B,'look_angle_start_20_ku').'); % in radians
                                data.HRM.look_ang_stop_surf=(ncread(filename_L1B,'look_angle_stop_20_ku').'); % in radians
                        end
                        
                end
            end
            
            
            
        elseif ~isempty(strfind(name,'LR'))
            % ----------- LR: low resolution data extraction --------------
        end
        
%         %------------------------------------------------------------------
%         %-------------------GLOBAL ATTRIBUTES -----------------------------
%         %------------------------------------------------------------------
%         % Data file relatged information
%         data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name=ncreadatt(filename_L1B,'/','altimeter_name');
%         data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name='Not available';
%         data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name='Not available';      
%         data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name='Not available';
%         data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=ncreadatt(filename_L1B,'/','first_meas_time');
%         data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=ncreadatt(filename_L1B,'/','last_meas_time');        
%         data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=ncreadatt(filename_L1B,'/','first_meas_time ');
%         data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit=ncreadatt(filename_L1B,'/','input_file_orbit');
%         data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO=ncreadatt(filename_L1B,'/','input_file_uso_drift');
%         data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1=strcat(ncreadatt(filename_L1B,'/','input_file_range_corr_internal_delay_cal_ku'),'; ',...
%                                                                    ncreadatt(filename_L1B,'/','input_file_altimeter_power_drift_ku'),'; ',...
%                                                                    ncreadatt(filename_L1B,'/','input_file_attenuator_table_ku'));
%         data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2=ncreadatt(filename_L1B,'/','input_file_transfer_function_ku ');
%         data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2=ncreadatt(filename_L1B,'/','input_file_transfer_function_c');
%         data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation=strcat(ncreadatt(filename_L1B,'/','input_file_characterization'),'; ',...
%                                                                        ncreadatt(filename_L1B,'/','input_file_characterization_array'));
%         data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=ncreadatt(filename_L1B,'/','semi_major_ellipsoid_axis');
%         data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=ncreadatt(filename_L1B,'/','ellipsoid_flattening');
%         % Orbit related information     
%         data.GLOBAL_ATT.ORBIT_INFO.orbit_phase_code='Not Available';
%         data.GLOBAL_ATT.ORBIT_INFO.orbit_cycle_num=ncreadatt(filename_L1B,'/','cycle_number');
%         data.GLOBAL_ATT.ORBIT_INFO.orbit_REL_Orbit=strcat('(Pass number in Cycle) ',ncreadatt(filename_L1B,'/','pass_number'));
%         data.GLOBAL_ATT.ORBIT_INFO.orbit_ABS_Orbit_Start=ncreadatt(filename_L1B,'/','first_meas_time');
%         data.GLOBAL_ATT.ORBIT_INFO.orbit_Rel_Time_ASC_Node_Start='Not Available';
%         data.GLOBAL_ATT.ORBIT_INFO.orbit_ABS_Orbit_Stop=ncreadatt(filename_L1B,'/','last_meas_time');
%         data.GLOBAL_ATT.ORBIT_INFO.orbit_Rel_Time_ASC_Node_Stop='Not Available';
%         data.GLOBAL_ATT.ORBIT_INFO.orbit_Equator_Cross_Time=ncreadatt(filename_L1B,'/','equator_time');
%         data.GLOBAL_ATT.ORBIT_INFO.Equator_Cross_Long=ncreadatt(filename_L1B,'/','equator_longitude');
%         data.GLOBAL_ATT.ORBIT_INFO.orbit_Ascending_Flag='Not available';

    case '.mat'
        %% ------------------------ Matlab data -------------------------------
        load(filename_L1B);
        s=size(wfm_cor_i2q2_sar_ku);
        N_samples=s(2);
        num_surfaces=s(1);
        data.N_samples=N_samples;
        data.N_records=num_surfaces;
        
        % -----------------------------------------------------------------
        % GEO: geographical information
        % -----------------------------------------------------------------
        data.GEO.TAI.total             =   time_surf;
        data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/sec_in_day_cst);
        data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*sec_in_day_cst);
        data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
        data.GEO.LAT                   =   lat_sat; %degrees
        data.GEO.LON                   =   lon_sat; %degrees defined between [-180,180]
        data.GEO.H_rate                =   alt_rate_sat; % m/s
        data.GEO.V                     =   sqrt(x_vel_sat.^2+y_vel_sat.^2+z_vel_sat.^2).'; %m/s
        data.GEO.H                     =   alt_sat; % m
        data.GEO.pitch                 =   pitch_sar; % rad
        data.GEO.roll                  =   roll_sar; % rad
        data.GEO.yaw                   =   yaw_sar; % rad
        
        % -----------------------------------------------------------------
        % MEA: measurements
        % -----------------------------------------------------------------
        data.MEA.win_delay = win_delay_surf;
        
        % ---------------------------------------------------------------------
        % COR: Geophysical corrections
        % ---------------------------------------------------------------------
        % Currently is not foreseen to be included in the L1B-product of
        % Sentinel-6
        
                % -----------------------------------------------------------------
        % HRM: High-resolution mode: Waveforms
        % -----------------------------------------------------------------
        % ------------------ Waveforms ------------------------------------
        if ~isempty(strfind(name,'HR'))
            % ------------ SAR mode ---------------------------------------
            data.HRM.power_wav=wfm_cor_i2q2_sar_ku.';
                        
            
            data.HRM.Neff       =   N_beams_start_stop; %effective number of beams that form the stack including possuible looks that are set entirely to zero
            data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
            data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
            %             data.HRM.ThN        =   zeros(1,num_surfaces);
            %             data.HRM.wfm_count  =   1:1:num_surfaces;
            
            
            % ----sigma0 scaling factor for conversion from Power to sigma0
            %units
            data.HRM.s0_sf=wfm_scaling_factor_sar_ku; %dB
            
            
            %--------------------------------------------------------------
            %------------- Stack characterization parameters --------------
            %--------------------------------------------------------------
            %----- Dopppler mask ------------------------------------------
            data.HRM.Doppler_mask   =   stack_mask_vector;
            data.HRM.pri_surf=pri_sar_nadir_surf; %in seconds
            %------------- Geometry-related parameters ------------------------
            switch cnf_p.looks_index_method
                case {'Doppler'}
                    %This option Does not provide the correct information but it is kept for comparison purposes
                    switch cnf_p.fd_method
                        case 'exact'
                            % exact method would require precise information for
                            % each beam within the stack (not currently developed for L1B products for SEOMs)
                            data.HRM.beam_ang_stack   =  beam_ang_surf.'; %TBD in radians
                            data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
                            data.HRM.V_stack          =  vel_norm_sat_beam_surf.';
                        case 'approximate'
                            data.HRM.pointing_ang_start_surf=start_pointing_angle; % in radians
                            data.HRM.pointing_ang_stop_surf=stop_pointing_angle; % in radians
                            data.HRM.doppler_angle_start_surf=start_doppler_angle; % in radians
                            data.HRM.doppler_angle_stop_surf=stop_doppler_angle; % in radians
                    end
                case {'Look_angle'}
                    switch cnf_p.look_ang_method
                        case 'exact'
                            % exact method would require precise information for
                            % each beam within the stack (not currently developed for L1B products for SEOMs)
                            data.HRM.look_ang_stack   =  look_ang_surf.'; % in radians
                            data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
                            data.HRM.V_stack          =  vel_norm_sat_beam_surf.'; % TBD norm of velocity for each beam within stack
                        case 'approximate'
                            data.HRM.look_ang_start_surf=start_look_angle; % in radians
                            data.HRM.look_ang_stop_surf=stop_look_angle; % in radians
                    end
                    
            end
            
            
            
        elseif ~isempty(strfind(name,'LR'))
            % ----------- LR: low resolution data extraction --------------
        end
        
        
        
        
    otherwise
        error(strcat('File extension ',cnf_p.mission,' is not currently contemplated or not valid'));
end

end

