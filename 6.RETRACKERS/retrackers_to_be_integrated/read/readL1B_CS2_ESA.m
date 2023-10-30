function data = readL1B_CS2_ESA (filename_L1B)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows for reading altimetry data from L1B ESA products processed 
%
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 13/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       filename_L1B    =   L1B filename with the fullpath information
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
global semi_major_axis_cst flat_coeff_cst;

%% ------------------------------------------------------------------------- 
% Loading data L1B
% ------------------------------------------------------------------------- 
[~,name,ext]=fileparts(filename_L1B);
switch ext
    case '.DBL'
        %% --------------------- Binary file ----------------------------------
        %----------------------------------------------------------------------
        [~,aux]=Cryo_L1b_read(filename_L1B);
        if ~isempty(strfind(name,'SAR'))
            s=size(aux.SAR.data);
        elseif ~isempty(strfind(name,'SIN'))
            s=size(aux.SIN.data);
        end
        
        num_surfaces=s(2)*s(3); %number of total records
        num_db=s(2); %number of records per block
        N_samples=s(1); %number of samples
        data.N_samples=N_samples;
        data.N_records=num_surfaces;
        
        % ---------------------------------------------------------------------
        % GEO: geographical information
        % ---------------------------------------------------------------------
        data.GEO.TAI.days         =   reshape(aux.GEO.TAI.days,[1,num_surfaces]);
        data.GEO.TAI.secs         =   reshape(aux.GEO.TAI.secs,[1,num_surfaces]);
        data.GEO.TAI.microsecs    =   reshape(aux.GEO.TAI.microsecs,[1,num_surfaces]);
        data.GEO.TAI.total        =   data.GEO.TAI.days*24*60*60+data.GEO.TAI.secs+data.GEO.TAI.microsecs*1e-6;
        data.GEO.LAT              =   reshape(aux.GEO.LAT,[1,num_surfaces]);
        data.GEO.LON              =   reshape(aux.GEO.LON,[1,num_surfaces]);
        data.GEO.H_rate           =   reshape(aux.GEO.H_rate,[1,num_surfaces]);
        data.GEO.V                =   reshape(aux.GEO.V.V,[1,num_surfaces]);
        data.GEO.H                =   reshape(aux.GEO.H,[1,num_surfaces]);
        data.GEO.pitch            =   reshape(aux.GEO.Antenna_Bench_Pitch,[1,num_surfaces]).*pi/180;
        data.GEO.roll             =   reshape(aux.GEO.Antenna_Bench_Roll,[1,num_surfaces]).*pi/180;
        data.GEO.yaw              =   reshape(aux.GEO.Antenna_Bench_Yaw,[1,num_surfaces]).*pi/180;
        
        % ---------------------------------------------------------------------
        % MEA: measurements
        % ---------------------------------------------------------------------
        data.MEA.win_delay = reshape(aux.MEA.win_delay,[1,num_surfaces]);
        
        % ---------------------------------------------------------------------
        % COR: Geophysical corrections
        % ---------------------------------------------------------------------
        %-------- load individual corrections ---------------------------------
        % will be replicated for each data block of 20 surfaces
        data.COR.dry_trop                   =   reshape((ones(num_db,1))*(aux.COR.dry_trop),[1,num_surfaces]);
        data.COR.wet_trop                   =   reshape((ones(num_db,1))*(aux.COR.wet_trop),[1,num_surfaces]);
        data.COR.inv_bar                    =   reshape((ones(num_db,1))*(aux.COR.inv_bar),[1,num_surfaces]);
        data.COR.dac                        =   reshape((ones(num_db,1))*(aux.COR.dac),[1,num_surfaces]);
        data.COR.gim_ion                    =   reshape((ones(num_db,1))*(aux.COR.gim_ion),[1,num_surfaces]);
        data.COR.model_ion                  =   reshape((ones(num_db,1))*(aux.COR.model_ion),[1,num_surfaces]);
        data.COR.ocean_equilibrium_tide     =   reshape((ones(num_db,1))*(aux.COR.ocean_equilibrium_tide),[1,num_surfaces]);
        data.COR.ocean_longperiod_tide      =   reshape((ones(num_db,1))*(aux.COR.ocean_longperiod_tide),[1,num_surfaces]);
        data.COR.ocean_loading_tide         =   reshape((ones(num_db,1))*(aux.COR.ocean_loading_tide),[1,num_surfaces]);
        data.COR.solidearth_tide            =   reshape((ones(num_db,1))*(aux.COR.solidearth_tide),[1,num_surfaces]);
        data.COR.geocentric_polar_tide      =   reshape((ones(num_db,1))*(aux.COR.geocentric_polar_tide),[1,num_surfaces]);
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
        data.surf_type_flag=reshape((ones(num_db,1))*(aux.COR.surf_type),[1,num_surfaces]);%{0: open ocean or semi-enclosed seas, 1: enclosed seas or lakes, 2: continental ice, 3: land/transponder}
        
        % ---------------------------------------------------------------------
        % HRM: High-resolution mode: Waveforms
        % ---------------------------------------------------------------------
        if ~isempty(strfind(name,'SAR'))
            % ------------ SAR mode -------------------------------------------
            i2q2_meas=reshape(aux.SAR.data,[N_samples,num_surfaces]);
            scale_power=reshape(aux.SAR.echo_scale_power,[1,num_surfaces]); %known as Scale power: power of 2
            scale_factor=reshape(aux.SAR.echo_scaling,[1,num_surfaces]); %known as Echo scale factor
            %apply scaling factor to waveforms
            data.HRM.power_wav=i2q2_meas.*repmat((2.^(scale_power)).*(1e-9*scale_factor),N_samples,1);
            clear i2q2_meas scale_factor scale_power;
            
            data.HRM.Neff       =   reshape(aux.SAR.N_averaged_echoes,[1,num_surfaces]);
            data.HRM.FLAG.mlQ   =   reshape(aux.SAR.FLAG.Multilook_Incomplete + aux.SAR.FLAG.Beam_Angle_Steering_Err,[1,num_surfaces]); % if 0 no error if 1 a error ocurred in the stack or multilook
            data.HRM.FLAG.pQ    =   reshape(aux.SAR.FLAG.AntiAliased_Power_Echo,[1,num_surfaces]); % if 1 then ok, 0 error in the power
%             data.HRM.ThN        =   zeros(1,num_surfaces);
%             data.HRM.wfm_count  =   1:1:num_surfaces;
        elseif ~isempty(strfind(name,'SIN'))
            % ------------ SARin mode -----------------------------------------
            i2q2_meas=reshape(aux.SIN.data,[N_samples,num_surfaces]);
            scale_power=reshape(aux.SIN.echo_scale_power,[1,num_surfaces]); %known as Scale power: power of 2
            scale_factor=reshape(aux.SIN.echo_scaling,[1,num_surfaces]); %known as Echo scale factor
            %apply scaling factor to waveforms
            data.HRM.power_wav=i2q2_meas.*repmat((2.^(scale_power)).*(1e-9*scale_factor),N_samples,1);
            clear i2q2_meas scale_factor scale_power;
            
            data.HRM.Neff       =   reshape(aux.SIN.N_averaged_echoes,[1,num_surfaces]);
            data.HRM.FLAG.mlQ   =   reshape(aux.SIN.FLAG.Multilook_Incomplete + aux.SIN.FLAG.Beam_Angle_Steering_Err,[1,num_surfaces]); % if 0 no error if 1 a error ocurred in the stack or multilook
            data.HRM.FLAG.pQ    =   reshape(aux.SIN.FLAG.AntiAliased_Power_Echo,[1,num_surfaces]); % if 1 then ok, 0 error in the power
%             data.HRM.ThN        =   zeros(1,num_surfaces);
%             data.HRM.wfm_count  =   1:1:num_surfaces;

        end
       
        %------------------------------------------------------------------
        %-------------------GLOBAL ATTRIBUTES -----------------------------
        %------------------------------------------------------------------
        %Read the header information
        sph=readSPH((filename_L1B));
        % Sensors names and reference files 
        % Data file relatged information
        data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name=deblank(sph.ins_info.ins_conf);
        %gnss sensor name
        if isempty(sph.gnss_info)
            data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name='Not available';
        else
            data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name=sph.gnss_info;
        end
        %doris sensor name
        if isempty(sph.gnss_info)
            data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name='Not available';
        else
            data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name=sph.gnss_info;
        end
        data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name=sph.acq_station;
        data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=sph.product_info.product_time_info.produc_start_time;
        data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=sph.product_info.product_time_info.produc_stop_time;
        
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=[sph.dsds(strncmp(strsplit(strtrim([sph.dsds.ds_name])),'SIRAL_LEVEL_0',length('SIRAL_LEVEL_0'))).filename];
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit=[sph.dsds(strncmp(strsplit(strtrim([sph.dsds.ds_name])),'ORBIT',length('ORBIT'))).filename];
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO=[sph.dsds(strncmp(strsplit(strtrim([sph.dsds.ds_name])),'DORIS',length('DORIS'))).filename];
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1=[sph.dsds(strncmp(strsplit(strtrim([sph.dsds.ds_name])),'CALIBRATION_TYPE_1',length('CALIBRATION_TYPE_1'))).filename];
        idx=find(strncmp(strsplit(strtrim([sph.dsds.ds_name])),'CALIBRATION_TYPE_2',length('CALIBRATION_TYPE_2'))~=0);
        if isempty(idx)
            data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2='Not available (not applied in the product)';
        else
            data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2=[sph.dsds(idx).filename];
        end
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2='Not available for CR2';
        data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation=[sph.dsds(strncmp(strsplit(strtrim([sph.dsds.ds_name])),'IPF_RA_DATABASE',length('IPF_RA_DATABASE'))).filename];
        data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=num2str(semi_major_axis_cst,15);
        data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=num2str(flat_coeff_cst,15); 
        % Orbit related information     
        data.GLOBAL_ATT.ORBIT_INFO.orbit_phase_code=sph.orbit_info.phase_code;
        data.GLOBAL_ATT.ORBIT_INFO.orbit_cycle_num=sph.orbit_info.cycle_num;
        data.GLOBAL_ATT.ORBIT_INFO.orbit_REL_Orbit=sph.orbit_info.rel_orbit;
        data.GLOBAL_ATT.ORBIT_INFO.orbit_ABS_Orbit_Start=sph.orbit_info.ABS_Orbit_Start;
        data.GLOBAL_ATT.ORBIT_INFO.orbit_Rel_Time_ASC_Node_Start=sph.orbit_info.Rel_Time_ASC_Node_Start;
        data.GLOBAL_ATT.ORBIT_INFO.orbit_ABS_Orbit_Stop=sph.orbit_info.ABS_Orbit_Stop;
        data.GLOBAL_ATT.ORBIT_INFO.orbit_Rel_Time_ASC_Node_Stop=sph.orbit_info.Rel_Time_ASC_Node_Stop;
        data.GLOBAL_ATT.ORBIT_INFO.orbit_Equator_Cross_Time=sph.orbit_info.Equator_Cross_Time;
        data.GLOBAL_ATT.ORBIT_INFO.orbit_Equator_Cross_Long=sph.orbit_info.Equator_Cross_Long;
        data.GLOBAL_ATT.ORBIT_INFO.orbit_Ascending_Flag=sph.orbit_info.Ascending_Flag;
        
        %COMPUTE THE SIGMA0 SCALING FACTOR
        [data]   = sigma0_scaling_factor (data);
        
        
        
    case '.nc'
        %% --------------------- NetCDF ---------------------------------------
        % ESA based approach & naming convetion
        % Further updates might be required for NetCDF from ESA
        
        
        
    otherwise
        error(strcat('File extension ',cnf_p.mission,' is not currently contemplated or not valid'));
end

end

