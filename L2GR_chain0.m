function filesBulk = L2GR_chain0(filesBulk, cnf, chd, cst)

%% 2. Input Files

% 2.1 L1B Products
filesBulkFields = fieldnames(filesBulk);
aux=struct2cell(filesBulk.auxFiles); aux=aux(1,:);
filesBulk.CNF_file_L2=[filesBulk.auxPath filesBulk.auxFiles(~cellfun(@isempty,strfind(aux,'cnf_file_L2'))).name];
run(filesBulk.CNF_file_L2); % cnf struct output

for i_file = 1:length(filesBulkFields)
    L1B=NaN;
    L2 = [];
    if strcmp('filename_L1B_LR',filesBulkFields(i_file))
        [L1B]= reading_L1B_LR (filesBulk, chd, cnf, cst); % Not implemented
        L1B.product_type='L1B_LR';
        display('L1B_LR found')
    elseif strcmp('filename_L1B_LRRMC',filesBulkFields(i_file))
        [L1B]= reading_L1B_LRRMC (filesBulk, chd, cnf, cst); % Not implemented
        L1B.product_type='L1B_LRRMC';
        display('L1B_LRRMC found')
    elseif strcmp('filename_L1B_LROS',filesBulkFields(i_file))
        [L1B]= reading_L1B_LROS (filesBulk, chd, cnf, cst); % Not implemented
        L1B.product_type='L1B_LROS';
        display('L1B_LROS found')
    elseif strcmp('filename_L1B_HR',filesBulkFields(i_file))
        [L1B]= reading_L1B_HR (filesBulk, chd, cnf, cst);
        L1B.product_type='L1B_HR';
        display('HR_1B found')
    elseif strcmp('filename_L1BFF_ML',filesBulkFields(i_file))
        [L1B]= reading_L1B_FF (filesBulk, chd, cnf, cst);
        L1B.product_type='L1B_FF';
        display('L1B_FF found')
    elseif isnan(L1B)
        display('L1B not found')
        continue;
    end
    
    % files in auxliar folder
    
    for i_file   =   1:length(aux)
        i_filename = aux{i_file};
        if(strfind(i_filename,'L1_pcnf'))
            filesBulk.L1_PCNF = [filesBulk.auxPath i_filename];
        elseif(strfind(i_filename,'GR_pcnf'))
            filesBulk.GR_PCNF = [filesBulk.auxPath i_filename];
        elseif(strfind(i_filename,'DEM'))
            filesBulk.DEM = [filesBulk.auxPath i_filename];
        elseif(strfind(i_filename,'geoid'))
            filesBulk.Geoid = [filesBulk.auxPath i_filename];
        elseif(strfind(i_filename,'mss'))
            filesBulk.MSS = [filesBulk.auxPath i_filename];
        elseif(strfind(i_filename, 'wet_tropospheric_correction'))
            filesBulk.wet_tropo_cor = [filesBulk.auxPath i_filename];
        elseif (strfind(i_filename, 'sea_ice_type'))
            filesBulk.sea_ice_type = [filesBulk.auxPath i_filename];
        elseif (strfind(i_filename, 'sea_ice_concentration'))
            filesBulk.sea_ice_concentration = [filesBulk.auxPath i_filename];
        end
    end
    % 2.2 Static auxiliary files
    if isfield(filesBulk, 'L1_PCNF')
        [L1_pcnf] = read_L1_PCNF(filesBulk);
    end
    if isfield(filesBulk, 'GR_PCNF')
        [GR_pcnf] = read_GR_pcnf(filesBulk);
    end
    if isfield(filesBulk, 'DEM') % canviar
        DEM = load('C:\Users\bballester\Documents\performance_tool_cristal\TLI1complexantenna_SARIn_OB_Ku\auxiliar\DEM_lat_lon_TLI1_SARIN_KU_CB.mat'); 
%         [DEM] = read_DEM(filesBulk);
    end
    if isfield(filesBulk, 'Geoid')
        [geoid] = read_geoid(filesBulk);
    end
    if isfield(filesBulk, 'MSS')
        [mss] = read_MSS(filesBulk);
    end
    
    % 2.3 Dynamic auxiliary files
    if isfield(filesBulk, 'wet_tropo_cor')
        [wet_tropo_cor] = read_wet_tropo_cor(filesBulk);
    end
    if isfield(filesBulk, 'sea_ice_type')
        [sea_ice_type] = read_sea_ice_type(filesBulk);
    end
    if isfield(filesBulk, 'sea_ice_concentration')
        [sea_ice_concentration] = read_sea_ice_concentration(filesBulk);
    end
    
    %% 3. Auxiliary data selection
%     
%     % 3.1 Geoid selection
%     [L2.geoid] = geoid_selection(L1B.lat, L1B.lon, geoid.h, geoid.lat, geoid.lon); % Not implemented yet
%     
%     % 3.2 MSS selection
%     [L2.mean_sea_surface] = mean_sea_surface_selection(L1B.lat, L1B.lon,  mss.h, mss.lat, mss.lon); % Not implemented yet
%     
%     % 3.3 DEM selection
%     %** cnf_gr
%     [length_swath, latitude_swath, longitude_swath, dem_height, dem_slope_h_cor, dem_height_swath] = ...
%         DEM_selection(L1B.lat, L1B.lon, L1B.interf_base_vector, DEM.h, DEM.lon, DEM.lat, DEM.step_lat,...
%         DEM.step_lon, cnf_L2.dem_interpolation_step);
%     
%     % 3.4 Sea ice concentration selection
%     [L2.sea_ice_concentration] = sea_ice_concentration_selection(L1B.time, L1B.lat, L1B.lon, cst.flat_coeff_cst, ...
%         cst.earth_radius_cst, sea_ice_concentration);
%     
%     % 3.5 Sea ice type selection
%     [L2.sea_ice_type] = sea_ice_type_selection(L1B.time, L1B.lat, L1B.lon, cst.flat_coeff_cst, cst.earth_radius_cst, ...
%         sea_ice_type.ice_type, sea_ice_type.xc, sea_ice_type.yc, sea_ice_type.lattrue, sea_ice_type.lonposy);
%     
%     % 3.6 Wet tropospheric correction selection
%     %** Input L2
%     [L2.rad_wet_tropo_cor, L2.rad_atm_cor_sig0] = rad_wet_tropo_cor_selection(L1B.time_tai, L2.time_tai, L2.rad_wet_tropo_cor, ...
%         L2.rad_wet_tropo_cor_qual, L2.rad_atm_cor_sig0_ku, L2.rad_atm_cor_sig0_ka, L2.rad_atm_cor_sig0_ku_qual, L2.rad_atm_cor_sig0_ka_qual);
    L2.rad_wet_tropo_cor = 0;
    L2.rad_atm_cor_sig0 = 0;
    L2.sea_ice_type = 0;
    L2.sea_ice_concentration = 0;
    L2.mean_sea_surface = 0;
   
    %% 4. Retrackers
%     
%     if strcmp('LR', L1B.product_type)
%         samples_aux = L1B.samples; %LR
%     else
%         samples_aux = L1B.samples_ov;%HR & FF
%     end
        
%% --------- 4.2 Waveform Characteristics
    [L2] = wfm_characteristics (L1B, L2, cnf_L2, chd); %**
        
%% ---------- 4.3 TFMRA Retracker
    if (cnf_L2.flag_tfmra_rtk_ku) %**
        [L2, validity_tfmra] = TFMRA_retracker (L1B, chd, cst, cnf_L2, L2);
    end
        
%% ---------- 4.3 OCOG Retracker
    if (cnf_L2.flag_ocog_rtk_ku) %**
        [L2, cog_ocog, width_ocog, validity_ocog] = OCOG_retracker(L1B, L2, chd, cst, cnf_L2);
    end
        
%% ---------- 4.3 TCOG Retracker %**
    if (cnf_L2.flag_tcog_rtk_ku)
        [L2, validity_tcog] = TCOG_retracker(L1B,L2, chd,cst, cnf_L2);
    end
    
%% ---------- 4.4 Physical Retracker

    if (cnf_L2.flag_physical_rtk_ku)
        L2.scaled_waveform = L2.scaled_waveform';
        [data, cst_p, chd_p, cnf_p] = read_adapt_cristal2retracker(L1B, L2, cnf_L2, cst, chd);
        filesBulk.filename_L1B=filesBulk.filename_L1B_HR;
        [L2.epoch, L2.range, L2.amplitude, L2.sig0, L2.swh, mqe, validity] =analytical_retracker_commonstruct...
            (data,L1B, L2.waveform_scale_factor, cnf_p, chd_p, cst_p, 'LUT_f0_file',filesBulk.LUT_f0_file,...
             'LUT_f1_file',filesBulk.LUT_f1_file,'path_Results',filesBulk.outputPath);
        L2.scaled_waveform = L2.scaled_waveform'; 
    end
%% ---------- 4.5 Swath Processing
    if strcmp('LR', L1B.product_type)==0
        if(cnf_L2.flag_swath_ku)
            L1B.scaled_waveforms = L2.scaled_waveform;
            L1B.range_oversampling_factor = cnf_L2.ZP;
            L1B.samples_ov = 1:1:size(L1B.time,1);
            L1B.phase_difference = L1B.phase_diff_meas_ku_l1b_echo;
            L1B.coherence = L1B.coherence_meas_ku_l1b_echo;
            L1B.range_ku_l1b_echo = L1B.tracker_range_calibrated;
            
            cnf_L2.zp_fact_range = cnf.zp_fact_range;
%             [SWATH,dem_height, dem_slope_h_cor, dem_height_swath, filesBulk]= ...
%             swath_processing (filesBulk,L1B, cnf_L2,chd,cst);
            DEM.height_m = double(DEM.height_m)';
            DEM.lat_deg = DEM.lat_deg';
%             DEM.lon_deg = DEM.lon_deg';
            SWATH=swath_processing(filesBulk,L1B,cnf_L2,chd,cst,DEM);

            if(cnf_L2.write_output)
                [filesBulk] = create_NetCDF_L2GR_Swath(filesBulk, N_bursts, cnf, chd, cst);
                [filesBulk] = write_NetCDF_L2GR_Swath(filesBulk, N_bursts, cnf, chd, cst);
            end
        end
    end
    
%% ---------- 5. Geophysical Parameters

%     [L2, range_snow, range_tfmra_snow]= geophysical_parameters(L1B, L2, cnf_L2, cst,  rad_wet_tropo_cor, dry_tropo_cor, ...
%         inv_bar_cor, dac, iono_cor, ocean_tide, long_period_tide, load_tide, solid_earth_tide, geocentric_polar_tide, ...
%         scaled_waveform, DEM.dem_slope_h_cor, SWATH);

    L2.load_tide = 0;
    cst.c_snow = 2.3e8; % I don't know if this is correct
    L1B.surface_classification_flag = 0; %needs to be read from the L1B (at the moment there is not this var)
    L1B.telemetry_type_flag = 0; %needs to be read from the L1B (at the moment there is not this var)
    
    [L2, range_snow, range_tfmra_snow]= geophysical_parameters(L1B, L2, cnf_L2, cst, chd,SWATH);
  
    %% Output Files
    if(cnf_L2.write_output)
        [filesBulk] = create_NetCDF_L2GR(filesBulk, N_bursts, cnf, chd, cst);
        [filesBulk] = write_NetCDF_L2GR(filesBulk, N_bursts, cnf, chd, cst);
    end
    
end
end
