%% ---------- GEOPHYSICAL CORRECTIONS AND FLAGS:ISD -----------------------------
ISD_corr_flags=(dec2bin(ncread(filename_L2_ISD,'flags_geo_corr_20_ku').'));

ISD_geo_corr.dry_trop=double(ncread(filename_L2_ISD,'dry_tropo_correction_20_ku')).';
ISD_geo_corr_flags.dry_trop=str2num(ISD_corr_flags(:,1));

ISD_geo_corr.wet_trop=double(ncread(filename_L2_ISD,'wet_tropo_correction_20_ku')).';
ISD_geo_corr_flags.wet_trop=str2num(ISD_corr_flags(:,2));

ISD_geo_corr.inv_baro=double(ncread(filename_L2_ISD,'inverse_baro_correction_20_ku')).';
ISD_geo_corr_flags.inv_baro=str2num(ISD_corr_flags(:,3));

ISD_geo_corr.dac=double(ncread(filename_L2_ISD,'Dynamic_atmospheric_correction_20_ku')).';
ISD_geo_corr_flags.dac=str2num(ISD_corr_flags(:,4));


ISD_geo_corr_flags.iono_corr=str2num(ISD_corr_flags(:,5));
if any(ISD_geo_corr_flags.iono_corr)
    ISD_geo_corr.iono_corr=double(ncread(filename_L2_ISD,'GIM_iono_correction_20_ku')).';    
else
    ISD_geo_corr.iono_corr=double(ncread(filename_L2_ISD,'model_iono_correction_20_ku')).';
    ISD_geo_corr_flags.iono_corr=str2num(ISD_corr_flags(:,6));
end

ISD_geo_corr.ocean_equilibrium_tide=double(ncread(filename_L2_ISD,'ocean_equilibrium_tide_20_ku')).';
ISD_geo_corr_flags.ocean_equilibrium_tide=str2num(ISD_corr_flags(:,7));

ISD_geo_corr.ocean_longperiod_tide=double(ncread(filename_L2_ISD,'long_period_tide_20_ku')).';
ISD_geo_corr_flags.ocean_longperiod_tide=str2num(ISD_corr_flags(:,8));

ISD_geo_corr.ocean_loading_tide=double(ncread(filename_L2_ISD,'ocean_loading_tide_20_ku')).';
ISD_geo_corr_flags.ocean_loading_tide=str2num(ISD_corr_flags(:,9));

ISD_geo_corr.solidearth_tide=double(ncread(filename_L2_ISD,'solid_earth_tide_20_ku')).';
ISD_geo_corr_flags.solidearth_tide=str2num(ISD_corr_flags(:,10));

ISD_geo_corr.geocentric_polar_tide=double(ncread(filename_L2_ISD,'geocentric_polar_tide_20_ku')).';
ISD_geo_corr_flags.geocentric_polar_tide=str2num(ISD_corr_flags(:,11));

%% ---------- GEOPHYSICAL CORRECTIONS AND FLAGS: ESA -----------------------------
ESA_geo_corr.dry_trop=reshape((ones(num_bursts_db,1))*(CS2.COR.dry_tropo.'),[1,ESA_num_surfaces]);
ESA_geo_corr_flags.dry_trop=reshape((CS2.MEA.Corr_Applic_flag_20Hz.corr_dry_tropo(:)),[1,ESA_num_surfaces]);

ESA_geo_corr.wet_trop=reshape((ones(num_bursts_db,1))*(CS2.COR.wet_tropo.'),[1,ESA_num_surfaces]);
ESA_geo_corr_flags.wet_trop=reshape((CS2.MEA.Corr_Applic_flag_20Hz.corr_wet_tropo(:)),[1,ESA_num_surfaces]);

ESA_geo_corr.inv_baro=reshape((ones(num_bursts_db,1))*(CS2.COR.inv_baro.'),[1,ESA_num_surfaces]);
ESA_geo_corr_flags.inv_baro=reshape((CS2.MEA.Corr_Applic_flag_20Hz.corr_inv_bar(:)),[1,ESA_num_surfaces]);

ESA_geo_corr.dac=reshape((ones(num_bursts_db,1))*(CS2.COR.dac.'),[1,ESA_num_surfaces]);
ESA_geo_corr_flags.dac=reshape((CS2.MEA.Corr_Applic_flag_20Hz.corr_hf_atmo(:)),[1,ESA_num_surfaces]);

ESA_geo_corr.iono_corr=reshape((ones(num_bursts_db,1))*(CS2.COR.iono.'),[1,ESA_num_surfaces]);
ESA_geo_corr_flags.iono_corr=reshape((CS2.MEA.Corr_Applic_flag_20Hz.corr_iono_gim(:)),[1,ESA_num_surfaces]); %assuming GIM is used
if ~any(ESA_geo_corr_flags.iono_corr)
    ESA_geo_corr_flags.iono_corr=reshape((CS2.MEA.Corr_Applic_flag_20Hz.corr_iono_model(:)),[1,ESA_num_surfaces]); %assuming GIM is used
end

ESA_geo_corr.ocean_equilibrium_tide=reshape((ones(num_bursts_db,1))*(CS2.COR.ocean_tide.'),[1,ESA_num_surfaces]);
ESA_geo_corr_flags.ocean_equilibrium_tide=reshape((CS2.MEA.Corr_Applic_flag_20Hz.corr_ocean_tide(:)),[1,ESA_num_surfaces]);

ESA_geo_corr.ocean_longperiod_tide=reshape((ones(num_bursts_db,1))*(CS2.COR.lpe_ocean.'),[1,ESA_num_surfaces]);
ESA_geo_corr_flags.ocean_longperiod_tide=reshape((CS2.MEA.Corr_Applic_flag_20Hz.corr_lpe_tide(:)),[1,ESA_num_surfaces]);

ESA_geo_corr.ocean_loading_tide=reshape((ones(num_bursts_db,1))*(CS2.COR.ocean_loading.'),[1,ESA_num_surfaces]);
ESA_geo_corr_flags.ocean_loading_tide=reshape((CS2.MEA.Corr_Applic_flag_20Hz.corr_ocean_loading(:)),[1,ESA_num_surfaces]);

ESA_geo_corr.solidearth_tide=reshape((ones(num_bursts_db,1))*(CS2.COR.solid_earth.'),[1,ESA_num_surfaces]);
ESA_geo_corr_flags.solidearth_tide=reshape((CS2.MEA.Corr_Applic_flag_20Hz.corr_se_tide(:)),[1,ESA_num_surfaces]);

ESA_geo_corr.geocentric_polar_tide=reshape((ones(num_bursts_db,1))*(CS2.COR.geoc_polar.'),[1,ESA_num_surfaces]);
ESA_geo_corr_flags.geocentric_polar_tide=reshape((CS2.MEA.Corr_Applic_flag_20Hz.corr_geo_polar_tide(:)),[1,ESA_num_surfaces]);


%% ------------------- aPPLY/COMPENSATION CORRECTIONS ---------------------
for i_surf=1:ISD_num_surfaces_filtered
    % dry tropo
    if ESA_geo_corr_flags.dry_trop(idx_int_ESA(i_surf))==1 && ISD_geo_corr_flags.dry_trop(idx_int_ISD(i_surf))==0
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))-ISD_geo_corr.dry_trop(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))-ISD_geo_corr.dry_trop(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))-ISD_geo_corr.dry_trop(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))-ISD_geo_corr.dry_trop(idx_int_ISD(i_surf));
        end
        
    elseif ESA_geo_corr_flags.dry_trop(idx_int_ESA(i_surf))==0 && ISD_geo_corr_flags.dry_trop(idx_int_ISD(i_surf))==1
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))+ISD_geo_corr.dry_trop(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))+ISD_geo_corr.dry_trop(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))+ISD_geo_corr.dry_trop(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))+ISD_geo_corr.dry_trop(idx_int_ISD(i_surf));
        end
    end
    % wet tropo
    if ESA_geo_corr_flags.wet_trop(idx_int_ESA(i_surf))==1 && ISD_geo_corr_flags.wet_trop(idx_int_ISD(i_surf))==0
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))-ISD_geo_corr.wet_trop(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))-ISD_geo_corr.wet_trop(idx_int_ISD(i_surf));
        end

        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))-ISD_geo_corr.wet_trop(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))-ISD_geo_corr.wet_trop(idx_int_ISD(i_surf));
        end        
    elseif ESA_geo_corr_flags.wet_trop(idx_int_ESA(i_surf))==0 && ISD_geo_corr_flags.wet_trop(idx_int_ISD(i_surf))==1
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))+ISD_geo_corr.wet_trop(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))+ISD_geo_corr.wet_trop(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))+ISD_geo_corr.wet_trop(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))+ISD_geo_corr.wet_trop(idx_int_ISD(i_surf));
        end        
    end
    % Inv barometric
    if ESA_geo_corr_flags.inv_baro(idx_int_ESA(i_surf))==1 && ISD_geo_corr_flags.inv_baro(idx_int_ISD(i_surf))==0
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))-ISD_geo_corr.inv_baro(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))-ISD_geo_corr.inv_baro(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))-ISD_geo_corr.inv_baro(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))-ISD_geo_corr.inv_baro(idx_int_ISD(i_surf));
        end         
    elseif ESA_geo_corr_flags.inv_baro(idx_int_ESA(i_surf))==0 && ISD_geo_corr_flags.inv_baro(idx_int_ISD(i_surf))==1
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))+ISD_geo_corr.inv_baro(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))+ISD_geo_corr.inv_baro(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))+ISD_geo_corr.inv_baro(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))+ISD_geo_corr.inv_baro(idx_int_ISD(i_surf));
        end
    end

    % DAC
    if ESA_geo_corr_flags.dac(idx_int_ESA(i_surf))==1 && ISD_geo_corr_flags.dac(idx_int_ISD(i_surf))==0
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))-ISD_geo_corr.dac(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))-ISD_geo_corr.dac(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))-ISD_geo_corr.dac(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))-ISD_geo_corr.dac(idx_int_ISD(i_surf));
        end
    elseif ESA_geo_corr_flags.dac(idx_int_ESA(i_surf))==0 && ISD_geo_corr_flags.dac(idx_int_ISD(i_surf))==1
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))+ISD_geo_corr.dac(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))+ISD_geo_corr.dac(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))+ISD_geo_corr.dac(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))+ISD_geo_corr.dac(idx_int_ISD(i_surf));
        end        
    end

    % Iono
    if ESA_geo_corr_flags.iono_corr(idx_int_ESA(i_surf))==1 && ISD_geo_corr_flags.iono_corr(idx_int_ISD(i_surf))==0
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))-ISD_geo_corr.iono_corr(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))-ISD_geo_corr.iono_corr(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))-ISD_geo_corr.iono_corr(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))-ISD_geo_corr.iono_corr(idx_int_ISD(i_surf));
        end        
    elseif ESA_geo_corr_flags.iono_corr(idx_int_ESA(i_surf))==0 && ISD_geo_corr_flags.iono_corr(idx_int_ISD(i_surf))==1
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))+ISD_geo_corr.iono_corr(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))+ISD_geo_corr.iono_corr(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))+ISD_geo_corr.iono_corr(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))+ISD_geo_corr.iono_corr(idx_int_ISD(i_surf));
        end
    end
    
    % ocean_equilibrium_tide
    if ESA_geo_corr_flags.ocean_equilibrium_tide(idx_int_ESA(i_surf))==1 && ISD_geo_corr_flags.ocean_equilibrium_tide(idx_int_ISD(i_surf))==0
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))-ISD_geo_corr.ocean_equilibrium_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))-ISD_geo_corr.ocean_equilibrium_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))-ISD_geo_corr.ocean_equilibrium_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))-ISD_geo_corr.ocean_equilibrium_tide(idx_int_ISD(i_surf));
        end
    elseif ESA_geo_corr_flags.ocean_equilibrium_tide(idx_int_ESA(i_surf))==0 && ISD_geo_corr_flags.ocean_equilibrium_tide(idx_int_ISD(i_surf))==1
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))+ISD_geo_corr.ocean_equilibrium_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))+ISD_geo_corr.ocean_equilibrium_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))+ISD_geo_corr.ocean_equilibrium_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))+ISD_geo_corr.ocean_equilibrium_tide(idx_int_ISD(i_surf));
        end        
    end    
    
    % ocean_longperiod_tide
    if ESA_geo_corr_flags.ocean_longperiod_tide(idx_int_ESA(i_surf))==1 && ISD_geo_corr_flags.ocean_longperiod_tide(idx_int_ISD(i_surf))==0
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))-ISD_geo_corr.ocean_longperiod_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))-ISD_geo_corr.ocean_longperiod_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))-ISD_geo_corr.ocean_longperiod_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))-ISD_geo_corr.ocean_longperiod_tide(idx_int_ISD(i_surf));
        end
    elseif ESA_geo_corr_flags.ocean_longperiod_tide(idx_int_ESA(i_surf))==0 && ISD_geo_corr_flags.ocean_longperiod_tide(idx_int_ISD(i_surf))==1
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))+ISD_geo_corr.ocean_longperiod_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))+ISD_geo_corr.ocean_longperiod_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))+ISD_geo_corr.ocean_longperiod_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))+ISD_geo_corr.ocean_longperiod_tide(idx_int_ISD(i_surf));
        end         
    end  
    
    % ocean_loading_tide
    if ESA_geo_corr_flags.ocean_loading_tide(idx_int_ESA(i_surf))==1 && ISD_geo_corr_flags.ocean_loading_tide(idx_int_ISD(i_surf))==0
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))-ISD_geo_corr.ocean_loading_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))-ISD_geo_corr.ocean_loading_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))-ISD_geo_corr.ocean_loading_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))-ISD_geo_corr.ocean_loading_tide(idx_int_ISD(i_surf));
        end        
    elseif ESA_geo_corr_flags.ocean_loading_tide(idx_int_ESA(i_surf))==0 && ISD_geo_corr_flags.ocean_loading_tide(idx_int_ISD(i_surf))==1
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))+ISD_geo_corr.ocean_loading_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))+ISD_geo_corr.ocean_loading_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))+ISD_geo_corr.ocean_loading_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))+ISD_geo_corr.ocean_loading_tide(idx_int_ISD(i_surf));
        end          
    end  

    % solidearth_tide
    if ESA_geo_corr_flags.solidearth_tide(idx_int_ESA(i_surf))==1 && ISD_geo_corr_flags.solidearth_tide(idx_int_ISD(i_surf))==0
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))-ISD_geo_corr.solidearth_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))-ISD_geo_corr.solidearth_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))-ISD_geo_corr.solidearth_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))-ISD_geo_corr.solidearth_tide(idx_int_ISD(i_surf));
        end         
    elseif ESA_geo_corr_flags.solidearth_tide(idx_int_ESA(i_surf))==0 && ISD_geo_corr_flags.solidearth_tide(idx_int_ISD(i_surf))==1
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))+ISD_geo_corr.solidearth_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))+ISD_geo_corr.solidearth_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))+ISD_geo_corr.solidearth_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))+ISD_geo_corr.solidearth_tide(idx_int_ISD(i_surf));
        end           
    end     

    
    % geocentric_polar_tide
    if ESA_geo_corr_flags.geocentric_polar_tide(idx_int_ESA(i_surf))==1 && ISD_geo_corr_flags.geocentric_polar_tide(idx_int_ISD(i_surf))==0
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))-ISD_geo_corr.geocentric_polar_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))-ISD_geo_corr.geocentric_polar_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))-ISD_geo_corr.geocentric_polar_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))-ISD_geo_corr.geocentric_polar_tide(idx_int_ISD(i_surf));
        end          
    elseif ESA_geo_corr_flags.geocentric_polar_tide(idx_int_ESA(i_surf))==0 && ISD_geo_corr_flags.geocentric_polar_tide(idx_int_ISD(i_surf))==1
        if ~isempty(idx_int_analytical_SWH_MSSfixed)
            ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_surf))+ISD_geo_corr.geocentric_polar_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_analytical_MSS_SWHfixed)
            ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))=ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_surf))+ISD_geo_corr.geocentric_polar_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_thres)
            ISD_L2_SSH_threshold(idx_int_ISD(i_surf))=ISD_L2_SSH_threshold(idx_int_ISD(i_surf))+ISD_geo_corr.geocentric_polar_tide(idx_int_ISD(i_surf));
        end
        if ~isempty(idx_int_ocog)
            ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))=ISD_L2_SSH_OCOG(idx_int_ISD(i_surf))+ISD_geo_corr.geocentric_polar_tide(idx_int_ISD(i_surf));
        end
    end
  
    
    
end