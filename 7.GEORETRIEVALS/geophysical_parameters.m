% function [L2, range_snow, range_tfmra_snow]= geophysical_parameters(L1B, L2, cnf_L2, cst,  rad_wet_tropo_cor, ...
%     dry_tropo_cor, inv_bar_cor, dac, iono_cor, ocean_tide, long_period_tide, load_tide, solid_earth_tide, ...
%     geocentric_polar_tide, scaled_waveform, dem_slope_h_cor, SWATH)

function [L2, range_snow, range_tfmra_snow]= geophysical_parameters(L1B, L2, cnf, cst, chd,SWATH)
%% 5.3 Snow depth
% Compute snow depth correction following Isobel publication
% Lawrence, I., Tsamados, M., Stroeve, J., Armitage, T., and Ridout, A.:
% Estimating snow depth over Arctic sea ice from calibrated dual-frequency % radar freeboards, The Cryosphere Discuss.,
% https://doi.org/10.5194/tc-2018-54, in review, 2018.
% Snow interface alignment shall be computed for Ku and Ka bands

% Increase Ku band range to align with the snow–ice interface
% if strcmp(chd.band , 'Ku')
%     L2.snow_interface_alignment = cnf_L2.snow_cor_slope*L2.peakiness+cnf_L2.snow_cor_intercept;
%     L2.range_tfmra_snow_ku=L2.range_tfmra_ku-L2.snow_interface_alignment_ku;
%     range_snow_ku=L2.range_ku-L2.snow_interface_alignment_ku;
% % Decrease Ka band range to align with the snow–air interface
% elseif strcmp(chd.band , 'Ka')
%     L2.snow_interface_alignment = cnf_L2.snow_cor_slope*L2.peakiness+cnf_L2.snow_cor_intercept;
%     L2.range_tfmra_snow_ka=L2.range_tfmra(cnf_L2.noise_first_sample_physical_rtk: cnf_L2.noise_last_sample_physical_rtk)-L2.snow_interface_alignment;
%     range_snow_ka=L2.range-L2.snow_interface_alignment;
% end
% % Computation of snow depth
% L2.snow_depth=(L2.range_tfmra_snow_ku-L2.range_tfmra_snow_ka)*cnf.c_snow/cst.c;
% % range_snow.range_snow_ku=range_snow_ku;
% % range_snow.range_snow_ka=range_snow_ka;
% % range_tfmra_snow.range_tfmra_snow_ku=range_tfmra_snow_ku;
% % range_tfmra_snow.range_tfmra_snow_ka=range_tfmra_snow_ka;
L2.snow_depth = 0; % CANVIAR
L2.range_tfmra_snow_ku = 0; % CANVIAR
range_snow = 0; % CANVIAR
range_tfmra_snow = 0; % CANVIAR
%% 5.4 Geophysical corrections
% Compute snow depth correction from snow depth in Ku band
% Set the snow depth correction to 0 in Ka band
if strcmp(chd.band, 'Ku')
    L2.snow_depth_correction = L2.snow_depth*(cst.c/cst.c_snow-1);
else % band = ‘Ka’
    L2.snow_depth_correction = 0;
end
%**
% Set to zero the corrections not computed by the SIRS prototype processor
ocean_tide = 0.;
long_period_tide = 0.;
load_tide = 0.;
solid_earth_tide=0.;
geocentric_polar_tide =0.;
dry_tropo_cor=0.;
iono_cor=0.;
inv_bar_cor=0.;
dynamic_atm=0.;
sea_state_bias=0.;
rad_wet_tropo_cor=0.;
dac=0.;
% Compute geophysical corrections
L2.seaice_geo_corrections = L2.snow_depth_correction + ocean_tide + long_period_tide + load_tide + solid_earth_tide+ geocentric_polar_tide + dry_tropo_cor + rad_wet_tropo_cor + iono_cor + inv_bar_cor;
L2.lead_geo_corrections = ocean_tide + long_period_tide + load_tide + solid_earth_tide+ geocentric_polar_tide + dry_tropo_cor + rad_wet_tropo_cor + iono_cor + inv_bar_cor;
L2.iceshelve_geo_corrections = L2.snow_depth_correction + ocean_tide + long_period_tide + load_tide + solid_earth_tide+ geocentric_polar_tide + dry_tropo_cor + rad_wet_tropo_cor + iono_cor + inv_bar_cor;
L2.landice_geo_corrections = L2.snow_depth_correction + load_tide+ solid_earth_tide+ geocentric_polar_tide + dry_tropo_cor + rad_wet_tropo_cor + inv_bar_cor;
L2.ocean_geo_corrections = ocean_tide + long_period_tide + load_tide+ solid_earth_tide+ geocentric_polar_tide + dry_tropo_cor + rad_wet_tropo_cor + iono_cor + dac + sea_state_bias;
L2.inland_geo_corrections = load_tide + solid_earth_tide+ geocentric_polar_tide + dry_tropo_cor + rad_wet_tropo_cor + iono_cor;


%% 5.5 Ocean Retrievals

L2.ssha=L1B.alt-(L2.range+L2.ocean_geo_corrections+L2.mean_sea_surface);
% L2.wind_speed_alt=model(L2.sig0,L2.swh); %dona error

%% 5.6 Sea Ice Retrievals

%     Waveform classification %FALTA
L2.waveform_classification_flag = 0;
% L2.waveform_classification_flag=classify_wvf(L1B.gaussian_fitting_std, L2.peakiness, L1B.surface_classification_flag,...
%     L2.sea_ice_concentration, cnf.up_peak_threshold_classif_ku, cnf.low_peak_threshold_classif_ku, ...
%     cnf.stack_std_threshold_classif_ku);

if cnf.flag_sea_ice_rtk % use range from the physical retracker
%     Sea surface removal
    ssha_lead=L1B.altitude-(L2.range+L2.lead_geo_corrections+L2.mean_sea_surface);
    ssha_lead_interpolated=interpolate_ssha_lead(ssha_lead,L1B.time_tai,L2.waveform_classification_flag);
    L2.ssh_lead_interpolated=ssha_lead_interpolated+L2.mean_sea_surface;
%     Sea ice freeboard
    L2.sea_ice_freeboard=L1B.altitude-(range_snow_ku+L2.seaice_geo_corrections)-L2.ssh_lead_interpolated;

else % use range from the TFMRA retracker
%     Sea surface removal
    ssha_lead=L1B.alt-(L2.range_tfmra+L2.lead_geo_corrections+L2.mean_sea_surface);
    % falta fer funcio interpolate_ssha_lead
    L2.ssha_lead_interpolated = 0;
%     ssha_lead_interpolated=interpolate_ssha_lead(ssha_lead,L1B.time,L2.waveform_classification_flag);
    L2.ssh_lead_interpolated=L2.ssha_lead_interpolated+L2.mean_sea_surface;
%     Sea ice freeboard
    L2.sea_ice_freeboard=L1B.alt-(L2.range_tfmra_snow_ku+L2.seaice_geo_corrections)-L2.ssh_lead_interpolated;
end
% falta fer funcio interpolate_ssha_lead
L2.ssha_lead_interpolated = 0;
% ssha_lead_interpolated=interpolate_ssha_lead(ssha_lead,L1B.time,L2.waveform_classification_flag);
L2.ssh_lead_interpolated=L2.ssha_lead_interpolated+L2.mean_sea_surface;
%     Sea ice thickness
L2.sea_ice_thickness=(L2.sea_ice_freeboard * cnf.ice_density + L2.snow_depth * cnf.snow_density) / ...
    (cnf.sea_water_density - cnf.sea_ice_density);

%% 5.7 Iceberg Detection
if cnf.iceberg_detection
    disp('------------Entering iceberg detection algorithm------------')
    [icebergs_properties_L1B]=iceberg_detection(filesBulk,L1B,cnf,cnf_L2,chd,cst);
else
    disp('------------Skipping iceberg detection algorithm------------')
end

%% 5.8 Land-ice retrievals

%land ice scenario, for land ice surfaces, or for the ice shelve scenario, for ice shelve surfaces. 
%The surface classification mask can be used to determine the appropriate geophysical correction to apply.

if (L1B.surface_classification_flag == 1) % land ice %**
    geo_corr=L2.landice_geo_corrections;
else % ice shelve
    geo_corr=L2.iceshelve_geo_corrections;
end
switch L1B.telemetry_type_flag
    case 'SARIn'
        L2.land_ice_elevation = (L1B.altitude - (L2.range_tfmra+geo_corr).*cos(SWATH.angle_of_arrival) + ...
            cst.earth_radius*(1-cos(L1B.altitude/cst.earth_radius.*SWATH.angle_of_arrival)).* ...
            cos(L1B.altitude/cst.earth_radius.*SWATH.angle_of_arrival));
    otherwise
        dem_slope_h_cor = 0; % Set to 0 in SAR mode within CRISTAL SIR
        L2.land_ice_elevation = (L1B.alt- (L2.range_tfmra+geo_corr+dem_slope_h_cor));
end
% latitude_surf & longitude_surf
switch cnf.mode
    case 'SAR'
        L2.latitude_surf = L1B.lat;
        L2.longitude_surf = L1B.lon;
    otherwise
        L2.latitude_surf = SWATH.lat;
        L2.longitude_surf = SWATH.lon;
end

%% 5.9 Inland waters retrievals
L2.water_level_height = L1B.alt-(L2.range+L2.inland_geo_corrections);

end
