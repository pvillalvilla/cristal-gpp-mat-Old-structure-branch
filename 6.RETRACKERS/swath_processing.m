% function [SWATH]=swath_processing(filesBulk,L1B,cnf,chd,cst)
function [SWATH,filesBulk]=swath_processing(filesBulk,L1B,cnf,chd,cst)

%MOVE TO THE CONFIG FILE
DEM_resolution = 0.00005; % TO DO: to be moved to cnf 0.00005== 5 meters resolution 0.0005==50 meters resolution , aprox 20 min
%DEM_resolution = 0.00002; % TO DO: to be moved to cnf 0.00005== 5 meters resolution 0.0005==50 meters resolution
coherence_threshold = 0.92;
coherence_threshold = 0.77;
coherence_threshold = 0.15;
coherence_threshold = 0.179;
coherence_threshold = 0.50;
power_threshold = -170;
%surf_distance       = 3; % distance in samples to consider a differencet surface from two consecutives high coherence periods
surf_distance       = 30; % distance in samples to consider a differencet surface from two consecutives high coherence periods
%surf_threshold      = 30; % minimum number of samples for a surface to be considered valid
surf_threshold      = 1; % minimum number of samples for a surface to be considered valid
peak_threshold      = 0; % The scaling factor has to be used, not yet included.
N_wraps             = 3; % wraps for the DEM of the plot comparison
coh_margin_beggining    = 20; % samples to avoid at the beggining of the coherence
%coh_margin_beggining    = 0; % samples to avoid at the beggining of the coherence
%coh_margin_beggining    = 242; % samples to avoid at the beggining of the coherence
 coh_margin_end          = 300; % samples to avoid a the end of the coherence
% coh_margin_end          = 738; % samples to avoid a the end of the coherence
%coh_margin_end          = 768; % samples to avoid a the end of the coherence


maskFlag            = 0;
plotting_swath_individuals  = 1;
plotting_swath_summary      = 1;

[~,outputPath_size] = size(filesBulk.outputPath);
plots_folder = [filesBulk.outputPath filesBulk.filename_L1B(outputPath_size+1:end-3) '/'];

if(~exist (plots_folder, 'dir'))
    mkdir(plots_folder);
end
if(plotting_swath_individuals) 
    font_size = 12;
    WaterColor = [143 226 255]/255;
    mida = get(0,'ScreenSize');
    %mida(3:4)=[1920,1080];
%     mida(3:4)=[1280,720];
    h=figure('Position', mida); subplot(3,4,3);
    
end
    
SWATH   = struct('lat',[],'lon',[],'elev',[],'N_selected',[],...
        'individual_dem_diff',[],'MAD',[],'average_dem_diff',[],'std_dem_diff',[],...
        'wvf_number',[],'sample',[]);

filesBulk.DEM_file= dir([filesBulk.auxPath 'DEM*']);
                
DEM = load([filesBulk.auxPath filesBulk.DEM_file.name]);
%build new DEM interpolated using 0.0005 -> 55 meters

%DEM struct reorganisation, in TLI1v6 Aresys changed names
% DEM.lon_deg=DEM.lon;
% DEM.lat_deg=DEM.lat;
% DEM.lat_step_deg=mean(diff(DEM.lat));
% DEM.lon_step_deg=mean(diff(DEM.lon));
% DEM.height_m=DEM.height;

 DEM.height_m=DEM.height_m.'; % In TLI 17 Aresys swapped lat/lon dimensions

[lon_grid,lat_grid]=meshgrid(min(DEM.lon_deg):DEM.lon_step_deg:max(DEM.lon_deg),min(DEM.lat_deg):DEM.lat_step_deg:max(DEM.lat_deg));
iDEMgood=~(isnan(DEM.height_m));
DEM.height_m = double(DEM.height_m).';
%F = scatteredInterpolant(lon_grid(iDEMgood.'),lat_grid(iDEMgood.'),double(DEM.height_m(iDEMgood.')));

elev_limit_plot= max(max(DEM.height_m));
elev_limit_min_plot= min(min(DEM.height_m));

% DEM scatter plot
%figure;scatter3(reshape(lon_grid,[1 5736*4281]),reshape(lat_grid, [1 5736*4281]), reshape(DEM.height_m, [1 5736*4281]))

DEM_reduced.lon_deg=DEM.lon_deg(1000:3000);
DEM_reduced.lat_deg=DEM.lat_deg(1500:4000);
DEM_reduced.height=DEM.height_m(1500:4000,1000:3000);
% DEM reduced scatter plot
[lon_grid_red,lat_grid_red]=meshgrid(min(DEM_reduced.lon_deg):DEM.lon_step_deg:max(DEM_reduced.lon_deg),min(DEM_reduced.lat_deg):DEM.lat_step_deg:max(DEM_reduced.lat_deg));
%figure;scatter3(reshape(lon_grid_red,[1 2501*2001]),reshape(lat_grid_red, [1 2501*2001]), reshape(DEM_reduced.height, [1 2501*2001]))


[lon_grid,lat_grid]=meshgrid(min(DEM.lon_deg):DEM.lon_step_deg:max(DEM.lon_deg),min(DEM.lat_deg):DEM.lat_step_deg:max(DEM.lat_deg));
iDEMgood=~(isnan(DEM.height_m));

DEM.height_m = double(DEM.height_m).'; % In 


create_dem_hd=0;
if create_dem_hd
    tic
    [lon_hd,lat_hd]= meshgrid(min(DEM.lon_deg):DEM_resolution:max(DEM.lon_deg),min(DEM.lat_deg):DEM_resolution:max(DEM.lat_deg));
    F = scatteredInterpolant(lon_grid(iDEMgood.'),lat_grid(iDEMgood.'),double(DEM.height_m(iDEMgood.')));
    
    DEM_hd = F(lon_hd,lat_hd);
    toc
end

B = abs(chd.y_cog_ant_2- chd.y_cog_ant);
cst.flat_coeff=0.003352810664747;
Eccentricity = sqrt(2*cst.flat_coeff-cst.flat_coeff^2);
geometric_factor = 1.113; % Geometric factor due to earth curvature

SURF                = [];
elev_average_error  = [];
ambiguity_wrap      = [];
elevation_std       = [];
surf_index          = [];
surf_samples_length = [];
init_sample         = [];

ellipsoid           =   referenceEllipsoid ('wgs84');
records_L1b=length(L1B.alt);

idx_counter=0;
%for index_L1b=1: records_L1b
%for index_L1b=50: records_L1b    
%for index_L1b=125:127  %JPLZ: do it for specific records TLI.1
%for index_L1b=124:128  %JPLZ: do it for specific records TLI.2
for index_L1b=125:127  %JPLZ: do it for specific records TLI.2
    idx_counter=idx_counter+1; % JPLZ: a idx for when we don't process the whole segment
    
    if(index_L1b>1)
        v =[L1B.lat(index_L1b),L1B.lon(index_L1b)]-[L1B.lat(index_L1b-1),L1B.lon(index_L1b-1)];
        [arclen_v,az_v] = distance('gc',L1B.lat(index_L1b-1),L1B.lon(index_L1b-1),L1B.lat(index_L1b),L1B.lon(index_L1b),[cst.semi_major_axis Eccentricity]);
    else
        v =[L1B.lat(index_L1b+1),L1B.lon(index_L1b+1)]-[L1B.lat(index_L1b),L1B.lon(index_L1b)];
        [arclen_v,az_v] = distance('gc',L1B.lat(index_L1b),L1B.lon(index_L1b),L1B.lat(index_L1b+1),L1B.lon(index_L1b+1),[cst.semi_major_axis Eccentricity]);
       
    end
   
    %% 2.COMPUTE X_AXIS
    top_elev        = L1B.alt(index_L1b)+chd.N_samples_sar/2*chd.T0_nom/cnf.zp_fact_range*cst.c/2;
    bottom_elev     = L1B.alt(index_L1b)-(chd.N_samples_sar/2)*chd.T0_nom/cnf.zp_fact_range*cst.c/2;
    far_Range       = L1B.range_ku_l1b_echo(index_L1b)+chd.N_samples_sar/2*chd.T0_nom/cnf.zp_fact_range*cst.c/2;
    short_Range     = L1B.range_ku_l1b_echo(index_L1b)-(chd.N_samples_sar/2)*chd.T0_nom/cnf.zp_fact_range*cst.c/2;
    
%     x_axisL1 = (top_elev-chd.T0_nom/cnf.zp_fact_range*cst.c/2:-chd.T0_nom/cnf.zp_fact_range*cst.c/2:bottom_elev);
    Range = (short_Range:chd.T0_nom/cnf.zp_fact_range*cst.c/2:far_Range-chd.T0_nom/cnf.zp_fact_range*cst.c/2);
  
    roll=0;
%    roll=0.001; %JPLZ: testing
    original_AoA = (chd.wv_length * (L1B.phase_diff_meas_ku_l1b_echo(index_L1b,:))./(2*pi*B)-roll)/pi*180;   %[deg] 
    AoA    = (smooth(chd.wv_length * unwrap(L1B.phase_diff_meas_ku_l1b_echo(index_L1b,:))./(2*pi*B)-roll,'moving')/pi*180).';   %[deg]
%      AoA    = (smooth(chd.wv_length * (L1B.phase_diff_meas_ku_l1b_echo(index_L1b,:))./(2*pi*B)-roll)/pi*180).';   %[deg]  
  
 % Antenna pattern phase difference compensation (TO BE IMPLEMENTED)
 filesBulk.antennapattern_file= dir([filesBulk.auxPath 'antenna_pattern*']);        
 antenna_pattern = load([filesBulk.auxPath filesBulk.antennapattern_file.name]);
  
    %%CHECK
    if AoA(350)>(chd.wv_length *pi./(2*pi*B))/pi*180
        AoA_shift = AoA(350)-original_AoA(350);
        AoA = AoA-AoA_shift;
        
    elseif AoA(350)<-(chd.wv_length * pi./(2*pi*B))/pi*180
        AoA_shift = AoA(350)-original_AoA(350);
        AoA = AoA-AoA_shift;
    end
        
        
        
    Xdist  =  sin(AoA /180*pi) .* Range;
    jump_AoA     = (chd.wv_length * 2*pi./(2*pi*B))/pi*180;   %[deg] 360deg to AoA
    
%     AoA_smooth = smooth(AoA);
%     th_AoA     = (chd.wv_length * 3/2*pi./(2*pi*B))/pi*180;   %[deg] theoretical Angle of Arrival, to check for jumps
%     
%     Slant_Range    = sqrt(Xdist.^2 + Range.^2);
    
%     Slant_Elev     = -Slant_Range+L1B.alt(index_L1b);
%     L2.Elev(index_L1b,:)   = cos(AoA.*pi./180).*Slant_Elev ; % need to be corrected for delta_Ha and delta_Hp: SEE TOM ARMITAGE PAPER https://www.researchgate.net/publication/260623282_Using_the_Interferometric_Capabilities_of_the_ESA_CryoSat-2_Mission_to_Improve_the_Accuracy_of_Sea_Ice_Freeboard_Retrievals
    
    if(index_L1b>1)
        v =[L1B.lat(index_L1b),L1B.lon(index_L1b)]-[L1B.lat(index_L1b-1),L1B.lon(index_L1b-1)];
        [arclen_v,az_v] = distance('gc',L1B.lat(index_L1b-1),L1B.lon(index_L1b-1),L1B.lat(index_L1b),L1B.lon(index_L1b),[cst.semi_major_axis Eccentricity]);
        [arclen_v,az_v] = distance('gc',L1B.lat(index_L1b-1),L1B.lon(index_L1b-1),L1B.lat(index_L1b),L1B.lon(index_L1b),wgs84Ellipsoid);       
    else
        v =[L1B.lat(index_L1b+1),L1B.lon(index_L1b+1)]-[L1B.lat(index_L1b),L1B.lon(index_L1b)];
        [arclen_v,az_v] = distance('gc',L1B.lat(index_L1b),L1B.lon(index_L1b),L1B.lat(index_L1b+1),L1B.lon(index_L1b+1),[cst.semi_major_axis Eccentricity]);
        [arclen_v,az_v] = distance('gc',L1B.lat(index_L1b),L1B.lon(index_L1b),L1B.lat(index_L1b+1),L1B.lon(index_L1b+1),wgs84Ellipsoid);
    end
    %az_v along track direction foreward
    %az_d always vector perpendicular to de track to the right side of the along track direction
    
    az_d=az_v-90;
    %az_d=az_v+90;
    %             if(v(1)>0) % ascending pass
    %                 az_d=az_v-90;
    %             elseif (v(1)<0) % descending pass
    %                 az_d=az_v-90;
    %             end
    
    [lats,lons] = reckon(L1B.lat(index_L1b),L1B.lon(index_L1b),Xdist,az_d,[cst.semi_major_axis Eccentricity]);
    AoA_points  = -jump_AoA*(N_wraps):jump_AoA:jump_AoA*(N_wraps);% [-jump_AoA jump_AoA] = case N_wraps = 1;
    %                 AoA_points = AoA_points(AoA_points~=0);
    phase_wraps     = -2*pi*N_wraps:2*pi:2*pi*N_wraps;
    phase_wraps     = phase_wraps(phase_wraps~=0);
    Ambiguities      = phase_wraps * chd.wv_length * L1B.range_ku_l1b_echo(index_L1b)./(2*pi*B);
    lats_nadir = zeros(1,length(Range))+ L1B.alt(index_L1b);
    lons_nadir = zeros(1,length(Range))+ L1B.lon(index_L1b);
    for i_wrap=1: length(AoA_points)
        Xdist_wrap(i_wrap,:)  =  sin((AoA+AoA_points(i_wrap)) /180*pi) .* Range;
        nadir_Range_wrap(i_wrap,:) = cos((AoA+AoA_points(i_wrap)) /180*pi).* Range;
        nadir_Elev_wrap(i_wrap,:) = (L1B.alt(index_L1b) - nadir_Range_wrap(i_wrap,:));
        %          R_Earth =sqrt((cst.semi_major_axis^4*cos(L1B.alt(index_L1b).*pi/180).^2+cst.semi_minor_axis^4*sin(L1B.alt(index_L1b)*pi/180).^2)./((cst.semi_major_axis*cos(L1B.alt(index_L1b).*pi/180)).^2+(cst.semi_minor_axis*sin(L1B.alt(index_L1b).*pi/180)).^2));
        %         alpha        =   1 + (nf_p.h)./(R_Earth);
        AoA_corr = antenna_pattern_compensation(AoA+AoA_points(i_wrap),antenna_pattern);
        % AoA_corr   =   AoA; %[deg]
        % AoA_corr   =   AoA+AoA_points(i_wrap); %[deg]
        
        [SURF_wrap(index_L1b).lat(i_wrap,:),SURF_wrap(index_L1b).lon(i_wrap,:),Elev_wrap(i_wrap,:)]   = aer2geodetic(az_d,AoA_corr-90,Range,L1B.lat(index_L1b),L1B.lon(index_L1b),L1B.alt(index_L1b),wgs84Ellipsoid);
        
%         %TLI.1 Glacier vertexes
%         TLI1_vertex1=[-0.48882730 -0.384333 +37]; %[lat lon] [rad]
%         TLI1_vertex2=[-0.48878721 -0.38389 -33]; %[lat lon] [rad]
%         TLI1_vertex3=[-0.48860515 -0.38391 -33]; %[lat lon] [rad]
%         TLI1_vertex4=[-0.48864414 -0.384356 +37]; %[lat lon] [rad]
%         %glacier azimuth size
%         TLI1_glacier_az_size=1160; %[m]
        
%         %TLI.1 Glacier vertexes comples antenna simulation
%         TLI1_vertex1=[1.4494266 -0.424305 +710]; %[lat lon alt] [rad, m]
%         TLI1_vertex2=[1.4517192 -0.363142 -651]; %[lat lon alt] [rad, m]
%         TLI1_vertex3=[1.4518122 -0.363372 -651]; %[lat lon alt] [rad, m]
%         TLI1_vertex4=[1.4495186 -0.424571 +710]; %[lat lon alt] [rad, m]
%         %glacier azimuth size
%         TLI1_glacier_az_size=620; %[m]
%
%         %TLI.1 v6 Glacier vertexes comples antenna simulation
%         TLI1_vertex1=[-0.53464893 -1.62505 -655]; %[lat lon alt] [rad, m]
%         TLI1_vertex2=[-0.53541752 -1.61605 +714]; %[lat lon alt] [rad, m]
%         TLI1_vertex3=[-0.53527015 -1.61603 +714]; %[lat lon alt] [rad, m]
%         TLI1_vertex4=[-0.53449828 -1.62503 -655]; %[lat lon alt] [rad, m]
%         %glacier azimuth size
%         TLI1_glacier_az_size=620; %[m]

%         %TLI.3-16 Glacier vertexes
%         TLI3_vertex1=[1.450612 -0.394397]; %[lat lon] [rad]
%         TLI3_vertex2=[1.450730 -0.391173]; %[lat lon] [rad]
%         TLI3_vertex3=[1.450872 -0.391534]; %[lat lon] [rad]
%         TLI3_vertex4=[1.450753 -0.394761]; %[lat lon] [rad]
%         %glacier azimuth size
%         TLI3_glacier_az_size=962; %[m]
        
%         %TLI.3-16 Glacier vertexes complex antenna simulation
%         TLI3_vertex1=[1.4485975 -0.443381 +2291.8]; %[lat lon alt] [rad, m]
%         TLI3_vertex2=[1.4523562 -0.343078 -2179.6]; %[lat lon alt] [rad, m]
%         TLI3_vertex3=[1.4524503 -0.343280 -2180.4]; %[lat lon alt] [rad, m]
%         TLI3_vertex4=[1.4486890 -0.443658 +2291.2]; %[lat lon alt] [rad, m]
%         %glacier azimuth size
%         TLI3_glacier_az_size=620; %[m]
        
        
%         %TLI.1 
%         inside= inpolygon(SURF_wrap(index_L1b).lat(i_wrap,:) ,SURF_wrap(index_L1b).lon(i_wrap,:) ,[TLI1_vertex1(1) TLI1_vertex2(1) TLI1_vertex3(1) TLI1_vertex4(1) TLI1_vertex1(1)]*180/pi,[TLI1_vertex1(2) TLI1_vertex2(2) TLI1_vertex3(2) TLI1_vertex4(2) TLI1_vertex1(2)]*180/pi); %TLI1 glacier vertexs       
%         for i=1:chd.N_samples_sar
%             if inside(i)==1
% %                dist_middle_glacier_upper(i)=distance('gc', ((TLI1_vertex1(1)+TLI1_vertex4(1))/2)*180/pi, ((TLI1_vertex1(2)+TLI1_vertex4(2))/2)*180/pi, SURF_wrap(index_L1b).lat(i_wrap,i), SURF_wrap(index_L1b).lon(i_wrap,i),[cst.semi_major_axis Eccentricity]);
% %                dist_upper(i)=min(nonzeros(dist_middle_glacier_upper));
% %                dist_upper_modifier(i)=dist_upper(i)*90/(TLI1_glacier_az_size/2); %TLI3  glacier 0.932km azimuth size
% %                dist_upper_factor(i)=cos(dist_upper_modifier(i)*pi/180);
%                 %dist_middle_glacier_lower(i)=distance('gc', 83.124774213359345, -22.422903847673062, SURF_wrap(index_L1b).lat(i_wrap,i), SURF_wrap(index_L1b).lon(i_wrap,i),[cst.semi_major_axis Eccentricity]);
%                  dist_beg_glacier_upper(i)=distance('gc', TLI1_vertex1(1)*180/pi, TLI1_vertex1(2)*180/pi, SURF_wrap(index_L1b).lat(i_wrap,i), SURF_wrap(index_L1b).lon(i_wrap,i),[cst.semi_major_axis Eccentricity]);
%                  dist_end_glacier_upper(i)=distance('gc', TLI1_vertex4(1)*180/pi, TLI1_vertex4(2)*180/pi, SURF_wrap(index_L1b).lat(i_wrap,i), SURF_wrap(index_L1b).lon(i_wrap,i),[cst.semi_major_axis Eccentricity]);                 
%                  if min(nonzeros(dist_beg_glacier_upper))>277 && min(nonzeros(dist_end_glacier_upper))>277 %277=Ku along track resolution
%                      Elev_wrap(i_wrap,i)=Elev_wrap(i_wrap,i)-(2.05-0.05); %surf ref height-snow_height
%                  end
%             end
%         end
%        clear inside        
%         %TLI.3
%         inside= inpolygon(SURF_wrap(index_L1b).lat(i_wrap,:) ,SURF_wrap(index_L1b).lon(i_wrap,:) ,[TLI3_vertex1(1) TLI3_vertex2(1) TLI3_vertex3(1) TLI3_vertex4(1) TLI1_vertex1(1)]*180/pi,[TLI3_vertex1(2) TLI3_vertex2(2) TLI3_vertex3(2) TLI3_vertex4(2) TLI3_vertex1(2)]*180/pi); %TLI3 glacier vertexs
%         for i=1:chd.N_samples_sar
%             if inside(i)==1
%                 dist_middle_glacier_upper(i)=distance('gc', 83.117984663487050, -22.607711384492511, SURF_wrap(index_L1b).lat(i_wrap,i), SURF_wrap(index_L1b).lon(i_wrap,i),[cst.semi_major_axis Eccentricity]);
%                 dist_upper(i)=min(nonzeros(dist_middle_glacier_upper));
%                 dist_upper_modifier(i)=dist_upper(i)*90/(TLI3_glacier_az_size/2); %TLI3  glacier 0.932km azimuth size
%                 dist_upper_factor(i)=cos(dist_upper_modifier(i)*pi/180);
%                 %dist_middle_glacier_lower(i)=distance('gc', 83.124774213359345, -22.422903847673062, SURF_wrap(index_L1b).lat(i_wrap,i), SURF_wrap(index_L1b).lon(i_wrap,i),[cst.semi_major_axis Eccentricity]);
%                 Elev_wrap(i_wrap,i)=Elev_wrap(i_wrap,i)-dist_upper_factor(i)*(2.05-0.05); %surf ref height-snow_height
%             end
%         end
%        clear inside

%         dist_upper=min(nonzeros(dist_middle_glacier_upper));
%         dist_mod_lower=min(nonzeros(dist_middle_glacier_lower));
%         dist_upper_modifier=dist_upper*90/466;
%         dist_upper_factor=cos(dist_upper_modifier*pi/180);

    
%         xyz_nadir =  lla2ecef([lats_nadir.' lons_nadir.' nadir_Elev_wrap(i_wrap,:).'],'wgs84');
%         ellipsoid_axis_nadir = sqrt(xyz_nadir(:,1).^2+xyz_nadir(:,2).^2+xyz_nadir(:,3).^2);
%         xyz =  lla2ecef([lats_nadir.' lons_nadir.' nadir_Elev_wrap(i_wrap,:).'],'wgs84');
%         ellipsoid_axis_nadir = sqrt(xyz(:,1).^2+xyz(:,2).^2+xyz(:,3).^2);
%         Elev_wrap2(i_wrap,:) = sqrt(Xdist_wrap(i_wrap,:).^2+(nadir_Elev_wrap(i_wrap,:)+R_Earth).^2)-R_Earth;
%         Slant_Range_wrap    = sqrt(Xdist_wrap(i_wrap,:).^2 + Range.^2);
%         Slant_Elev_wrap(i_wrap,:)     = -Slant_Range_wrap+ L1B.range_ku_l1b_echo(index_L1b)+L1B.alt(index_L1b);
%         Elev_wrap2(i_wrap,:)   = cos((AoA+AoA_points(i_wrap)).*pi./180).*Slant_Elev_wrap(i_wrap,:)-L1B.range_ku_l1b_echo(index_L1b); % need to be corrected for delta_Ha and delta_Hp: SEE TOM ARMITAGE PAPER https://www.researchgate.net/publication/260623282_Using_the_Interferometric_Capabilities_of_the_ESA_CryoSat-2_Mission_to_Improve_the_Accuracy_of_Sea_Ice_Freeboard_Retrievals
%         
%         Elev_wrap2(i_wrap,:) = sqrt(Xdist_wrap(i_wrap,:).^2+(nadir_Elev_wrap(i_wrap,:)+ ellipsoid_axis_nadir.').^2)-ellipsoid_axis_nadir.';

        
    end
    
    Coh_pos = find((L1B.coherence_meas_ku_l1b_echo(index_L1b,1:end-coh_margin_end)>coherence_threshold));
    Coh_pos=Coh_pos(Coh_pos>coh_margin_beggining);
    surf_trans_pos    = find(diff(Coh_pos)>surf_distance);

    SURF(index_L1b).N_surfs 	= length(surf_trans_pos)+1;
    clear  surf_samples point_color N_selected mean_DEM_error wrap_selected 
    SURF(index_L1b).peakness_surf = zeros (1,SURF(index_L1b).N_surfs);
    
    i_surf_index=0;
    SURF(index_L1b).N_surfs_OK=0;

    if(plotting_swath_individuals)
        subplot(4,3,[2 3 5 6]); hold on
        dem_2D=imagesc(DEM.lon_deg,DEM.lat_deg,DEM.height_m);
%         [distance_axis,angle_dist]=distance('gc',L1B.lat(index_L1b),0,lat_dem_vector,lon_dem_vector-L1B.lon(index_L1b),[cst.semi_major_axis Eccentricity]);
%         distance_axis=distance_axis.*sign(cosd(angle_dist));
%         dem_2D=imagesc(-distance_axis,DEM.lat_deg,DEM.height_m);
        demcmap(DEM.height_m);
        colorbar;
        uistack(dem_2D,'bottom');
        set(gca,'YDir','normal');
    end

    
    for i_surf=1: SURF(index_L1b).N_surfs
        i_surf_index = i_surf_index+1;
        SURF(index_L1b).POCA=0;
        SURF(index_L1b).Wrap_flag(i_surf_index)=0;
        if(i_surf==1)
            SURF(index_L1b).surf_ini_sample(i_surf_index)    = Coh_pos(1);
        else
            SURF(index_L1b).surf_ini_sample(i_surf_index)    = Coh_pos(surf_trans_pos(i_surf-1)+1);
            %SURF(index_L1b).surf_ini_sample(i_surf_index)    = Coh_pos(1); %JPLZ: change it, dunno why we take this length
        end
        
        if(i_surf==SURF(index_L1b).N_surfs)
            SURF(index_L1b).surf_end_sample(i_surf_index)    = Coh_pos(end);
        else
            SURF(index_L1b).surf_end_sample(i_surf_index)    = Coh_pos(surf_trans_pos(i_surf));
            %SURF(index_L1b).surf_end_sample(i_surf_index)    = Coh_pos(end); %JPLZ: change it, dunno why we take this length
        end
        
        surf_samples        = SURF(index_L1b).surf_ini_sample(i_surf_index):SURF(index_L1b).surf_end_sample(i_surf_index);
        if (length(surf_samples)>surf_threshold)
            SURF(index_L1b).N_surfs_OK = SURF(index_L1b).N_surfs_OK+1;
            clear DEM_diff MAD
            for i_wrap=1: length(AoA_points)
                clear point_color
%                 % locations of the section: SURF_wrap(index_L1b).lat(i_wrap,surf_samples),SURF_wrap(index_L1b).lon(i_wrap,surf_samples)
%                 % elevatiosn plot_Elev_wrap(i_wrap,surf_samples)
                 lat_dem_vector = SURF_wrap(index_L1b).lat(i_wrap,surf_samples);
                 [lat_diff,lat_dem_index]= min(abs((DEM.lat_deg-lat_dem_vector(1))));
                 lon_dem_vector = SURF_wrap(index_L1b).lon(i_wrap,surf_samples);
%                 %[lon_diff,lon_dem_index]= min(abs((DEM.lon_deg-lon_dem_vector(1)))); %JPLZ
%                 DEM_slope=polyfit(DEM.lon_deg,DEM.height_m(lat_dem_index,:),1);
%                 %DEM_slope2=polyfit(DEM.lat_deg,DEM.height_m(:,lon_dem_index),1);%JPLZ
%                  h_DEM = polyval(DEM_slope,lon_dem_vector);
%                  h_DEM2 = interp1(DEM.lon_deg,DEM.height_m(lat_dem_index,:), lon_dem_vector.','linear').'; %interpolation of the swath longitudes and the DEM
if create_dem_hd==1
    [lat_diff_hd,lat_dem_index_hd]= min(abs((lat_hd(:,1)-lat_dem_vector(1))));    
end

%build new DEM interpolated using 0.0005 -> 55 meters
%  [lon_grid,lat_grid]=meshgrid(min(DEM.lon_deg):DEM.lon_step_deg:max(DEM.lon_deg),min(DEM.lat_deg):DEM.lat_step_deg:max(DEM.lat_deg));
% iDEMgood=~(isnan(DEM.height_m));
% DEM.height_m = double(DEM.height_m).';
% F = scatteredInterpolant(lon_grid(iDEMgood.'),lat_grid(iDEMgood.'),double(DEM.height_m(iDEMgood.')));
% create_dem_hd=0;
% if create_dem_hd
%     [lon_hd,lat_hd]= meshgrid(min(DEM.lon_deg):DEM_resolution:max(DEM.lon_deg),min(DEM.lat_deg):DEM_resolution:max(DEM.lat_deg));
%     F = scatteredInterpolant(lon_grid(iDEMgood.'),lat_grid(iDEMgood.'),double(DEM.height_m(iDEMgood.')));
%     
%     tic
%     DEM_hd = F(lon_hd,lat_hd);
%     toc
% end

% create_dem_hd=0;
% if create_dem_hd && index_L1b==108 && i_wrap==5
%     for i=1:length(DEM.lon_deg)
%         min_lon(i)=abs(min(lon_dem_vector)- DEM.lon_deg(i));
%         max_lon(i)=abs(max(lon_dem_vector)- DEM.lon_deg(i));
%     end
%     for i=1:length(DEM.lat_deg)
%         min_lat(i)=abs(min(lat_dem_vector)- DEM.lat_deg(i));
%         max_lat(i)=abs(max(lat_dem_vector)- DEM.lat_deg(i));
%     end
%     [val_min_lon,idx_min_lon]=min(min_lon);
%     [val_max_lon,idx_max_lon]=min(max_lon);
%     [val_min_lat,idx_min_lat]=min(min_lat);
%     [val_max_lat,idx_max_lat]=min(max_lat);
%     
%     %figure;mesh(DEM.lon_deg(idx_min_lon-100:idx_max_lon+100),DEM.lat_deg(idx_min_lat-100:idx_max_lat+100) ,DEM.height_m(idx_min_lat-100:idx_max_lat+100,idx_min_lon-100:idx_max_lon+100))
%     %scatter3(SURF_wrap(index_L1b).lat(i_wrap,surf_samples),SURF_wrap(index_L1b).lon(i_wrap,surf_samples),Elev_wrap(i_wrap,surf_samples))
%     
%     lon_grid=DEM.lon_deg(idx_min_lon-100):DEM.lon_step_deg:DEM.lon_deg(idx_max_lon+100);
%     lat_grid=DEM.lat_deg(idx_min_lat-100):DEM.lat_step_deg:DEM.lat_deg(idx_max_lat+100);
%     
%     [lon_grid,lat_grid]=meshgrid(DEM.lon_deg(idx_min_lon-100):DEM.lon_step_deg:DEM.lon_deg(idx_max_lon+100),DEM.lat_deg(idx_min_lat-100):DEM.lat_step_deg:DEM.lat_deg(idx_max_lat+100));
%     iDEMgood=~(isnan(DEM.height_m(idx_min_lat-100:idx_max_lat+100,idx_min_lon-100:idx_max_lon+100)));
%     %DEM.height_m = double(DEM.height_m).';
%     [lon_hd,lat_hd]= meshgrid(DEM.lon_deg(idx_min_lon-100):DEM_resolution:DEM.lon_deg(idx_max_lon+100),DEM.lat_deg(idx_min_lat-100):DEM_resolution:DEM.lat_deg(idx_max_lat+100));
%     %F = scatteredInterpolant(lat_grid(iDEMgood.'),lon_grid(iDEMgood.'),double(DEM.height_m(idx_min_lat-100:idx_max_lat+100,idx_min_lon-100:idx_max_lon+100)));
%     F = scatteredInterpolant(lat_grid,lon_grid,DEM.height_m(idx_min_lat-100:idx_max_lat+100,idx_min_lon-100:idx_max_lon+100));
% 
%     
%     tic
%     DEM_hd = F(lat_hd,lon_hd);
%     toc
% end
% if create_dem_hd
%     h_DEM_hd = F(SURF_wrap(index_L1b).lat(i_wrap,surf_samples),SURF_wrap(index_L1b).lon(i_wrap,surf_samples));
% end

%                 DEM_diff(i_wrap,:) = Elev_wrap(i_wrap,surf_samples)-h_DEM;
%figure;mesh(DEM.lon_deg(:),DEM.lat_deg(2000:2500),DEM.height_m(2000:2500,:));hold on;
%figure;mesh(DEM.lon_deg(:),DEM.lat_deg(2000:2500),DEM.height_m(2000:2500,:));hold on;

%plot3(L1B.lon(90:120),L1B.lat(90:120),L1B.alt(90:120)-L1B.range_ku_l1b_echo(90:120),'-k');
% if index_L1b >= 105 && i_wrap==4
%                  scatter3(lon_dem_vector(:), lat_dem_vector(:),h_DEM(:),'w')
%                  scatter3(SURF_wrap(index_L1b).lon(i_wrap,surf_samples), SURF_wrap(index_L1b).lat(i_wrap,surf_samples), Elev_wrap(i_wrap,surf_samples), 'c')
 %                 scatter3(SURF_wrap(index_L1b).lat(i_wrap,surf_samples), SURF_wrap(index_L1b).lon(i_wrap,surf_samples), Elev_wrap(i_wrap,surf_samples))
% end
                for ipoint=1:length(lat_dem_vector)
                    [~,lat_dem_index(ipoint)]= min(abs((lat_grid(:,1)-lat_dem_vector(ipoint))));
                    %[~,lat_dem_index]= min(abs((DEM.lat_deg(:,1)-lat_dem_vector(ipoint)))); %JPLZ
                    [~,lon_dem_index(ipoint)]= min(abs((lon_grid(1,:)-lon_dem_vector(ipoint))));
                    %[~,lon_dem_index]= min(abs((DEM.lon_deg(1,:)-lon_dem_vector(ipoint)))); %JPLZ
                    h_DEM(ipoint)= (DEM.height_m(lat_dem_index(ipoint),lon_dem_index(ipoint)));
                    
                    
                    %DEM_hd(ipoint)= (DEM.height_m(lat_dem_index,lon_dem_index));
                    %h_DEM(ipoint)= (DEM.height_m(lat_dem_index,lon_dem_index)); %JPLZ
                    %DEM_diff(i_wrap,ipoint) = Elev_wrap(i_wrap,surf_samples(ipoint))-h_DEM(ipoint);
                    DEM_diff(i_wrap,ipoint) = Elev_wrap(i_wrap,surf_samples(ipoint))-h_DEM(ipoint)-2.05; % JPLZ: glacier is 2.05m above DEM

                    

                    
                    if create_dem_hd==1
                        [~,lat_dem_index_hd(ipoint)]= min(abs((lat_hd(:,1)-lat_dem_vector(ipoint))));
                        [~,lon_dem_index_hd(ipoint)]= min(abs((lon_hd(1,:)-lon_dem_vector(ipoint))));
                        h_DEM_hd(ipoint)= (DEM_hd(lat_dem_index_hd(ipoint),lon_dem_index_hd(ipoint)));
                        DEM_diff_hd(i_wrap,ipoint) = Elev_wrap(i_wrap,surf_samples(ipoint))-h_DEM_hd(ipoint);
                    end
                    
                    %set colors
                    if(abs(DEM_diff(i_wrap,ipoint))<10)
                        point_color(ipoint,:) = [0;1;0].';
                    elseif(abs(DEM_diff(i_wrap,ipoint))<50)
                        point_color(ipoint,:) = [1;1;0].';
                    elseif(abs(DEM_diff(i_wrap,ipoint))>50)
                        point_color(ipoint,:) = [1;0;0].';
                    else
                        point_color(ipoint,:) = [0;0;0].'; % Case DEM is NaN
                    end
                end
                
%                 lat_dem_index=unique(lat_dem_index);
%                 [lon_dem_index,ia,ic]=unique(lon_dem_index);
%                 h_DEM= (DEM.height_m(lat_dem_index.',lon_dem_index));
                %h_DEM_interp1 = spline(lon_grid(1,lon_dem_index),h_DEM(ia),lon_dem_vector);
                 %DEM_diff(i_wrap,:) = Elev_wrap(i_wrap,surf_samples)-h_DEM_interp1;
%                 MAD(i_wrap) = std(DEM_diff(i_wrap,:)-std(DEM_diff(i_wrap,:)));
                MAD(i_wrap) = median(abs(DEM_diff(i_wrap,:)-median(DEM_diff(i_wrap,:))));
   %index_L1b==104 && i_wrap==8 
%        if(plotting_swath_individuals && i_wrap==length(AoA_points))
%         subplot(4,3,[2 3 5 6]); hold on
%         %dem_2D=imagesc(DEM.lon_deg,DEM.lat_deg,DEM.height_m);
%         [distance_axis,angle_dist]=distance('gc',L1B.lat(index_L1b),0,lat_dem_vector,lon_dem_vector-L1B.lon(index_L1b),[cst.semi_major_axis Eccentricity]);
%         distance_axis=distance_axis.*sign(cosd(angle_dist));
%         dem_2D=imagesc(distance_axis,DEM.lat_deg,DEM.height_m);
%         demcmap(DEM.height_m);
%         colorbar;
%         uistack(dem_2D,'bottom');
%         set(gca,'YDir','normal');
%     end
   
                if(plotting_swath_individuals)
                    subplot(4,3,1);plot(1:chd.N_samples_sar,L1B.scaled_waveforms(index_L1b,:),'k','LineWidth',2);hold on;
                    subplot(4,3,4);plot(1:chd.N_samples_sar,L1B.coherence_meas_ku_l1b_echo(index_L1b,:),'k','LineWidth',2);hold on;
                    %                 plot(1:N_samples_sar_chd,plot_Coh_smooth,'b','LineWidth',1);
                    subplot(4,3,7); plot(1:chd.N_samples_sar,AoA,'og','LineWidth',2);hold on;
                                    plot(1:chd.N_samples_sar,original_AoA,'.k','LineWidth',2);
                    
                    
                    subplot(4,3,[2 3 5 6]); 
                    p_wrap(i_wrap,i_surf_index)=plot3(SURF_wrap(index_L1b).lon(i_wrap,surf_samples),SURF_wrap(index_L1b).lat(i_wrap,surf_samples),Elev_wrap(i_wrap,surf_samples),'o');
                    %p_wrap(i_wrap,i_surf_index)=plot3(Xdist_wrap(i_wrap,surf_samples),SURF_wrap(index_L1b).lat(i_wrap,surf_samples),Elev_wrap(i_wrap,surf_samples),'o');
                    subplot(4,3,8:9);
                    plot(SURF_wrap(index_L1b).lon(i_wrap,surf_samples),Elev_wrap(i_wrap,surf_samples),'o'); hold on;
                    %plot(Xdist_wrap(i_wrap,surf_samples),Elev_wrap(i_wrap,surf_samples),'o'); hold on;
                    subplot(4,3,11:12);
                    plot(SURF_wrap(index_L1b).lon(i_wrap,surf_samples),DEM_diff(i_wrap,:),'o'); hold on;
                    %plot(Xdist_wrap(i_wrap,surf_samples),DEM_diff(i_wrap,:),'o'); hold on;
                    
                    
                end

            end
            [mean_DEM_error(i_surf_index),wrap_selected(i_surf_index)]=min(abs(mean(DEM_diff.')));
            [mean_MAD(i_surf_index),wrap_selected(i_surf_index)]=min(MAD);
            N_selected(i_surf_index)= (-wrap_selected(i_surf_index)+N_wraps+1);
            elev_average_error=[elev_average_error; mean_DEM_error(i_surf_index)];
            surf_index=[surf_index; i_surf_index];
            ambiguity_wrap=[ambiguity_wrap; N_selected(i_surf_index)];
            elevation_std=[elevation_std; std(DEM_diff(wrap_selected(i_surf_index),:))];
            surf_samples_length=[surf_samples_length;length(surf_samples)];
            init_sample=[init_sample;surf_samples(1)];

            SWATH.lat               =[SWATH.lat SURF_wrap(index_L1b).lat(wrap_selected(i_surf_index),surf_samples)];
            SWATH.lon               =[SWATH.lon SURF_wrap(index_L1b).lon(wrap_selected(i_surf_index),surf_samples)];
            SWATH.elev              =[SWATH.elev Elev_wrap(wrap_selected(i_surf_index),surf_samples)];
            SWATH.N_selected        =[SWATH.N_selected N_selected(i_surf_index)+zeros(1,length(surf_samples))];
            SWATH.individual_dem_diff   =[SWATH.individual_dem_diff DEM_diff(wrap_selected(i_surf_index),:)];
            SWATH.MAD               =[SWATH.MAD mean_MAD(i_surf_index)];
            SWATH.average_dem_diff  = [SWATH.average_dem_diff mean_DEM_error(i_surf_index)];
            SWATH.std_dem_diff      = [SWATH.std_dem_diff std(DEM_diff(wrap_selected(i_surf_index),:))];
            SWATH.wvf_number        = [SWATH.wvf_number index_L1b+zeros(1,length(surf_samples))];
            SWATH.sample            = [SWATH.sample surf_samples];
        end
        if(plotting_swath_individuals)
            
            subplot(4,3,1);
            plot(surf_samples,L1B.scaled_waveforms(index_L1b,surf_samples),'Color',[(i_surf_index-1)/SURF(index_L1b).N_surfs 1-(i_surf_index-1)/SURF(index_L1b).N_surfs 0],'LineWidth',1);
            plot((SURF(index_L1b).surf_ini_sample(i_surf_index)),L1B.scaled_waveforms(index_L1b,SURF(index_L1b).surf_ini_sample(i_surf_index)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
            plot((SURF(index_L1b).surf_end_sample(i_surf_index)),L1B.scaled_waveforms(index_L1b,SURF(index_L1b).surf_ini_sample(i_surf_index)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
            
            subplot(4,3,4);
            plot(surf_samples,L1B.coherence_meas_ku_l1b_echo(index_L1b,surf_samples),'Color',[(i_surf_index-1)/SURF(index_L1b).N_surfs 1-(i_surf_index-1)/SURF(index_L1b).N_surfs 0],'LineWidth',1);
            plot(SURF(index_L1b).surf_ini_sample(i_surf_index),L1B.coherence_meas_ku_l1b_echo(index_L1b,SURF(index_L1b).surf_ini_sample(i_surf_index)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
            plot(SURF(index_L1b).surf_end_sample(i_surf_index),L1B.coherence_meas_ku_l1b_echo(index_L1b,SURF(index_L1b).surf_end_sample(i_surf_index)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');

            subplot(4,3,7); 
            plot(surf_samples,AoA(surf_samples),'Color',[(i_surf_index-1)/SURF(index_L1b).N_surfs 1-(i_surf_index-1)/SURF(index_L1b).N_surfs 0],'LineWidth',1);
            plot(SURF(index_L1b).surf_ini_sample(i_surf_index),AoA(SURF(index_L1b).surf_ini_sample(i_surf_index)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
            plot(SURF(index_L1b).surf_end_sample(i_surf_index),AoA(SURF(index_L1b).surf_end_sample(i_surf_index)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
        end
    end
    
    
%% PLOTING

if(plotting_swath_individuals) 
    subplot(4,3,1); hold off;
    figlabels('Elevation [m]','Power [dB]','',['#',num2str(index_L1b,'%03d'),' L1B Waveform, #',num2str(SURF(index_L1b).N_surfs_OK,'%d'),' Reflections retracked'] ,font_size); set(gca,'XLim',[1 chd.N_samples_sar],'FontSize',font_size);
    subplot(4,3,4); hold off;
    figlabels('samples','','',['#',num2str(index_L1b,'%03d'),' Coherence'],font_size); set(gca,'XLim',[1 chd.N_samples_sar],'FontSize',font_size);set(gca,'YLim',[0 1],'FontSize',font_size);
    subplot(4,3,7); hold off;
    figlabels('samples','AoA [deg]','',['#',num2str(index_L1b,'%03d'), ' Angle of Arrival - Roll'],font_size); set(gca,'XLim',[1 chd.N_samples_sar],'FontSize',font_size); set(gca,'YLim',[-0.6*3 0.6*3],'FontSize',font_size);
    legend('local unwrapped AoA', 'original AoA','Location','southeast');
    AoA_points          = -jump_AoA*(N_wraps+1/2):jump_AoA/2:jump_AoA*(N_wraps+1/2);% [-jump_AoA/2 0 jump_AoA/2]
    top_swath_Xdist     = ((AoA_points)/180*pi-roll) * (L1B.range_ku_l1b_echo(index_L1b)/cst.c-chd.N_samples_sar/2*chd.T0_nom/cnf.zp_fact_range)*cst.c/2;
    bottom_swath_Xdist  = ((AoA_points)/180*pi-roll) * (L1B.range_ku_l1b_echo(index_L1b)/cst.c+(chd.N_samples_sar/2-1)*chd.T0_nom/cnf.zp_fact_range)*cst.c/2;
    top_Slant_Range     = sqrt(top_swath_Xdist.^2 + (L1B.range_ku_l1b_echo(index_L1b)/cst.c*cst.c/2)^2);
    bottom_Slant_Range  = sqrt(bottom_swath_Xdist.^2 + (L1B.range_ku_l1b_echo(index_L1b)/2)^2);
    top_Slant_Elev      = top_Slant_Range- L1B.range_ku_l1b_echo(index_L1b)/2+top_elev;
    top_Elev            = cos((AoA_points).*pi./180 -roll).*top_Slant_Elev; % need to be corrected for delta_Ha and delta_Hp: SEE TOM ARMITAGE PAPER https://www.researchgate.net/publication/260623282_Using_the_Interferometric_Capabilities_of_the_ESA_CryoSat-2_Mission_to_Improve_the_Accuracy_of_Sea_Ice_Freeboard_Retrievals
    bottom_Slant_Elev   = bottom_Slant_Range- L1B.range_ku_l1b_echo(index_L1b)/2+bottom_elev;
    bottom_Elev         = cos((AoA_points).*pi./180-roll).*bottom_Slant_Elev; % need to be corrected for delta_Ha and delta_Hp: SEE TOM ARMITAGE PAPER https://www.researchgate.net/publication/260623282_Using_the_Interferometric_Capabilities_of_the_ESA_CryoSat-2_Mission_to_Improve_the_Accuracy_of_Sea_Ice_Freeboard_Retrievals
    [lat_top_elev,lon_top_elev]         = reckon(L1B.lat(index_L1b),L1B.lon(index_L1b),top_swath_Xdist,az_d,[cst.semi_major_axis Eccentricity]);
    [lat_bottom_elev,lon_bottom_elev]   = reckon(L1B.lat(index_L1b),L1B.lon(index_L1b),bottom_swath_Xdist,az_d,[cst.semi_major_axis Eccentricity]);
    lon_bottom_elev_pairs               = lon_bottom_elev(1:2:end);
    lon_top_elev_pairs                  = lon_top_elev(1:2:end);
    lat_bottom_elev_pairs               = lat_bottom_elev(1:2:end);
    lat_top_elev_pairs                  = lat_top_elev(1:2:end);
    bottom_Elev_pairs                   = bottom_Elev(1:2:end);
    top_Elev_pairs                      = top_Elev(1:2:end);
    [lat_elev_DEM_north,lon_elev_DEM_north] = reckon(lat_bottom_elev_pairs,lon_bottom_elev_pairs,100e3,az_v,[cst.semi_major_axis Eccentricity]);
    [lat_elev_DEM_south,lon_elev_DEM_south] = reckon(lat_bottom_elev_pairs,lon_bottom_elev_pairs,100e3,az_v+180,[cst.semi_major_axis Eccentricity]);
    
    for i_wrap=1: length(top_Elev_pairs)
        lon_elev_wrap(i_wrap,:)     = [lon_bottom_elev_pairs(i_wrap) lon_top_elev_pairs(i_wrap)];
        lat_elev_wrap(i_wrap,:)     = [lat_bottom_elev_pairs(i_wrap) lat_top_elev_pairs(i_wrap)];
        lon_elev_wrap_DEM(i_wrap,:) = [lon_elev_DEM_north(i_wrap) lon_elev_DEM_south(i_wrap)];
        lat_elev_wrap_DEM(i_wrap,:) = [lat_elev_DEM_north(i_wrap) lat_elev_DEM_south(i_wrap)];
        elev_wrap(i_wrap,:)         = [elev_limit_min_plot elev_limit_plot];
        str(i_wrap,:)               = ['N =' num2str(i_wrap-length(top_Elev_pairs)/2,'%+i ')];
        elev_wrapdiff(i_wrap,:)     = [-100 300];
    end
    
    str(length(top_Elev_pairs)/2,:)='NADIR';
    [lat_dem1,lon_dem1] = reckon(L1B.lat(index_L1b),L1B.lon(index_L1b),(8000*N_wraps*2),az_d,[cst.semi_major_axis Eccentricity]);  % N_wraps warps
    [lat_dem2,lon_dem2] = reckon(L1B.lat(index_L1b),L1B.lon(index_L1b),(8000*N_wraps*2),az_d+180,[cst.semi_major_axis Eccentricity]);  %N_wraps 2 wraps
    lat_res =(max(lat_dem2,lat_dem1)-(min(lat_dem2,lat_dem1)))/1000;
    lon_res =(max(max(DEM.lon_deg),min(DEM.lon_deg))-(min(max(DEM.lon_deg),min(DEM.lon_deg))))/1000;
    
    if(lon_dem1<lon_dem2)
        lon_dem_vector = min(DEM.lon_deg):lon_res:max(DEM.lon_deg);
    elseif(lon_dem1>lon_dem2)
        lon_dem_vector = min(DEM.lon_deg):lon_res:max(DEM.lon_deg);
    end
    lat_dem_vector = lat_dem1:lat_res:lat_dem2;
    if (lat_res==0)
        lat_dem_vector = zeros(1,length(lon_dem_vector))+lat_dem1;
    end
    for ipoint=1:length(lat_dem_vector)
        [~,lat_dem_index]= min(abs((DEM.lat_deg-lat_dem_vector(ipoint))));
        [~,lon_dem_index]= min(abs((DEM.lon_deg-lon_dem_vector(ipoint))));
        h_DEM_strip(ipoint)= double(DEM.height_m(lat_dem_index,lon_dem_index));
        
    end
    
    %                 h_DEM_strip(~h_DEM_strip)= NaN;
    subplot(4,3,[2 3 5 6]); 
     %plot(lon_dem_vector,lat_dem_vector,':k','LineWidth',6);
     %plot(lon_elev_wrap_DEM',lat_elev_wrap_DEM','--k','LineWidth',1,'MarkerFaceColor','k','MarkerEdgeColor','k');
    long_shift=+0.025;
    if(v(1)<0) %descending
        long_shift=-0.025;
    end
    if (lat_res==0)
        long_shift=0;   
    end
    text(lon_top_elev(end-1:-2:2)+long_shift,lat_top_elev(end-1:-2:2).'+0.5,cellstr(str(1:end-1,:)),'Color','black','FontSize',12,'HorizontalAlignment','center');
    colorbar('Location','northoutside');hold off;
    figlabels('Longitude [degrees]','Latitude [degrees]','Elevation [m]',['# ' num2str(index_L1b,'%04d') ' record: Top view DEM + Swaths computed for the different Nwraps'] ,font_size);
    set(gca,'XLim',[min(DEM.lon_deg) max(DEM.lon_deg)],'FontSize',font_size);
    %figlabels('X distance [m]','Latitude [degrees]','Elevation [m]',['# ' num2str(index_L1b,'%04d') ' record: Top view DEM + Swaths computed for the different Nwraps'] ,font_size);
    %set(gca,'XLim',[-5e4 +5e4],'FontSize',font_size);
    set(gca,'YLim',[min(DEM.lat_deg) max(DEM.lat_deg)],'FontSize',font_size);
    set(gca,'ZLim',[elev_limit_min_plot elev_limit_plot],'FontSize',font_size);
    
    subplot(4,3,8:9);
    %                 x_range=length(lon_dem_vector(~isnan(h_DEM_strip)));
    im_area=imagesc(lon_dem_vector(~isnan(h_DEM_strip)),elev_limit_min_plot:1:elev_limit_plot,(elev_limit_min_plot:1:elev_limit_plot).');hold on
    set(gca,'YDir','normal');
    dem_line=area(lon_dem_vector(~isnan(h_DEM_strip)),h_DEM_strip(~isnan(h_DEM_strip)),'LineWidth',1,'BaseValue',elev_limit_plot,'FaceColor',WaterColor);

    %testing plots with distance across track
    %[distance_axis,angle_dist]=distance('gc',L1B.lat(index_L1b),0,lat_dem_vector,lon_dem_vector-L1B.lon(index_L1b),[cst.semi_major_axis Eccentricity]);
    %distance_axis=distance_axis.*sign(cosd(angle_dist));
    %im_area=imagesc(distance_axis(~isnan(h_DEM_strip)),elev_limit_min_plot:1:elev_limit_plot,(elev_limit_min_plot:1:elev_limit_plot).');hold on
    %distance_axis_mod=distance_axis;
    %small_vec=-6953:105:2072;
    %distance_axis_mod(548:633)=small_vec;
    %dem_line=area(distance_axis_mod(~isnan(h_DEM_strip)),-h_DEM_strip(~isnan(h_DEM_strip)),'LineWidth',1,'BaseValue',elev_limit_plot,'FaceColor',WaterColor);
  
   
    
    %     demcmap(DEM_hd);
    %                 colormap(colors_new);
    uistack(dem_line,'bottom');
    uistack(im_area,'bottom');
    plot(lon_elev_wrap',elev_wrap','--k','LineWidth',1,'MarkerFaceColor','k','MarkerEdgeColor','k');
    text(lon_top_elev(end-1:-2:2),elev_wrapdiff(1:end-1,1).'.*0+500,cellstr(str(1:end-1,:)),'Color','black','FontSize',12,'HorizontalAlignment','center')
    set(gca,'box','on');
    set(gca,'Color',WaterColor);
    set(gca,'XLim',[min(DEM.lon_deg) max(DEM.lon_deg)],'FontSize',font_size);
    %set(gca,'XLim',[min(Xdist_wrap(i_wrap,surf_samples)) max(Xdist_wrap(i_wrap,surf_samples))],'FontSize',font_size);
    set(gca,'YLim',[elev_limit_min_plot elev_limit_plot],'FontSize',font_size);
    subplot(4,3,11:12);
    %                 plot(lon_bottom_elev,bottom_Elev,'--k','LineWidth',1)
    %                 plot(lon_top_elev,top_Elev,'--k','LineWidth',1)
    %                 plot(lon_bottom_elev(1:2:end),bottom_Elev(1:2:end),'--k','LineWidth',1,'MarkerFaceColor','r','MarkerEdgeColor','k');
    %                 plot(lon_top_elev(1:2:end),top_Elev(1:2:end),'--k','LineWidth',1,'MarkerFaceColor','r','MarkerEdgeColor','k');
    plot(lon_elev_wrap',elev_wrapdiff','--k','LineWidth',1,'MarkerFaceColor','k','MarkerEdgeColor','k');
    %                 set(gca,'Color',WaterColor);
    
    text(lon_top_elev(end-1:-2:2),elev_wrapdiff(1:end-1,1).'.*0-70,cellstr(str(1:end-1,:)),'Color','black','FontSize',12,'HorizontalAlignment','center')
    set(gca,'box','on');
    
     
%     if(~isempty(SWATH.lon))
%         scatter3(SWATH.lon,SWATH.lat,SWATH.elev,30,SWATH.elev,'filled','MarkerEdgeColor','k','LineWidth',0.5);
%     end
    
    %                 view(-10,75);
    
    
    
    subplot(4,3,8:9);hold off;
    figlabels('Longitude [degrees]','Elevation [m]','',['# ' num2str(index_L1b,'%04d') ' record: Elevation profile DEM + Swaths computed for the different Nwraps #',num2str(index_L1b,'%04d')] ,font_size);
    set(gca,'XLim',[min(DEM.lon_deg) max(DEM.lon_deg)],'FontSize',font_size); %JPLZ: commented to do some plots
    %figlabels('X dist [m]','Elevation [m]','',['# ' num2str(index_L1b,'%04d') ' record: Elevation profile DEM + Swaths computed for the different Nwraps #',num2str(index_L1b,'%04d')] ,font_size);
    %set(gca,'XLim',[-5e4 +5e4],'FontSize',font_size);
    set(gca,'YLim',[elev_limit_min_plot elev_limit_plot],'FontSize',font_size);
    
    
    subplot(4,3,11:12);hold off;
    figlabels('Longitude [degrees]','DEM diff [m]','',['# ' num2str(index_L1b,'%04d') ' record: DEM difference for the different Nwraps #',num2str(index_L1b,'%04d') ] ,font_size);
    set(gca,'XLim',[min(DEM.lon_deg) max(DEM.lon_deg)],'FontSize',font_size); %JPLZ: commented to do some plots
    %figlabels('X dist [m]','DEM diff [m]','',['# ' num2str(index_L1b,'%04d') ' record: DEM difference for the different Nwraps #',num2str(index_L1b,'%04d') ] ,font_size);
    %set(gca,'XLim',[-5e4 +5e4],'FontSize',font_size);
    set(gca,'YLim',[-100 300],'FontSize',font_size);
    
    subp=subplot(4,3,10);
    set(subp,'Visible','off');
    dim = get(subp,'position');
    columnname =   {'init','Ns', 'Best Wrap', 'Mean Elev. diff','MAD'};
    columnformat = {'numeric','numeric' , 'numeric', 'numeric', 'numeric'};
    dat={};
    %dat=horzcat((init_sample(index_L1b)),(surf_samples_length(index_L1b)),(ambiguity_wrap(index_L1b)),(elev_average_error(index_L1b)),(mean_MAD(i_surf_index)));
    %dat=horzcat((init_sample(index_L1b)),(surf_samples_length(index_L1b)),(ambiguity_wrap(index_L1b)),(elev_average_error(index_L1b)));
    dat = horzcat((init_sample(idx_counter)),(surf_samples_length(idx_counter)),(ambiguity_wrap(idx_counter)),(elev_average_error(idx_counter)),(SWATH.MAD(idx_counter)));
    t = uitable('Units','normalized','Position',dim, 'Data', (dat),...
        'ColumnName', columnname,...
        'ColumnFormat', columnformat,'ColumnWidth', {30,30,70,90,90},...
        'RowName',[],'FontName','Arial','Units','normalized');
    t_pos =[0.1300 0.0149 0.1809 0.2525];
    set(t,'Position',t_pos);
    axis off;
    
    
    
    img = getframe(h);
%     imwrite(img.cdata, [plots_folder num2str(index_L1b,'%04d'),'_DEM_and_swath.png']);
    set(h,'PaperPositionMode','auto');
    saveas(h,[plots_folder num2str(index_L1b,'%04d'),'_DEM_and_swath.png'])
%     fig = gcf;
% h.PaperPositionMode = 'auto';
% print('ScreenSizeFigure','-dpng','-r0')
%     print('ScreenSizeFigure','-dpng','-r0');
    
end




end

if(plotting_swath_summary)
    font_size=10;
    record_plot=107;
    central_points_track=(floor(records_L1b/2)-50:floor(records_L1b/2)+50);

    mida = get(0,'ScreenSize');
    mida(3:4)=[1920,1080];
%     mida(3:4)=[1280,720];
    h=figure('Position', mida); 
    
    subplot(2,2,1);
    scatter3(SWATH.lon, SWATH.lat, SWATH.elev,40, SWATH.elev);
    caxis([-500 900]);view(192,44);
    set(gca,'XLim',[min(DEM.lon_deg) max(DEM.lon_deg)]);
    set(gca,'YLim',[min(DEM.lat_deg) max(DEM.lat_deg)]);
    set(gca,'ZLim',[min(min(DEM.height_m)) max(max(DEM.height_m))]);
    colormap(parula); colorbar('Location','southoutside','FontSize',font_size); hold on;
    plot3(L1B.lon,L1B.lat,L1B.alt-L1B.range_ku_l1b_echo,'-k');
    %legend('L2 Swath', 'L1 ground track', 'Location','northeastoutside')
    legend('L2 Swath', 'L1 ground track', 'Location','northeast')
    set(gca, 'FontSize',font_size);
    figlabels('Longitude[deg]','Latitude [deg]','Elevation[m]',{'Swath processing output' ;[' Coherence: theshold ' num2str(coherence_threshold)];['Along track records: ' num2str(records_L1b) ' / Total elevations measured : ' num2str(length(SWATH.lon))]},font_size);
    
    subplot(2,2,2);
    %mesh(DEM.lon_deg,DEM.lat_deg,DEM.height_m); caxis([-500 900]);view(128,34); hold all;
    imagesc(DEM.lon_deg,DEM.lat_deg,DEM.height_m); hold all;
    colorbar('Location','southoutside','FontSize',font_size); hold on;set(gca, 'FontSize',font_size);
    %hold all;  plot3(L1B.lon,L1B.lat,L1B.alt-L1B.range_ku_l1b_echo,'-k');
    hold all;  plot(L1B.lon,L1B.lat,'-k');
    %legend('reference DEM','L1 ground track',  'Location','northeastoutside');
    legend('reference DEM','L1 ground track',  'Location','northeast');
    set(gca,'XLim',[min(DEM.lon_deg) max(DEM.lon_deg)]);
    set(gca,'YLim',[min(DEM.lat_deg) max(DEM.lat_deg)]);
    %set(gca,'ZLim',[min(min(DEM.height_m)) max(max(DEM.height_m))]);
    figlabels('Longitude[deg]','Latitude [deg]','Elevation[m]','reference DEM' ,font_size);

    subplot(2,2,3);
    scatter3(SWATH.lon, SWATH.lat, SWATH.elev,40, SWATH.individual_dem_diff);
    caxis([-2 2]);view(192,44);
    colormap(jet); colorbar('Location','southoutside'); hold on;
    set(gca,'XLim',[min(DEM.lon_deg) max(DEM.lon_deg)]);
    set(gca,'YLim',[min(DEM.lat_deg) max(DEM.lat_deg)]);
    set(gca,'ZLim',[min(min(DEM.height_m)) max(max(DEM.height_m))]);
    plot3(L1B.lon,L1B.lat,L1B.alt-L1B.range_ku_l1b_echo,'-k');
    %legend('DEM difference', 'L1 ground track', 'Location','northeastoutside'); set(gca, 'FontSize',font_size);
    legend('DEM difference', 'L1 ground track', 'Location','northeast'); set(gca, 'FontSize',font_size);
    figlabels('Longitude[deg]','Latitude [deg]','Elevation[m]',{['Difference with DEM'];['mean difference with DEM: ' num2str(mean(SWATH.average_dem_diff(central_points_track))) ' meters,  mean MAD: ' num2str(mean(SWATH.MAD(central_points_track))) ' meters']} ,font_size);
    
    subplot(2,2,4);
    %plot(L1B.lat,SWATH.average_dem_diff,'b');
    plot(SWATH.average_dem_diff,'b');
    %figlabels('Latitude[deg]','[meters]','',{'SWATH results';['mean difference with DEM: ' num2str(mean(SWATH.average_dem_diff(central_points_track))) ' meters,  mean MAD: ' num2str(mean(SWATH.MAD(central_points_track))) ' meters']} ,font_size);
    figlabels('record L1B','[meters]','',{'SWATH results';['mean difference with DEM: ' num2str(mean(SWATH.average_dem_diff(central_points_track))) ' meters,  mean MAD: ' num2str(mean(SWATH.MAD(central_points_track))) ' meters']} ,font_size);
    hold all;
    %plot(L1B.lat,SWATH.MAD,'r'); set(gca, 'FontSize',font_size);
    plot(SWATH.MAD,'r'); set(gca, 'FontSize',font_size);
    plot(L1B.lat(central_points_track),zeros(1,length(central_points_track)) + mean(SWATH.average_dem_diff(central_points_track)),'-.b');
    plot(L1B.lat(central_points_track),zeros(1,length(central_points_track)) + mean(SWATH.MAD(central_points_track)),'-.r');
    legend('mean difference with DEM' , 'mean MAD ', 'Location','northeast');

    img = getframe(h);
    set(h,'PaperPositionMode','auto');
    saveas(h,[plots_folder num2str(coherence_threshold) '_Coh_DEM_and_swath.png'])
    close all

    h=figure('Position', mida); 
    
    subplot(2,2,2);
    scatter(SWATH.lon(SWATH.wvf_number==record_plot), SWATH.elev(SWATH.wvf_number==record_plot),30, SWATH.elev(SWATH.wvf_number==record_plot)); hold on;
    caxis([-500 900]); hold on
    set(gca,'XLim',[min(SWATH.lon(SWATH.wvf_number==record_plot)) max(SWATH.lon(SWATH.wvf_number==record_plot))]);
    set(gca,'YLim',[min(SWATH.elev(SWATH.wvf_number==record_plot)) max(SWATH.elev(SWATH.wvf_number==record_plot))]);
    colormap(parula); colorbar('Location','eastoutside','FontSize',font_size); hold on;
    set(gca, 'FontSize',font_size);
    figlabels('Longitude[deg]','Elevation[m]','',['Swath processing output, record #',num2str(record_plot,'%03d')] ,font_size);
    plot(SWATH.lon(SWATH.wvf_number==record_plot),SWATH.individual_dem_diff(SWATH.wvf_number==record_plot)+SWATH.elev(SWATH.wvf_number==record_plot),'-k','LineWidth',2); caxis([-500 900]); hold all;
    scatter(L1B.lon(record_plot), L1B.alt(record_plot)-L1B.range_ku_l1b_echo(record_plot),30, 'k', 'fill');
    legend('L2 Swath', 'reference DEM','L1b location' , 'Location','southeast')
    
%     scatter(SWATH.lat(SWATH.wvf_number==record_plot), SWATH.elev(SWATH.wvf_number==record_plot),30, SWATH.elev(SWATH.wvf_number==record_plot)); hold on;
%     caxis([-500 900]); hold on
%     set(gca,'XLim',[min(SWATH.lat(SWATH.wvf_number==record_plot)) max(SWATH.lat(SWATH.wvf_number==record_plot))]);
%     set(gca,'YLim',[min(SWATH.elev(SWATH.wvf_number==record_plot)) max(SWATH.elev(SWATH.wvf_number==record_plot))]);
%     colormap(parula); colorbar('Location','southoutside','FontSize',font_size); hold on;
%     set(gca, 'FontSize',font_size);
%     figlabels('Latitude[deg]','Elevation[m]','','Swath processing output, record #100' ,font_size);
%     plot(SWATH.lat(SWATH.wvf_number==record_plot),SWATH.individual_dem_diff(SWATH.wvf_number==record_plot)+SWATH.elev(SWATH.wvf_number==record_plot),'-k','LineWidth',2); caxis([-500 900]); hold all;
%     scatter(L1B.lon(record_plot), L1B.alt(record_plot)-L1B.range_ku_l1b_echo(record_plot),30, 'k', 'fill');
%     legend('L2 Swath', 'reference DEM','L1b location' , 'Location','northeastoutside')
%     
%     scatter3(SWATH.lon(SWATH.wvf_number==record_plot), SWATH.lat(SWATH.wvf_number==record_plot),SWATH.elev(SWATH.wvf_number==record_plot),30, SWATH.elev(SWATH.wvf_number==record_plot)); hold on;

    subplot(2,2,4);
    scatter(SWATH.lon(SWATH.wvf_number==record_plot), SWATH.individual_dem_diff(SWATH.wvf_number==record_plot),30, SWATH.individual_dem_diff(SWATH.wvf_number==record_plot),'fill');
    caxis([-2 2]); colormap(jet); colorbar('Location','eastoutside'); hold on;
    set(gca,'XLim',[min(SWATH.lon(SWATH.wvf_number==record_plot)) max(SWATH.lon(SWATH.wvf_number==record_plot))]);
    legend('DEM difference', 'Location','northeast'); set(gca, 'FontSize',font_size);
    figlabels('Longitude[deg]','DEM diff[m]','',{['Difference with DEM, record #',num2str(record_plot,'%03d') ]; ['mean DEM diff: ' num2str(mean(SWATH.average_dem_diff(record_plot))) ' meters MAD: ' num2str(mean(SWATH.MAD(record_plot))) '  meters']} ,font_size);
    
    surf_samples        = SURF(record_plot).surf_ini_sample(end):SURF(record_plot).surf_end_sample(end);
    
    subplot(3,2,1);
    plot(1:chd.N_samples_sar,10*log10(L1B.scaled_waveforms(record_plot,:)),'k','LineWidth',4); hold on;
    plot(surf_samples,10*log10(L1B.scaled_waveforms(record_plot,surf_samples)),'g','LineWidth',1);
    figlabels('sample','Power [dB]','',['#',num2str(record_plot,'%03d'),' L1B Waveform, #',num2str(SURF(record_plot).N_surfs_OK,'%d'),' Reflections retracked'] ,font_size); set(gca,'XLim',[1 chd.N_samples_sar],'FontSize',font_size);
    
    subplot(3,2,3);
    plot(1:chd.N_samples_sar,L1B.coherence_meas_ku_l1b_echo(record_plot,:),'k','LineWidth',4);hold on;
    plot(surf_samples,L1B.coherence_meas_ku_l1b_echo(record_plot,surf_samples),'g','LineWidth',1);
    figlabels('samples','','',['#',num2str(record_plot,'%03d'),' Coherence: theshold ' num2str(coherence_threshold)],font_size); set(gca,'XLim',[1 chd.N_samples_sar],'FontSize',font_size);set(gca,'YLim',[0 1],'FontSize',font_size);
    
    subplot(3,2,5);     
    AoA    = ((chd.wv_length * unwrap(L1B.phase_diff_meas_ku_l1b_echo(record_plot,:))./(2*pi*B)-roll)/pi*180);   %[deg]
    plot(1:chd.N_samples_sar,AoA,'k','LineWidth',4);hold on;
    plot(surf_samples,AoA,'g','LineWidth',1);    
    figlabels('samples','AoA [deg]','',['#',num2str(record_plot,'%03d'), ' Angle of Arrival - Roll'],font_size); set(gca,'XLim',[1 chd.N_samples_sar],'FontSize',font_size); set(gca,'YLim',[-0.6*3 0.6*3],'FontSize',font_size);
    img = getframe(h);
    set(h,'PaperPositionMode','auto');
    saveas(h,[plots_folder num2str(coherence_threshold) '_Coh_' num2str(record_plot,'%03d') 'record_DEM_and_swath.png'])
             close all
    
    h=figure('Position', mida);         
    subplot(2,1,1)
    histogram(SWATH.individual_dem_diff,-2:0.1:2);
    figlabels('','#','',{'Swath processing output' ;[' Coherence: theshold ' num2str(coherence_threshold)];['Along track records: ' num2str(records_L1b) ' / Total elevations measured : ' num2str(length(SWATH.lon))]},font_size);
    set(gca, 'FontSize',font_size);
    img = getframe(h);
    set(h,'PaperPositionMode','auto');
    saveas(h,[plots_folder num2str(coherence_threshold) '_Coh_' num2str(record_plot,'%03d') 'histogram_individual_dem_diffs.png'])
    close all
    
    



end

filesBulk.SURF  = SURF;

end

%% EXTRA visualisation plots

% %TLI1v6 vertexes:
% %      - Glacier vertexes
%         TLI1_glacier_vertex1=[-0.53464893 -1.62505 -655]; %[lat lon alt] [rad, m]
%         TLI1_glacier_vertex2=[-0.53541752 -1.61605 +714]; %[lat lon alt] [rad, m]
%         TLI1_glacier_vertex3=[-0.53527015 -1.61603 +714]; %[lat lon alt] [rad, m]
%         TLI1_glacier_vertex4=[-0.53449828 -1.62503 -655]; %[lat lon alt] [rad, m]
% %      - Left of glacier mountain vertexes
%         TLI1_leftmount_vertex1=[-0.5409118  -1.62579 -653.46]; %[lat lon alt] [rad, m]
%         TLI1_leftmount_vertex2=[-0.54168306 -1.61676  715.24]; %[lat lon alt] [rad, m]
%         TLI1_leftmount_vertex3=[-0.5354391  -1.61605  715.90]; %[lat lon alt] [rad, m]
%         TLI1_leftmount_vertex4=[-0.53465205 -1.62505 -650.01]; %[lat lon alt] [rad, m]
% %      - Right of glacier mountain vertexes
%         TLI1_rightmount_vertex1=[-0.53449083 -1.62502 -648.89]; %[lat lon alt] [rad, m]
%         TLI1_rightmount_vertex2=[-0.53526465 -1.61603  714.60]; %[lat lon alt] [rad, m]
%         TLI1_rightmount_vertex3=[-0.52908592 -1.61534  712.92]; %[lat lon alt] [rad, m]
%         TLI1_rightmount_vertex4=[-0.52831206 -1.6243  -655.02]; %[lat lon alt] [rad, m]

% %TLI2v2 vertexes:
% %      - Glacier vertexes
%         TLI1_glacier_vertex1=[-0.53464893 -1.62505 -216.6]; %[lat lon alt] [rad, m]
%         TLI1_glacier_vertex2=[-0.53541752 -1.61605 +240.1]; %[lat lon alt] [rad, m]
%         TLI1_glacier_vertex3=[-0.53527015 -1.61603 +240.0]; %[lat lon alt] [rad, m]
%         TLI1_glacier_vertex4=[-0.53449828 -1.62503 -216.4]; %[lat lon alt] [rad, m]
% %      - Left of glacier mountain vertexes
%         TLI1_leftmount_vertex1=[-0.5409118  -1.62579 -216.1]; %[lat lon alt] [rad, m]
%         TLI1_leftmount_vertex2=[-0.54168306 -1.61676  238.4]; %[lat lon alt] [rad, m]
%         TLI1_leftmount_vertex3=[-0.5354391  -1.61605  243.1]; %[lat lon alt] [rad, m]
%         TLI1_leftmount_vertex4=[-0.53465205 -1.62505 -214.3]; %[lat lon alt] [rad, m]
% %      - Right of glacier mountain vertexes
%         TLI1_rightmount_vertex1=[-0.53449083 -1.62502 -214.1]; %[lat lon alt] [rad, m]
%         TLI1_rightmount_vertex2=[-0.53526465 -1.61603  243.1]; %[lat lon alt] [rad, m]
%         TLI1_rightmount_vertex3=[-0.52908592 -1.61534  245]; %[lat lon alt] [rad, m]
%         TLI1_rightmount_vertex4=[-0.52831206 -1.6243  -212]; %[lat lon alt] [rad, m]

% %TLI3v2 vertexes:
% %      - Glacier vertexes
%         TLI1_glacier_vertex1=[-0.53487716 -1.62245 -521.4 ]; %[lat lon alt] [rad, m]
%         TLI1_glacier_vertex2=[-0.53586218 -1.61073 +3046.2]; %[lat lon alt] [rad, m]
%         TLI1_glacier_vertex3=[-0.53571207 -1.61072 +3046.9]; %[lat lon alt] [rad, m]
%         TLI1_glacier_vertex4=[-0.53472346 -1.62243 -521.4 ]; %[lat lon alt] [rad, m]
% %      - Left of glacier mountain vertexes
%         TLI1_leftmount_vertex1=[-0.54114247 -1.62317 -515.1 ]; %[lat lon alt] [rad, m]
%         TLI1_leftmount_vertex2=[-0.54212994 -1.61142 +3050.8]; %[lat lon alt] [rad, m]
%         TLI1_leftmount_vertex3=[-0.53586805 -1.61073 +3048.2]; %[lat lon alt] [rad, m]
%         TLI1_leftmount_vertex4=[-0.53488015 -1.62245 -516.7 ]; %[lat lon alt] [rad, m]
% %      - Right of glacier mountain vertexes
%         TLI1_rightmount_vertex1=[-0.53471782 -1.62243 -519.1]; %[lat lon alt] [rad, m]
%         TLI1_rightmount_vertex2=[-0.5357087  -1.61072 3045.0]; %[lat lon alt] [rad, m]
%         TLI1_rightmount_vertex3=[-0.53586805 -1.61073 +3048.2]; %[lat lon alt] [rad, m]
%         TLI1_rightmount_vertex4=[-0.53488015 -1.62245 -516.7 ]; %[lat lon alt] [rad, m]

% TLI17 vertexes:
%      - Mountain vertexes
%         TLI1_leftmount_vertex1=[-0.54095 -1.62615 -707.8]; %[lat lon alt] [rad, m]
%         TLI1_leftmount_vertex2=[-0.54177 -1.61676  715.8]; %[lat lon alt] [rad, m]
%         TLI1_leftmount_vertex3=[-0.52908 -1.61534  714.9]; %[lat lon alt] [rad, m]
%         TLI1_leftmount_vertex4=[-0.52828 -1.62465 -708.1]; %[lat lon alt] [rad, m]
% % 
% % 
% figure;
% Plot glacier and mountain slope
% patch([TLI1_leftmount_vertex1(1) TLI1_leftmount_vertex2(1) TLI1_rightmount_vertex3(1) TLI1_rightmount_vertex4(1)]*180/pi, [TLI1_leftmount_vertex1(2) TLI1_leftmount_vertex2(2) TLI1_rightmount_vertex3(2) TLI1_rightmount_vertex4(2)]*180/pi, [TLI1_leftmount_vertex1(3) TLI1_leftmount_vertex2(3) TLI1_rightmount_vertex3(3) TLI1_rightmount_vertex4(3)],'green');
% 
% %Just for TLI.17, plot the whole slope
% patch([TLI1_leftmount_vertex1(1) TLI1_leftmount_vertex2(1) TLI1_leftmount_vertex3(1) TLI1_leftmount_vertex4(1)]*180/pi, [TLI1_leftmount_vertex1(2) TLI1_leftmount_vertex2(2) TLI1_leftmount_vertex3(2) TLI1_leftmount_vertex4(2)]*180/pi, [TLI1_leftmount_vertex1(3) TLI1_leftmount_vertex2(3) TLI1_leftmount_vertex3(3) TLI1_leftmount_vertex4(3)], 'green')

% 
% % Use these 3. 1st plot one mountain side and the glacier, and after all wraps the other mountain side, so we don't get 2 mountain data tags
% patch([TLI1_glacier_vertex1(1) TLI1_glacier_vertex2(1) TLI1_glacier_vertex3(1) TLI1_glacier_vertex4(1)]*180/pi,           [TLI1_glacier_vertex1(2) TLI1_glacier_vertex2(2) TLI1_glacier_vertex3(2) TLI1_glacier_vertex4(2)]*180/pi,           [TLI1_glacier_vertex1(3) TLI1_glacier_vertex2(3) TLI1_glacier_vertex3(3) TLI1_glacier_vertex4(3)],           'cyan')
% patch([TLI1_leftmount_vertex1(1) TLI1_leftmount_vertex2(1) TLI1_glacier_vertex2(1) TLI1_glacier_vertex1(1)]*180/pi, [TLI1_leftmount_vertex1(2) TLI1_leftmount_vertex2(2) TLI1_glacier_vertex2(2) TLI1_glacier_vertex1(2)]*180/pi, [TLI1_leftmount_vertex1(3) TLI1_leftmount_vertex2(3) TLI1_glacier_vertex2(3) TLI1_glacier_vertex1(3)], 'green')
% patch([TLI1_glacier_vertex4(1) TLI1_glacier_vertex3(1) TLI1_rightmount_vertex3(1) TLI1_rightmount_vertex4(1)]*180/pi, [TLI1_glacier_vertex4(2) TLI1_glacier_vertex3(2) TLI1_rightmount_vertex3(2) TLI1_rightmount_vertex4(2)]*180/pi, [TLI1_glacier_vertex4(3) TLI1_glacier_vertex3(3) TLI1_rightmount_vertex3(3) TLI1_rightmount_vertex4(3)], 'green')
% 
% hold on
% %Plot satellite track
% plot3(L1B.lat,L1B.lon,L1B.alt-L1B.range_ku_l1b_echo+-850,'-k');
% legend('Mountain', 'Glacier', 'Satellite Track') %legend without tracks
% legend('Glacier', 'Mountain', 'Satellite Track', 'L1B record 124','L1B record 125','L1B record 126','L1B record 127','L1B record 128')
% %Plot the current record and wrap
%  scatter3(SURF_wrap(index_L1b).lat(i_wrap,surf_samples),SURF_wrap(index_L1b).lon(i_wrap,surf_samples),Elev_wrap(i_wrap,surf_samples))
% % % 
% % 

