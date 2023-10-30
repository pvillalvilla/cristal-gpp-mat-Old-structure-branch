%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L.
% --------------------------------------------------------
% CryoSat-2
% This compute the SWATH L2 for SARIn L1B Cryosat porducts
% ---------------------------------------------------------
% Objective:
%
% INPUTs : - CR2 L1B
% OUTPUTs: - Video, stadistics
%
% ----------------------------------------------------------
% Author:    Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (11/09/2013)
% v1.0  First version
% v2.0  Roll sign and variable, CR2 Baseline C, unwrapping within same surf
% v3.0  Make it a function
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MultiPeak_retracker_func('F:\MoGLA\Central_Asia\', 0, 1, 1)
% MultiPeak_retracker_func('F:\CR2+Forest\', 0, 1, 1)
% MultiPeak_retracker_func('/data/M xOGLA/DATA/Antarctic_peninsula',2010, 0, 1, 0)

function MultiPeak_Swath_Processing(Path,inputyear, plotting, kml_flag, ambiguity_flag)

swath_processing_flag=1;
plotting_DEM=1;
roll_bias_TRP=0.0072; %[deg]
elev_limit_plot= 7500;
elev_limit_min_plot= 2000;

cd (Path);
separators=find(Path=='/');
if(isempty(separators))
    separators=find(Path=='\');
end


inputdata=Path(separators(end)+1:end);
% kml_flag=1;
% plotting = 1;
% ambiguity_flag      = 0; % generate SURF_lef and SURF_right with one wrap.



font_size = 8;
coherence_threshold = 600;
power_threshold=-170;

geometric_factor = 1.113; % Geometric factor due to earth curvature


surf_distance       = 3; % distance in samples to consider a differencet surface from two consecutives high coherence periods
surf_threshold      = 2; % minimum number of samples for a surface to be considered valid
peak_threshold      = 0; % The scaling factor has to be used, not yet included.
N_wraps             = 3; % wraps for the DEM of the plot comparison
zp_factor           = 2;
maskFlag            = 0;
c_cst               = 299792458; %[meters]
T0_chd              = 1/320e6;
N_samples_sar_chd   = 1024;
ellipsoid           =   referenceEllipsoid ('wgs84');
flat_coeff_cst      = 0.003352810664747;
semi_major_axis_cst = 6378137;
wv_length_ku        = 0.022084;  %[meters]
B                   = 1.1676;    %[meters]


WaterColor = [143 226 255]/255;
mida = get(0,'ScreenSize');
mida(3:4)=[1920,1080];
Baltoro=[35.736389, 76.380833];
files=dir(['./inputs_' inputyear '/*.DBL']);
N_files=length(files);
maskfilename = dir('./mask/*.kml');

if(~isempty(maskfilename))
    mask  = kml2lla(['./mask/' maskfilename.name]);
    maskFlag=1;
end

for i_file=2:N_files
    processing_time=tic;
    L1B_FILE=['./inputs_' inputyear '/' files(i_file).name];
    %disp(L1B_FILE(29:end-4));
    [lat,lon, ~ ,records_num] = auto_readL1b_all_modes_C_lat_lon(L1B_FILE,1,0);
    if(maskFlag)
        
%         [lat, lon] = adapaptC2L1_lat_lon(out_L1B_,records_num);
        index_inside = find(inpolygon(lon,lat,mask.coord(:,1),mask.coord(:,2))); % records inside the given mask
        if(isempty(index_inside))
            fclose all;
            movefile(L1B_FILE,['./skip/']);
            continue;
        else
            init_record = max(1,(floor(index_inside(1)/20)));
            end_record = min(ceil(max(index_inside)/20),records_num);
        end
        
    else
        init_record = 1;
        end_record = (records_num);
    end
    aux=dir(['./results/' L1B_FILE(end-50:end-4) '*.mat']);
    if(~isempty(aux))
        if(strcmp(aux.name,[L1B_FILE(end-50:end-4) '.mat']))
            load(['./results/' L1B_FILE(end-50:end-4) '.mat']);
        else
            [out_L1B, records_num_] = auto_readL1b_all_modes_C(L1B_FILE,1,0, init_record,end_record);
            
            save(['./results/' L1B_FILE(end-50:end-4) '.mat'],'out_L1B','init_record','end_record');
        end
    else
        [out_L1B, records_num_] = auto_readL1b_all_modes_C(L1B_FILE,1,0, init_record,end_record);
        
        save(['./results/' L1B_FILE(end-50:end-4) '.mat'],'out_L1B','init_record','end_record');
        
    end
    
    if(plotting)
        h=figure('Position',mida); subplot(3,4,3);
    end
    
    records_num=end_record-init_record;
    clear orbit_track ground_track index index2 geo_corr alt scale_factor plot_Wvf
    clear plot_Coh
    SURF = [];
    SURF_left = [];
    SURF_right = [];
    lat=[];
    lon=[];
    fineH=[];
    stdAoA = [];
    peakness = [];
    slopeAoA = [];
    slopeElev = [];
    lon_dem_index= [];
    lat_dem_index=[];
    SWATH=struct('lat',[],'lon',[],'elev',[],'N_selected',[],...
        'individual_dem_diff',[],'average_dem_diff',[],'std_dem_diff',[],...
        'wvf_number',[],'sample',[]);
    
                                   
    if(plotting)
        mkdir(['./results/' L1B_FILE(end-50:end-4)]);
    end
    
    for i_record=1: records_num
        %% 1.EXTRACT L1B Parameters
        b = out_L1B(i_record).corrections_group;
        geo_corr_aux = geo_corrections_isard(b);
        for k_record = 1:20
            index           =out_L1B(i_record).time_orbit_group(1,k_record).burst_counter;
            index2          =((i_record-1)*20+k_record);
            geo_corr(index2)= geo_corr_aux; % computed for land ice
            clear plot_Wvf plot_AoA plot_Coh 
            elev_average_error=[];
            ambiguity_wrap=[];
            elevation_std=[];
            surf_index=[];
            surf_samples_length=[];
            init_sample=[];
            alt(index2)     = out_L1B(1,i_record).time_orbit_group(1,k_record).altitude_COG*1e-3-out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12*c_cst/2-geo_corr(index2)*1e-3;
            scale_factor    = out_L1B(1,i_record).waveform_group(1,k_record).echo_scale_factor;
            scale_power     = out_L1B(1,i_record).waveform_group(1,k_record).echo_scale_power;
            plot_Wvf        = 10.*log10(out_L1B(1,i_record).waveform_group(1,k_record).averaged_power_echo_waveform.*(scale_factor.*10^-9).*2.^scale_power);
            plot_Wvf_lin    = (out_L1B(1,i_record).waveform_group(1,k_record).averaged_power_echo_waveform.*(scale_factor.*10^-9).*2.^scale_power)./max(out_L1B(1,i_record).waveform_group(1,k_record).averaged_power_echo_waveform.*(scale_factor.*10^-9).*2.^scale_power);
            plot_Coh        = out_L1B(1,i_record).waveform_group(1,k_record).coherence;
            plot_Coh(plot_Coh>1000)= 0;
            plot_Coh_smooth = smooth(plot_Coh);
            orbit_track(index2,1).Geometry  ='Point';
            orbit_track(index2,1).Lon       = out_L1B(1,i_record).time_orbit_group(1,k_record).longitude*1e-7;
            orbit_track(index2,1).Lat       = out_L1B(1,i_record).time_orbit_group(1,k_record).latitude*1e-7;
            orbit_track(index2,1).Alt       = out_L1B(1,i_record).time_orbit_group(1,k_record).altitude_COG*1e-3;
            ground_track(index2,1).Geometry ='Point';
            if(index2==1)
                ground_track(index2+1,1).Geometry ='Point';
                ground_track(index2+1,1).Lon      = out_L1B(1,i_record).time_orbit_group(1,k_record+1).longitude*1e-7;
                ground_track(index2+1,1).Lat      = out_L1B(1,i_record).time_orbit_group(1,k_record+1).latitude*1e-7;
                ground_track(index2+1,1).Alt      = out_L1B(1,i_record).time_orbit_group(1,k_record+1).altitude_COG*1e-3-out_L1B(1,i_record).measurements_group(1,k_record+1).window_delay*1e-12*c_cst/2;
            end
            ground_track(index2,1).Lon      = out_L1B(1,i_record).time_orbit_group(1,k_record).longitude*1e-7;
            ground_track(index2,1).Lat      = out_L1B(1,i_record).time_orbit_group(1,k_record).latitude*1e-7;
            ground_track(index2,1).Alt      = out_L1B(1,i_record).time_orbit_group(1,k_record).altitude_COG*1e-3-out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12*c_cst/2;
            
            % ROLL for baseline B, bias not included
            %             roll          = -atan(out_L1B(i_record).time_orbit_group(k_record).inferometer_baseline(1,1)/out_L1B(i_record).time_orbit_group(k_record).inferometer_baseline(1,3)); % [rad]
            roll            = -(out_L1B(i_record).time_orbit_group(k_record).bench_roll*1e-7+roll_bias_TRP)/180*pi; %[10^-1 microdeg to radians]
            yaw             = atan(out_L1B(i_record).time_orbit_group(k_record).inferometer_baseline(1,2)/out_L1B(i_record).time_orbit_group(k_record).inferometer_baseline(1,3)); % [rad]
            ib              = out_L1B(i_record).time_orbit_group(k_record).inferometer_baseline*1e-6;
            rb              = out_L1B(i_record).time_orbit_group(k_record).real_beam_direction*1e-6;
            
            
            %% 2.COMPUTE X_AXIS
            top_elev  = alt(index2)+N_samples_sar_chd/2*T0_chd/zp_factor*c_cst/2;
            bottom_elev    = alt(index2)-(N_samples_sar_chd/2)*T0_chd/zp_factor*c_cst/2;
            far_Range = out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12*c_cst/2+N_samples_sar_chd/2*T0_chd/zp_factor*c_cst/2;
            short_Range = out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12*c_cst/2-(N_samples_sar_chd/2)*T0_chd/zp_factor*c_cst/2;
            
            x_axisL1 = (top_elev-T0_chd/zp_factor*c_cst/2:-T0_chd/zp_factor*c_cst/2:bottom_elev);
            plot_Range = (short_Range:T0_chd/zp_factor*c_cst/2:far_Range-T0_chd/zp_factor*c_cst/2);
            phase_difference =out_L1B(1,i_record).waveform_group(1,k_record).phase_difference/1e6;
            plot_AoA    = (wv_length_ku * phase_difference./(2*pi*B)-roll)/pi*180;   %[deg]
            
            plot_AoA_smooth = smooth(plot_AoA);
            plot_Xdist  =  plot_AoA /180*pi .* plot_Range;
            jump_AoA     = (wv_length_ku * 2*pi./(2*pi*B))/pi*180;   %[deg] 360º to AoA
            
            plot_AoA_smooth = smooth(plot_AoA);
            th_AoA     = (wv_length_ku * 3/2*pi./(2*pi*B))/pi*180;   %[deg] theoretical Angle of Arrival, to check for jumps
            
            %% COMPUTE RANGE and uncorrected ELEVATIONS
            plot_Slant_Range    = sqrt(plot_Xdist.^2 + plot_Range.^2);
            
            plot_Slant_Elev     = -plot_Slant_Range+ out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12*c_cst/2+alt(index2);
            plot_Elev   = cos(plot_AoA.*pi./180).*plot_Slant_Elev; % need to be corrected for delta_Ha and delta_Hp: SEE TOM ARMITAGE PAPER https://www.researchgate.net/publication/260623282_Using_the_Interferometric_Capabilities_of_the_ESA_CryoSat-2_Mission_to_Improve_the_Accuracy_of_Sea_Ice_Freeboard_Retrievals
           
            if(index2>1)
                v =[ground_track(index2,1).Lat,ground_track(index2,1).Lon]-[ground_track(index2-1,1).Lat,ground_track(index2-1,1).Lon];
                [arclen_v,az_v] = distance('gc',ground_track(index2-1,1).Lat,ground_track(index2-1,1).Lon,ground_track(index2,1).Lat,ground_track(index2,1).Lon,[ellipsoid.SemimajorAxis ellipsoid.Eccentricity]);
            else
                v =[ground_track(index2+1,1).Lat,ground_track(index2+1,1).Lon]-[ground_track(index2,1).Lat,ground_track(index2,1).Lon];
                [arclen_v,az_v] = distance('gc',ground_track(index2,1).Lat,ground_track(index2,1).Lon,ground_track(index2+1,1).Lat,ground_track(index2+1,1).Lon,[ellipsoid.SemimajorAxis ellipsoid.Eccentricity]);
                
            end
            %az_v along track direction foreward
            %az_d always vector perpendicular to de track to the right side of the along track direction
            
            az_d=az_v-90;
            
%             if(v(1)>0) % ascending pass
%                 az_d=az_v-90;
%             elseif (v(1)<0) % descending pass
%                 az_d=az_v-90;
%             end

            [lats,lons] = reckon(ground_track(index2,1).Lat,ground_track(index2,1).Lon,plot_Xdist,az_d,[ellipsoid.SemimajorAxis ellipsoid.Eccentricity]);
            if(swath_processing_flag)
                AoA_points= -jump_AoA*(N_wraps):jump_AoA:jump_AoA*(N_wraps);% [-jump_AoA jump_AoA] = case N_wraps = 1;
%                 AoA_points = AoA_points(AoA_points~=0);
                phase_wraps=-2*pi*N_wraps:2*pi:2*pi*N_wraps;
                phase_wraps = phase_wraps(phase_wraps~=0);
                Ambiguities      = phase_wraps * wv_length_ku * out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12*c_cst/2./(2*pi*B);
                for i_wrap=1: length(AoA_points)
                    plot_Xdist_wrap  =  (plot_AoA+AoA_points(i_wrap)) /180*pi .* plot_Range;
%                     plot_Slant_Range_wrap    = sqrt(plot_Xdist_wrap.^2 + plot_Range.^2);
                    plot_Slant_Elev_wrap(i_wrap,:)     = -plot_Range+ out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12*c_cst/2+alt(index2);
                    plot_Elev_wrap(i_wrap,:)   = cos((plot_AoA+AoA_points(i_wrap)).*pi./180).*plot_Slant_Elev_wrap(i_wrap,:); % need to be corrected for delta_Ha and delta_Hp: SEE TOM ARMITAGE PAPER https://www.researchgate.net/publication/260623282_Using_the_Interferometric_Capabilities_of_the_ESA_CryoSat-2_Mission_to_Improve_the_Accuracy_of_Sea_Ice_Freeboard_Retrievals
                                                   
                    [SURF_wrap(index2).lat(i_wrap,:),SURF_wrap(index2).lon(i_wrap,:)]   = reckon(ground_track(index2,1).Lat,ground_track(index2,1).Lon, plot_Xdist_wrap,az_d,[ellipsoid.SemimajorAxis ellipsoid.Eccentricity]);
                    
                end
            end
            
            AoA_pos = find(abs(plot_AoA)<0.01);
            AoA_reflex = find(abs(diff(diff(plot_AoA)))<0.005);
            
            if(plotting)
                subplot(4,3,1);plot(x_axisL1,plot_Wvf,'k','LineWidth',2);hold on;
                subplot(4,3,4);plot(1:N_samples_sar_chd,plot_Coh,'k','LineWidth',2);hold on;
%                 plot(1:N_samples_sar_chd,plot_Coh_smooth,'b','LineWidth',1);
                subplot(4,3,7);plot(1:N_samples_sar_chd,plot_AoA,'k','LineWidth',2);hold on;
                
                
            end
            if(swath_processing_flag)
                
                min_lat=min(min(SURF_wrap(index2).lat(:,:)))-0.25;
                max_lat=max(max(SURF_wrap(index2).lat(:,:)))+0.25;
                min_lon=min(min(SURF_wrap(index2).lon(:,:)));
                max_lon=max(max(SURF_wrap(index2).lon(:,:)));
                % Compute the DEM diff between all the points to select the best ambiguity wrap N
                DEM=readhgt([min_lat,max_lat,min_lon,max_lon]);
                lon_step=mean(diff(DEM.lon));
                lat_step=mean(diff(DEM.lat));
                [lon_grid,lat_grid]=meshgrid(min(DEM.lon):lon_step:max(DEM.lon),min(DEM.lat):lat_step:max(DEM.lat));
                DEM.z= double(DEM.z);
                DEM.z((DEM.z==-32768))=NaN;
                iDEMgood=~(isnan(DEM.z));
                %build new DEM interpolated using 0.0005 -> 55 meters
                DEM_resolution = 0.0005;
                [lon_hd,lat_hd]= meshgrid(min(DEM.lon):DEM_resolution:max(DEM.lon),min(DEM.lat):DEM_resolution:max(DEM.lat));
%                 DEM_hd= griddatta( lon_grid(iDEMgood),lat_grid(iDEMgood),DEM.z(iDEMgood), lon_hd, lat_hd) ;
                F = scatteredInterpolant( lon_grid(iDEMgood),lat_grid(iDEMgood),DEM.z(iDEMgood)) ;
                DEM_hd = F(lon_hd,lat_hd);
                 
                if(plotting_DEM)
                    %                         hh=figure;
                    subplot(4,3,[2 3 5 6]); 
                    dem_2D=imagesc(lon_hd(1,:),(lat_hd(:,1)),DEM_hd);
                    demcmap(DEM_hd);
                    uistack(dem_2D,'bottom');
                    set(gca,'YDir','normal');
                    hold on;
                    
%                     set(gca,'Color',WaterColor);
%                     set(k, 'edgecolor','none');set(k, 'FaceAlpha',0.75);
                    
                    
                end
                
            end
            
            Coh_pos = find((plot_Coh>coherence_threshold).*(plot_Coh<1010));
            [max_values,max_samples,width]=findpeaks(plot_Coh,1:1024,'WidthReference','halfheight');
            if(~isempty(Coh_pos))
                
                plot_diff_Coh = diff(plot_Coh_smooth);
                
                %surf_trans_pos find the jumps between different surfaces
                surf_trans_pos    = find(diff(Coh_pos)>surf_distance);
                %             if(max(Coh_pos)==N_samples_sar_chd)
                %                 [~,surf_trans_pos(length(surf_trans_pos)+1)]=max(Coh_pos);
                %             end
                
                % the number of different surfaces is then the number of jumps+1
                SURF(index2).N_surfs 	= length(surf_trans_pos)+1;
                clear  surf_samples point_color N_selected mean_DEM_error wrap_selected
                SURF(index2).peakness_surf = zeros (1,SURF(index2).N_surfs);
                
                i_surf_index=0;
                SURF(index2).N_surfs_OK=0;
                
                
                for i_surf=1: SURF(index2).N_surfs
                    i_surf_index = i_surf_index+1;
                    SURF(index2).POCA=0;
                    SURF(index2).Wrap_flag(i_surf_index)=0;
                    if(i_surf==1)
                        SURF(index2).surf_ini_sample(i_surf_index)    = Coh_pos(1);
                    else
                        SURF(index2).surf_ini_sample(i_surf_index)    = Coh_pos(surf_trans_pos(i_surf-1)+1);
                    end
                    
                    if(i_surf==SURF(index2).N_surfs)
                        SURF(index2).surf_end_sample(i_surf_index)    = Coh_pos(end);
                    else
                        SURF(index2).surf_end_sample(i_surf_index)    = Coh_pos(surf_trans_pos(i_surf));
                    end
                    
                    surf_samples        = SURF(index2).surf_ini_sample(i_surf_index):SURF(index2).surf_end_sample(i_surf_index);
                    
                     
                    if(sum(isfinite(plot_Wvf(surf_samples)))>0) %avoid infinites from the log10 conversion
                        [wvf_min,wvf_min_pos] = min(plot_Wvf_lin(surf_samples(isfinite(plot_Wvf_lin(surf_samples)))));
                        [surf_max,surf_max_pos] = max(plot_Wvf_lin(surf_samples(isfinite(plot_Wvf_lin(surf_samples)))));
                        SURF(index2).peakness_surf(i_surf_index) = surf_max/sum(plot_Wvf_lin(surf_samples(isfinite(plot_Wvf_lin(surf_samples)))));
                        
                        
                        %% CONDITIONS
                        %                         if(SURF(index2).peakness_surf(i_surf_index)> peak_threshold && (length(surf_samples)>surf_threshold))
                        if((length(surf_samples)>surf_threshold)&&(max(plot_Wvf(surf_samples))>power_threshold))
                            SURF(index2).N_surfs_OK = SURF(index2).N_surfs_OK+1;
                            if(SURF(index2).N_surfs_OK==1)
                                SURF(index2).POCA=i_surf_index;
                            end
                            Power_surf_pondered=0;
                            
                            for kk=SURF(index2).surf_ini_sample(i_surf_index):SURF(index2).surf_end_sample(i_surf_index);
                                if(isfinite(plot_Wvf((kk))))
                                    Power_surf_pondered = Power_surf_pondered + kk*plot_Wvf((kk));
                                end
                            end
                            %Check jumps in the phase related with warps
                            
                            jumps_phase = abs(diff(plot_AoA_smooth(surf_samples)))>th_AoA;
                            if(~isempty(find(jumps_phase)))
                                [min_pks,min_locs]=findpeaks(-cumsum(plot_AoA_smooth(surf_samples)));
                                [max_pks,max_locs]=findpeaks(cumsum(plot_AoA_smooth(surf_samples)));
                                plot_AoA_unwrapped(surf_samples) = plot_AoA(surf_samples);
                                offset_samples=surf_samples(1)-1;
                                
                                if(~isempty(max_locs)&&~isempty(min_locs))%multiple jumps due to noise
                                    if(length(min_locs)>length(max_locs))% decreasing phase jumping to positives values
                                        for i_jump=1:length(max_locs)
                                            samples_to_unwrap = offset_samples+(min_locs(i_jump)+1):offset_samples+max_locs(i_jump);
                                            plot_AoA_unwrapped(samples_to_unwrap) = -jump_AoA+plot_AoA(samples_to_unwrap);
                                        end
                                        samples_to_unwrap = (offset_samples+min_locs(i_jump+1)+1):surf_samples(end);
                                        plot_AoA_unwrapped(samples_to_unwrap) = -jump_AoA+plot_AoA(samples_to_unwrap);
                                        
                                    elseif(length(min_locs)<length(max_locs)) % increasing phase jumping to negative values
                                        
                                        for i_jump=1:length(min_locs)
                                            samples_to_unwrap = offset_samples+(max_locs(i_jump)+1):offset_samples+min_locs(i_jump);
                                            plot_AoA_unwrapped(samples_to_unwrap) = jump_AoA+plot_AoA(samples_to_unwrap);
                                        end
                                        samples_to_unwrap = (offset_samples+max_locs(i_jump+1)+1):surf_samples(end);
                                        plot_AoA_unwrapped(samples_to_unwrap) = jump_AoA+plot_AoA(samples_to_unwrap);
                                    else
                                        
                                        
                                    end
                                else % single jump
                                    if(~isempty(max_locs))% decreasing phase jumping to positives values
                                        
                                        samples_to_unwrap = (offset_samples+max_locs+1):surf_samples(end);
                                        plot_AoA_unwrapped(samples_to_unwrap) = jump_AoA+plot_AoA(samples_to_unwrap);
                                        
                                    elseif(~isempty(min_locs)) % increasing phase jumping to negative values
                                        
                                        
                                        samples_to_unwrap = (offset_samples+min_locs+1):surf_samples(end);
                                        plot_AoA_unwrapped(samples_to_unwrap) = -jump_AoA+plot_AoA(samples_to_unwrap);
                                        
                                    end
                                end
                                SURF(index2).Wrap_flag(i_surf_index)= 1;
                                plot_AoA(surf_samples) = plot_AoA_unwrapped(surf_samples);
                                plot_Xdist(surf_samples)  =  plot_AoA(surf_samples) /180*pi * out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12*c_cst/2;
                            end
                            Power_surf_sum = sum(plot_Wvf(surf_samples(isfinite(plot_Wvf(surf_samples)))));
                            
                            SURF(index2).CoG_pos(i_surf_index) = Power_surf_pondered/Power_surf_sum;
                            [~,SURF(index2).Peak_pos(i_surf_index)]= max(plot_Wvf(surf_samples));
                            SURF(index2).Peak_pos(i_surf_index) = SURF(index2).surf_ini_sample(i_surf_index) + SURF(index2).Peak_pos(i_surf_index)-1;
                            SURF(index2).Peak_AoA(i_surf_index) = plot_AoA(SURF(index2).Peak_pos(i_surf_index));
                            SURF(index2).Peak_Coh(i_surf_index) = plot_Coh(SURF(index2).Peak_pos(i_surf_index));
                            SURF(index2).Peak_Wvf(i_surf_index) = plot_Wvf(SURF(index2).Peak_pos(i_surf_index));
                            SURF(index2).std_AoA(i_surf_index)  = std(plot_AoA(surf_samples));
                            SURF(index2).mean_AoA(i_surf_index)  = mean(plot_AoA(surf_samples));
                            coefficients = polyfit(surf_samples, plot_AoA(surf_samples), 1);
                            SURF(index2).slopeAoA(i_surf_index) = coefficients(1);
                            SURF(index2).mean_Coh(i_surf_index)  = mean(plot_Coh(surf_samples));
                            SURF(index2).wvf_number = index;
                            
                            SURF(index2).Xdistance(i_surf_index)        = plot_Xdist(round(SURF(index2).CoG_pos(i_surf_index)));
                            SURF(index2).h_COG(i_surf_index)            = plot_Elev(round(SURF(index2).CoG_pos(i_surf_index)));
                            coefficients2                           = polyfit(plot_Xdist(surf_samples),plot_Elev(surf_samples), 1);
                            SURF(index2).slopeElev(i_surf_index)    = coefficients2(1);
                            
                            
                            [SURF(index2).lat(i_surf_index),SURF(index2).lon(i_surf_index)] = reckon(ground_track(index2,1).Lat,ground_track(index2,1).Lon,plot_Xdist(round(SURF(index2).CoG_pos(i_surf_index))),az_d,[ellipsoid.SemimajorAxis ellipsoid.Eccentricity]);
                            if(kml_flag)
                                lat = [lat SURF(index2).lat(i_surf_index)];
                                lon = [lon SURF(index2).lon(i_surf_index)];
                                fineH = [fineH SURF(index2).h_COG(i_surf_index)];
                                stdAoA = [stdAoA SURF(index2).std_AoA(i_surf_index)];
                                peakness = [peakness SURF(index2).peakness_surf(i_surf_index)];
                                slopeAoA = [slopeAoA SURF(index2).slopeAoA(i_surf_index)];
                                slopeElev = [slopeElev SURF(index2).slopeElev(i_surf_index)];
                            end
                         
                            if(swath_processing_flag)
                                clear DEM_diff
                                for i_wrap=1: length(AoA_points)
                                    clear point_color 
                                    % locations of the section: SURF_wrap(index2).lat(i_wrap,surf_samples),SURF_wrap(index2).lon(i_wrap,surf_samples)
                                    % elevatiosn plot_Elev_wrap(i_wrap,surf_samples)
                                    lat_dem_vector = SURF_wrap(index2).lat(i_wrap,surf_samples);
                                    lon_dem_vector = SURF_wrap(index2).lon(i_wrap,surf_samples);
                                    for ipoint=1:length(lat_dem_vector)
                                        [~,lat_dem_index]= min(abs((lat_hd(:,1)-lat_dem_vector(ipoint))));
                                        [~,lon_dem_index]= min(abs((lon_hd(1,:)-lon_dem_vector(ipoint))));
                                        h_DEM(ipoint)= (DEM_hd(lat_dem_index,lon_dem_index));
                                        DEM_diff(i_wrap,ipoint) = plot_Elev_wrap(i_wrap,surf_samples(ipoint))-h_DEM(ipoint);
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
                                    if(plotting_DEM)
                                        subplot(4,3,[2 3 5 6]); 
                                        p_wrap(i_wrap,i_surf_index)=scatter3(SURF_wrap(index2).lon(i_wrap,surf_samples),SURF_wrap(index2).lat(i_wrap,surf_samples),plot_Elev_wrap(i_wrap,surf_samples),40,point_color,'filled','MarkerEdgeColor','none');
                                        subplot(4,3,8:9);
                                        scatter(SURF_wrap(index2).lon(i_wrap,surf_samples),plot_Elev_wrap(i_wrap,surf_samples),40,point_color,'filled','MarkerEdgeColor','k'); hold on;
                                        subplot(4,3,11:12);
                                        scatter(SURF_wrap(index2).lon(i_wrap,surf_samples),DEM_diff(i_wrap,:),40,point_color,'filled','MarkerEdgeColor','k'); hold on;
                                        
                                    end
                                    %                     h_DEM=double(h_DEM);
                                    %                     h_DEM(~h_DEM)= NaN;
                                    
                                end
                                [mean_DEM_error(i_surf_index),wrap_selected(i_surf_index)]=min(abs(mean(DEM_diff.')));
                                N_selected(i_surf_index)= (-wrap_selected(i_surf_index)+N_wraps+1);
                                elev_average_error=[elev_average_error; mean_DEM_error(i_surf_index)];
                                surf_index=[surf_index; i_surf_index];
                                ambiguity_wrap=[ambiguity_wrap; N_selected(i_surf_index)];
                                elevation_std=[elevation_std; std(DEM_diff(wrap_selected(i_surf_index),:))];
                                surf_samples_length=[surf_samples_length;length(surf_samples)];
                                init_sample=[init_sample;surf_samples(1)];
                                
                                SWATH.lat               =[SWATH.lat SURF_wrap(index2).lat(wrap_selected(i_surf_index),surf_samples)];
                                SWATH.lon               =[SWATH.lon SURF_wrap(index2).lon(wrap_selected(i_surf_index),surf_samples)];
                                SWATH.elev              =[SWATH.elev plot_Elev_wrap(wrap_selected(i_surf_index),surf_samples)];
                                SWATH.N_selected        = [SWATH.N_selected N_selected(i_surf_index)+zeros(1,length(surf_samples))];
                                SWATH.individual_dem_diff   =[SWATH.individual_dem_diff DEM_diff(wrap_selected(i_surf_index),:)];
                                SWATH.average_dem_diff  = [SWATH.average_dem_diff mean_DEM_error(i_surf_index)];
                                SWATH.std_dem_diff      = [SWATH.std_dem_diff std(DEM_diff(wrap_selected(i_surf_index),:))];
                                SWATH.wvf_number        = [SWATH.wvf_number index+zeros(1,length(surf_samples))];
                                SWATH.sample            = [SWATH.sample surf_samples];
                                plot(SURF_wrap(index2).lon(wrap_selected(i_surf_index),surf_samples),DEM_diff(wrap_selected(i_surf_index),:),'o','MarkerSize',8,'MarkerFaceColor',[(i_surf_index-1)/SURF(index2).N_surfs 1-(i_surf_index-1)/SURF(index2).N_surfs 0],'MarkerEdgeColor','k');
                                subplot(4,3,8:9);
                                plot(SURF_wrap(index2).lon(wrap_selected(i_surf_index),surf_samples),plot_Elev_wrap(wrap_selected(i_surf_index),surf_samples),'o','MarkerSize',8,'MarkerFaceColor',[(i_surf_index-1)/SURF(index2).N_surfs 1-(i_surf_index-1)/SURF(index2).N_surfs 0],'MarkerEdgeColor','k');
                                
                            end
                            
                            if(plotting)
                                subplot(4,3,1);hold on;
                                plot(x_axisL1(surf_samples),plot_Wvf(surf_samples),'Color',[(i_surf_index-1)/SURF(index2).N_surfs 1-(i_surf_index-1)/SURF(index2).N_surfs 0],'LineWidth',1);
                                plot(x_axisL1(SURF(index2).surf_ini_sample(i_surf_index)),plot_Wvf(SURF(index2).surf_ini_sample(i_surf_index)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
                                plot(x_axisL1(SURF(index2).surf_end_sample(i_surf_index)),plot_Wvf(SURF(index2).surf_end_sample(i_surf_index)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
%                                 plot(x_axisL1(round(SURF(index2).CoG_pos(i_surf_index))),plot_Wvf(round(SURF(index2).CoG_pos(i_surf_index))),'o','MarkerSize',6,'MarkerFaceColor','y','MarkerEdgeColor','k');
%                                 plot(x_axisL1(SURF(index2).Peak_pos(i_surf_index)),plot_Wvf(SURF(index2).Peak_pos(i_surf_index)),'o','MarkerSize',6,'MarkerFaceColor','c','MarkerEdgeColor','k');
                                
                                subplot(4,3,4);hold on;
                                plot(surf_samples,plot_Coh(surf_samples),'Color',[(i_surf_index-1)/SURF(index2).N_surfs 1-(i_surf_index-1)/SURF(index2).N_surfs 0],'LineWidth',1);
                                plot(SURF(index2).surf_ini_sample(i_surf_index),plot_Coh(SURF(index2).surf_ini_sample(i_surf_index)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
                                plot(SURF(index2).surf_end_sample(i_surf_index),plot_Coh(SURF(index2).surf_end_sample(i_surf_index)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
%                                 plot(round(SURF(index2).CoG_pos((i_surf_index))),plot_Coh(round(SURF(index2).CoG_pos(i_surf_index))),'o','MarkerSize',6,'MarkerFaceColor','y','MarkerEdgeColor','k');
%                                 plot(SURF(index2).Peak_pos(i_surf_index),plot_Coh(SURF(index2).Peak_pos(i_surf_index)),'o','MarkerSize',6,'MarkerFaceColor','c','MarkerEdgeColor','k');
                                
                                subplot(4,3,7); hold on;
                                plot(surf_samples,plot_AoA(surf_samples),'Color',[(i_surf_index-1)/SURF(index2).N_surfs 1-(i_surf_index-1)/SURF(index2).N_surfs 0],'LineWidth',1);
                                plot(SURF(index2).surf_ini_sample(i_surf_index),plot_AoA(SURF(index2).surf_ini_sample(i_surf_index)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
                                plot(SURF(index2).surf_end_sample(i_surf_index),plot_AoA(SURF(index2).surf_end_sample(i_surf_index)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k');
%                                 plot(round(SURF(index2).CoG_pos(i_surf_index)),plot_AoA(round(SURF(index2).CoG_pos(i_surf_index))),'o','MarkerSize',6,'MarkerFaceColor','y','MarkerEdgeColor','k');
                            end
                        else
                            i_surf_index = i_surf_index-1;
                        end
                    end
                end
                                
                
            else
                SURF(index2).N_surfs_OK=0;
            end
            
            %Retreived info for the SURF
            
            if(plotting)
                %subplot(3,4,3);hold on;geoshow(orbit_track); figlabels('','Ground Track [m]','',['Loc ',num2str(index)],font_size); hold off
                subplot(4,3,1); hold off;
                figlabels('Elevation [m]','Power [dB]','',['#',num2str(index,'%04d'),' L1B Waveform, #',num2str(SURF(index2).N_surfs_OK,'%d'),' Reflections retracked'] ,font_size); set(gca,'XLim',[min(x_axisL1) max(x_axisL1)],'FontSize',font_size);
                set(gca, 'xdir','reverse');
                subplot(4,3,4); hold off;
                figlabels('samples','','',['#',num2str(index,'%04d'),' Coherence'],font_size); set(gca,'XLim',[1 N_samples_sar_chd],'FontSize',font_size);set(gca,'YLim',[0 1000],'FontSize',font_size);
                subplot(4,3,7); hold off;
                figlabels('samples','AoA [deg]','',['#',num2str(index,'%04d'), ' Angle of Arrival - Roll'],font_size); set(gca,'XLim',[1 N_samples_sar_chd],'FontSize',font_size); set(gca,'YLim',[-0.6 0.6],'FontSize',font_size);
                   


            end
            if(plotting_DEM)
                AoA_points= -jump_AoA*(N_wraps+1/2):jump_AoA/2:jump_AoA*(N_wraps+1/2);% [-jump_AoA/2 0 jump_AoA/2]
                top_swath_Xdist= ((AoA_points)/180*pi-roll) * (out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12-N_samples_sar_chd/2*T0_chd/zp_factor)*c_cst/2;
                bottom_swath_Xdist  = ((AoA_points)/180*pi-roll) * (out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12+(N_samples_sar_chd/2-1)*T0_chd/zp_factor)*c_cst/2;
                top_Slant_Range    = sqrt(top_swath_Xdist.^2 + (out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12*c_cst/2)^2);
                bottom_Slant_Range    = sqrt(bottom_swath_Xdist.^2 + (out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12*c_cst/2)^2);
                top_Slant_Elev     = top_Slant_Range- out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12*c_cst/2+top_elev;
                top_Elev   = cos((AoA_points).*pi./180 -roll).*top_Slant_Elev; % need to be corrected for delta_Ha and delta_Hp: SEE TOM ARMITAGE PAPER https://www.researchgate.net/publication/260623282_Using_the_Interferometric_Capabilities_of_the_ESA_CryoSat-2_Mission_to_Improve_the_Accuracy_of_Sea_Ice_Freeboard_Retrievals
                bottom_Slant_Elev     = bottom_Slant_Range- out_L1B(1,i_record).measurements_group(1,k_record).window_delay*1e-12*c_cst/2+bottom_elev;
                bottom_Elev   = cos((AoA_points).*pi./180-roll).*bottom_Slant_Elev; % need to be corrected for delta_Ha and delta_Hp: SEE TOM ARMITAGE PAPER https://www.researchgate.net/publication/260623282_Using_the_Interferometric_Capabilities_of_the_ESA_CryoSat-2_Mission_to_Improve_the_Accuracy_of_Sea_Ice_Freeboard_Retrievals
                [lat_top_elev,lon_top_elev] = reckon(ground_track(index2,1).Lat,ground_track(index2,1).Lon,top_swath_Xdist,az_d,[ellipsoid.SemimajorAxis ellipsoid.Eccentricity]);  
                [lat_bottom_elev,lon_bottom_elev] = reckon(ground_track(index2,1).Lat,ground_track(index2,1).Lon,bottom_swath_Xdist,az_d,[ellipsoid.SemimajorAxis ellipsoid.Eccentricity]); 
                lon_bottom_elev_pairs=lon_bottom_elev(1:2:end);
                lon_top_elev_pairs=lon_top_elev(1:2:end);
                lat_bottom_elev_pairs=lat_bottom_elev(1:2:end);
                lat_top_elev_pairs=lat_top_elev(1:2:end);
                bottom_Elev_pairs=bottom_Elev(1:2:end);
                top_Elev_pairs=top_Elev(1:2:end);
                [lat_elev_DEM_north,lon_elev_DEM_north] = reckon(lat_bottom_elev_pairs,lon_bottom_elev_pairs,100e3,az_v,[ellipsoid.SemimajorAxis ellipsoid.Eccentricity]); 
                [lat_elev_DEM_south,lon_elev_DEM_south] = reckon(lat_bottom_elev_pairs,lon_bottom_elev_pairs,100e3,az_v+180,[ellipsoid.SemimajorAxis ellipsoid.Eccentricity]); 

                for i_wrap=1: length(top_Elev_pairs)
                    lon_elev_wrap(i_wrap,:)= [lon_bottom_elev_pairs(i_wrap) lon_top_elev_pairs(i_wrap)];
                    lat_elev_wrap(i_wrap,:)= [lat_bottom_elev_pairs(i_wrap) lat_top_elev_pairs(i_wrap)];
                    lon_elev_wrap_DEM(i_wrap,:)= [lon_elev_DEM_north(i_wrap) lon_elev_DEM_south(i_wrap)];
                    lat_elev_wrap_DEM(i_wrap,:)= [lat_elev_DEM_north(i_wrap) lat_elev_DEM_south(i_wrap)];
                    elev_wrap(i_wrap,:)= [0 elev_limit_plot];
                    str(i_wrap,:)=['N =' num2str(i_wrap-length(top_Elev_pairs)/2,'%+i ')];
                    elev_wrapdiff(i_wrap,:)= [-100 300];
                end
                
                str(length(top_Elev_pairs)/2,:)='NADIR';
                [lat_dem1,lon_dem1] = reckon(ground_track(index2,1).Lat,ground_track(index2,1).Lon,(8000*N_wraps*2),az_d,[ellipsoid.SemimajorAxis ellipsoid.Eccentricity]);  % N_wraps warps
                [lat_dem2,lon_dem2] = reckon(ground_track(index2,1).Lat,ground_track(index2,1).Lon,(8000*N_wraps*2),az_d+180,[ellipsoid.SemimajorAxis ellipsoid.Eccentricity]);  %N_wraps 2 wraps
                lat_res =(max(lat_dem2,lat_dem1)-(min(lat_dem2,lat_dem1)))/1000;
                lon_res =(max(max(DEM.lon),min(DEM.lon))-(min(max(DEM.lon),min(DEM.lon))))/1000;
                lat_dem_vector = lat_dem1:lat_res:lat_dem2;
                if(lon_dem1<lon_dem2)
                    lon_dem_vector = min(DEM.lon):lon_res:max(DEM.lon);
                elseif(lon_dem1>lon_dem2)
                    lon_dem_vector = min(DEM.lon):lon_res:max(DEM.lon);
                end
                for ipoint=1:length(lat_dem_vector)
                    [~,lat_dem_index]= min(abs((DEM.lat-lat_dem_vector(ipoint))));
                    [~,lon_dem_index]= min(abs((DEM.lon-lon_dem_vector(ipoint))));
                    h_DEM_strip(ipoint)= double(DEM.z(lat_dem_index,lon_dem_index));
                    
                end
                
%                 h_DEM_strip(~h_DEM_strip)= NaN;
                subplot(4,3,[2 3 5 6]); 
                plot(lon_dem_vector,lat_dem_vector,':k','LineWidth',6);
                plot(lon_elev_wrap_DEM',lat_elev_wrap_DEM','--k','LineWidth',1,'MarkerFaceColor','k','MarkerEdgeColor','k');
                long_shift=+0.025;
                if(v(1)<0) %descending
                    long_shift=-0.025;
                end
                text(lon_top_elev(end-1:-2:2)+long_shift,lat_top_elev(end-1:-2:2).'-0.2,cellstr(str(1:end-1,:)),'Color','black','FontSize',12,'HorizontalAlignment','center');

                subplot(4,3,8:9);
%                 x_range=length(lon_dem_vector(~isnan(h_DEM_strip)));
                im_area=imagesc(lon_dem_vector(~isnan(h_DEM_strip)),0:1:max(max(DEM.z)),(0:1:max(max(DEM.z))).');hold on
                set(gca,'YDir','normal');
                dem_line=area(lon_dem_vector(~isnan(h_DEM_strip)),h_DEM_strip(~isnan(h_DEM_strip)),'LineWidth',1,'BaseValue',elev_limit_plot,'FaceColor',WaterColor);
                demcmap(DEM_hd);
%                 colormap(colors_new);
                uistack(dem_line,'bottom');
                uistack(im_area,'bottom');
                plot(lon_elev_wrap',elev_wrap','--k','LineWidth',1,'MarkerFaceColor','k','MarkerEdgeColor','k');
                text(lon_top_elev(end-1:-2:2),elev_wrapdiff(1:end-1,1).'.*0+2000,cellstr(str(1:end-1,:)),'Color','white','FontSize',12,'HorizontalAlignment','center')
                set(gca,'box','on');
                set(gca,'Color',WaterColor);
                set(gca,'XLim',[min(DEM.lon) max(DEM.lon)],'FontSize',font_size);
                set(gca,'YLim',[0 elev_limit_plot],'FontSize',font_size);
                subplot(4,3,11:12);    
%                 plot(lon_bottom_elev,bottom_Elev,'--k','LineWidth',1)
%                 plot(lon_top_elev,top_Elev,'--k','LineWidth',1)
%                 plot(lon_bottom_elev(1:2:end),bottom_Elev(1:2:end),'--k','LineWidth',1,'MarkerFaceColor','r','MarkerEdgeColor','k');
%                 plot(lon_top_elev(1:2:end),top_Elev(1:2:end),'--k','LineWidth',1,'MarkerFaceColor','r','MarkerEdgeColor','k');
                plot(lon_elev_wrap',elev_wrapdiff','--k','LineWidth',1,'MarkerFaceColor','k','MarkerEdgeColor','k');
%                 set(gca,'Color',WaterColor);

                text(lon_top_elev(end-1:-2:2),elev_wrapdiff(1:end-1,1).'.*0-70,cellstr(str(1:end-1,:)),'Color','black','FontSize',12,'HorizontalAlignment','center')
                set(gca,'box','on');
            
                subplot(4,3,[2 3 5 6]); 
                if(~isempty(SWATH.lon))
                    scatter3(SWATH.lon,SWATH.lat,SWATH.elev,30,SWATH.elev,'filled','MarkerEdgeColor','k','LineWidth',0.5);
                end
                colorbar('Location','northoutside');hold off;
%                 view(-10,75);
                figlabels('Longitude [degrees]','Latitude [degrees]','Elevation [m]',['2D DEM + Unwrapping verification #',num2str(index,'%04d'), ': Green: Ns for each surf (Lowest DEM diff), Red other Ns checked', ] ,font_size); 
                set(gca,'XLim',[min(DEM.lon) max(DEM.lon)],'FontSize',font_size);
                set(gca,'YLim',[min(DEM.lat) max(DEM.lat)],'FontSize',font_size);
                set(gca,'ZLim',[0 elev_limit_plot],'FontSize',font_size);
                
                
                subplot(4,3,8:9);hold off;
                figlabels('Longitude [degrees]','Elevation [m]','',['1D DEM + Unwrapping verification #',num2str(index,'%04d'), ': Green: Ns for each surf (Lowest DEM diff), Red other Ns checked', ] ,font_size); 
                set(gca,'XLim',[min(DEM.lon) max(DEM.lon)],'FontSize',font_size);
                set(gca,'YLim',[0 elev_limit_plot],'FontSize',font_size);


                subplot(4,3,11:12);hold off;
                figlabels('Longitude [degrees]','DEM diff [m]','',['DEM diff #',num2str(index,'%04d'),'', ': Green: Ns for each surf (Lowest DEM diff), Red other Ns checked', ] ,font_size); 
                set(gca,'XLim',[min(DEM.lon) max(DEM.lon)],'FontSize',font_size);
                set(gca,'YLim',[-100 300],'FontSize',font_size);

                subp=subplot(4,3,10);
                set(subp,'Visible','off');
                dim = get(subp,'position');
                columnname =   {'init','Ns', 'Best Wrap', 'Mean Elev. diff','std Elev. diff'};
                columnformat = {'numeric','numeric' , 'numeric', 'numeric', 'numeric'};
                dat={};
                dat=horzcat((init_sample),(surf_samples_length),(ambiguity_wrap),(elev_average_error),(elevation_std));
                
                t = uitable('Units','normalized','Position',dim, 'Data', (dat),... 
                'ColumnName', columnname,...
                'ColumnFormat', columnformat,'ColumnWidth', {30,30,70,90,90},...
                'RowName',[],'FontName','Arial','Units','normalized');
%                 t_pos = get(t,'Position');
%                 t_pos(3)=0.5357142857142857;
                t_pos =[0.1300 0.0149 0.1809 0.2525];
                set(t,'Position',t_pos);
%                 t_pos(2)=0.0636;
%                 if(length(init_sample)>7)
%                     t_pos(4)=t_pos(4)+(length(init_sample)-7)*0.0188;
%                 end
%                 dim(4)=0.185;
                axis off;
                
                
%                 str_text=[[{'Elevation diff [m]',num2str(elev_average_error)}].' [{'Elevation std [m]',num2str(elevation_std)}].'].';
%                 dim=[0.4595068027210884 0.11 0.05008503401360548 0.21573529411764708];
%                 uitable('Reflexion index', surf_index,'Best N phase ambiguity' ,ambiguity_wrap, 'Mean Elevation diff wrt DEM', elev_average_error, 'std Elevation diff wrt DEM',elevation_std);
%                 annotation('textbox',dim,'String',str_text,'FitBoxToText','on','Tag' , 'results');
                
%                 set(gca,'YLim',[min_Zaxis max_Zaxis],'FontSize',font_size);

            end
            if(plotting_DEM)
%                 saveas (h, ['./results/' L1B_FILE(end-50:end-4) '/',num2str(index,'%04d'),'_DEM_and_swath.fig']);
                img = getframe(h);
                imwrite(img.cdata, ['./results/' L1B_FILE(end-50:end-4) '/',num2str(index,'%04d'),'_DEM_and_swath.png']);
%                 subplot(4,3,[2 3 5 6]); 
                
%                 delete(p_wrap);
%                 delete(dem_2D);
%                 delete(cross_track_line);
%                 delete(ambiguity_lines);
%                 delete(ambiguity_text);
%                 drawnow();
%                 hold(ax2356,'on');
%                 close(h);
%                 clear DEM h_DEM;
            end
            %         hold on; plot(AoA_pos,plot_AoA(AoA_pos),'ro','MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor','k');hold off;
            %         hold on; plot(AoA_reflex,plot_AoA(AoA_reflex),'yo','MarkerSize',4,'MarkerFaceColor','y','MarkerEdgeColor','k');hold off;legend('AoA', 'Nadir', 'diff < 0.01' )
            %
            %         hold on; plot(surf_samples,plot_Wvf(surf_samples),'g','LineWidth',3);hold off
            %
            %         subplot(4,3,4); plot(plot_Coh); figlabels('samples','','',['Coherence ',num2str(index)],18); set(gca,'XLim',[1 N_samples_sar_chd],'FontSize',18);set(gca,'YLim',[0 1000],'FontSize',18);
            %         hold on; plot(Coh_pos,plot_Coh(Coh_pos),'o','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');hold off
            %         subplot(4,3,7);plot(plot_diff_Coh,'r');
            %         hold on; plot(Coh_pos,plot_diff_Coh(Coh_pos),'o','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');hold off
            %
            %         subplot(3,4,[3 6]); plot(1:index,alt(1:index)','.k'); figlabels('samples','','',['Tracking Window Elev ',num2str(ground_track(index,1).Alt)],18);hold off;
            %         subplot(4,3,7); plot(plot_AoA); figlabels('samples','','',['Angle of Arrival - Roll ',num2str(index)],18); set(gca,'XLim',[1 N_samples_sar_chd],'FontSize',18); set(gca,'YLim',[-0.6 0.6],'FontSize',18);
            %         hold on; plot(AoA_pos,plot_AoA(AoA_pos),'ro','MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor','k');hold off;
            %         hold on; plot(AoA_reflex,plot_AoA(AoA_reflex),'yo','MarkerSize',4,'MarkerFaceColor','y','MarkerEdgeColor','k');hold off;legend('AoA', 'Nadir', 'diff < 0.01' )
            
            %         saveas (h, ['./',num2str(index,'%04d'),'_Coherence_greenland.fig']);
        end
        
    end
    if(kml_flag)
        lla2kml(['./results/' L1B_FILE(end-50:end-4) '_L2_Multisurf.kml'],lat,lon,fineH,'.');
    end
    save(['./results/L2I/' L1B_FILE(end-50:end-4) '.mat'],'SURF');
    fclose all;
    %     movefile(L1B_FILE,['./inputs_' inputyear '/done/']);
%     if(ambiguity_flag)
%         save(['./results/' L1B_FILE(30:end-4) '_ambiguous.mat'],'SURF_left', 'SURF_right')
%     end
    if(plotting)
        close(h);
    end
    time(i_file) = toc(processing_time);
    minutes= floor(time(i_file)/60);
    secs = time(i_file)-minutes*60;
    %strcat('Processing L1B: ', L1B_FILE(29:end-4));
    disp(['Processing L1B: ' L1B_FILE(end-50:end-4) ' ' num2str(i_file) '/' num2str(N_files)  ' in: ' num2str(minutes) ' minutes, ' num2str(secs) ' seconds; ' num2str(sum(time)/60) ' minutes passed'])
    
end
end
