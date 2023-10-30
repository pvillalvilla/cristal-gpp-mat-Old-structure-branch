%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: Computation of the surface locations over which the Fully
% Focused Processing in time domain will be performed
%
% Calling: 
%
% INPUTs:
%  time_sar_ku_pulse:       time/datation of each pulse
%  alt_rate_sar_sat_pulse:  height rate of the satellite per pulse
%  x/y/z_vel_sat_sar_pulse: x/y/z satellite velocity comp. per pulse 
%  roll_sar_pulse:          roll of the satellite per pulse
%  pitch_sar_pulse:         pitch of the satellite per pulse
%  yaw_sar_pulse:           yaw of the satellite per pulse
%  x/y/z_sar_sat_pulse:     x/y/z positions of the satellite per pulse 
%  lat_sar_sat_pulse:       latitude of satellite per pulse 
%  lon_sar_sat_pulse:       longitude of satellite per pulse
%  alt_sar_sat_pulse:       altitude of the satellite per pulse
%  win_delay_sar_ku_pulse:  window delay of each pulse
%  pri_sar_pre_dat_pulse:   pulse repetition interval per pulse
%  T0_sar_pre_dat_pulse:    Sampling period per pulse
%  cnf:                     Configuration parameters structure
%  chd:                     Characterization parameters structure
%  cst:                     Constant parameters structure
%
% OUTPUTs:
%
%  x/y/z_surf:             x/y/z positions of surfaces
%  lat_surf:               latitude of surfaces
%  lon_surf:               longitude of surfaces
%  alt_surf:               altitude of surfaces
%  win_delay_surf:         window delay associated to surfaces
%  time_surf:              datation time of surfaces
%  x/y/z_vel_sat_surf:     x/y/z satellite velocity over the surfaces
%  x/y/z_sat_surf:         x/y/z positions of satellite over the surfaces
%  alt_rate_sat_surf:      height rate of satellite over surfaces
%  alt_sat_surf:           altitude of the satellite over surfaces
%  pitch_sat_surf:         pitch of satellite over surfaces
%  roll_sat_surf:          roll of satellite over surfaces
%  yaw_sat_surf:           yaw of satellite over surfaces        
%  TRP_surf_idx:           index of the surface corresponding to TRP
%
%
% COMMENTS: Current version is based on projecting the desired spatial
% separation of the surfaces into a temporal step used to construct the
% datation of the surface from the first pulse datation and for all the
% surfaces within the temporal span of the acquisition: then related info
% is interpolated from the pulse information with the given time of the
% surface and the geo-locations computed. For transponder case the
% along-track surfaces are centered around the transponder and the total
% number of output surfaces is decided based on a spatial distance from
% transponder per configuration
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%            Albert Garcia / isardSAT
%            Ferran Gibert / isardSAT
% ----------------------------------------------------------
% v1.1  Simplified version
% 		-Forced case 0 and removed case 1 (BI-DIMENSIONAL GRID FOR A GIVEN SPACING ALONG-TRACK BASED ON TIME)
% v1.2 2019/06/20 - Enable possibility of forcing processing over desired and 
%                  coordinates (cnf.FFt.force_POI)
%				  - Added initialisation of variable 'TRP_surf_idx'
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [x_surf, y_surf, z_surf,...
        lat_surf, lon_surf, alt_surf, win_delay_surf,...
        time_surf,...
        x_vel_sat_surf,y_vel_sat_surf,z_vel_sat_surf,...
        x_sat_surf,y_sat_surf,z_sat_surf,...
        alt_rate_sat_surf,alt_sat_surf,...
        pitch_sat_surf,roll_sat_surf,yaw_sat_surf,...        
        TRP_surf_idx]=...
        focusing_grid_computation(time_sar_ku_pulse,...
        alt_rate_sar_sat_pulse, x_vel_sat_sar_pulse, y_vel_sat_sar_pulse, z_vel_sat_sar_pulse,...        
       roll_sar_pulse, pitch_sar_pulse, yaw_sar_pulse,...
       x_sar_sat_pulse, y_sar_sat_pulse, z_sar_sat_pulse,...       
       lat_sar_sat_pulse,lon_sar_sat_pulse,alt_sar_sat_pulse,...
       win_delay_sar_ku_pulse,...    
       pri_sar_pre_dat_pulse,T0_sar_pre_dat_pulse,...
       cnf,chd,cst)
   
   %forced until any new potential approach is available
   cnf.FFt.grid_computation_option = 0;
   %% -------- Computation of output surfaces where to focus --------------
   TRP_surf_idx = 0;
   switch cnf.FFt.grid_computation_option
       
       case 0
           %% FOR A GIVEN SPACING ALONG-TRACK BASED ON TIME
           % convert from space to time using the ground velocity of the
           % first pulse (approximation) timing same for all the surfaces
           % but separation probably not the same 
           
           %ground velocity
           R_Earth =sqrt((cst.semi_major_axis^4*cos(lat_sar_sat_pulse(1).*pi/180).^2+cst.semi_minor_axis^4*sin(lat_sar_sat_pulse(1)*pi/180).^2)./...
               ((cst.semi_major_axis*cos(lat_sar_sat_pulse(1).*pi/180)).^2+(cst.semi_minor_axis*sin(lat_sar_sat_pulse(1).*pi/180)).^2));
           alpha_r = 1 + alt_sar_sat_pulse(1)/R_Earth; % we should modify earth curvature with lat and longitude
           vg = sqrt(x_vel_sat_sar_pulse(1).^2+y_vel_sat_sar_pulse(1).^2+...
                    z_vel_sat_sar_pulse(1).^2)/alpha_r;
           
           %delta time to locate the surfaces along the track
           delta_time = cnf.FFt.spacing_al_surfaces/vg;
           
           N_surfaces = ceil((time_sar_ku_pulse(end)-time_sar_ku_pulse(1))/delta_time);
           
           %time of surfaces of interest along-track
           time_surf = time_sar_ku_pulse(1)+(0:1:N_surfaces).*delta_time;
           
           %window delay
           win_delay_surf = spline(time_sar_ku_pulse,win_delay_sar_ku_pulse,time_surf);
           
           
           % LLA
           lat_surf = spline(time_sar_ku_pulse,lat_sar_sat_pulse,time_surf);
           lon_surf = spline(time_sar_ku_pulse,lon_sar_sat_pulse,time_surf);
           alt_surf = spline(time_sar_ku_pulse,alt_sar_sat_pulse,time_surf)-win_delay_surf*cst.c/2;
           
           %ECEF
           xyz_surf = lla2ecef([lat_surf.',lon_surf.',alt_surf.'],cst.flat_coeff,cst.semi_major_axis);
           x_surf = xyz_surf(:,1).';
           y_surf = xyz_surf(:,2).';
           z_surf = xyz_surf(:,3).';
           clear xyz_surf
           
           %SATELLITE related info
           x_vel_sat_surf = spline(time_sar_ku_pulse,x_vel_sat_sar_pulse,time_surf);
           y_vel_sat_surf = spline(time_sar_ku_pulse,y_vel_sat_sar_pulse,time_surf);
           z_vel_sat_surf = spline(time_sar_ku_pulse,z_vel_sat_sar_pulse,time_surf);
           
           alt_rate_sat_surf = spline(time_sar_ku_pulse,alt_rate_sar_sat_pulse,time_surf);
           alt_sat_surf = spline(time_sar_ku_pulse,alt_sar_sat_pulse,time_surf);
           
           %xyz_sat_surf = lla2ecef([lat_surf.',lon_surf.',alt_sat_surf.'],cst.flat_coeff,cst.semi_major_axis);
           x_sat_surf = spline(time_sar_ku_pulse,x_sar_sat_pulse,time_surf);
           y_sat_surf = spline(time_sar_ku_pulse,y_sar_sat_pulse,time_surf);
           z_sat_surf = spline(time_sar_ku_pulse,z_sar_sat_pulse,time_surf);
            
           
           pitch_sat_surf = spline(time_sar_ku_pulse,pitch_sar_pulse,time_surf);
           roll_sat_surf  = spline(time_sar_ku_pulse,roll_sar_pulse,time_surf);
           yaw_sat_surf   = spline(time_sar_ku_pulse,yaw_sar_pulse,time_surf);
           
           
           if cnf.trp_flag==1 %JPLZ: deleted the '|| cnf.FFt.force_POI==1', don't know why it was here
               % TRANSPONDER MODE
               % TO IMAGE 2D IRF MOVE THE WHOLE TRACK OF PULSES CENTERED AT POSITION OF TRANSPONDER
               xyz_TRP = [chd.x_trp, chd.y_trp, chd.z_trp];
               xyz_sar_surf = [x_surf.',y_surf.',z_surf.'];
               
               Distances_TRP=sqrt(sum((xyz_sar_surf -...
                                    (ones(length(xyz_sar_surf(:,1)),1)*xyz_TRP)).^2,2));
               [TRP_surf_dist,TRP_surf_idx] = min(Distances_TRP);
               
               % ---- Move projected surface to center it to TRP location -
               xyz_over_TRP = (xyz_sar_surf-...
                   (ones(length(xyz_sar_surf(:,1)),1)*xyz_sar_surf(TRP_surf_idx(1),:)))...
                   + (ones(length(xyz_sar_surf(:,1)),1)*xyz_TRP);
               
               
               %>> Cartesian to geodetic
               latlon_sar_surf = ecef2lla([xyz_over_TRP(:,1),xyz_over_TRP(:,2),xyz_over_TRP(:,3)],cst.flat_coeff,cst.semi_major_axis);%
               
%                               %----------- New datation time of the displaced surfaces----
%                %As original surfaces moved around the location of TRP: the related
%                %datation shall be re-computed
%                %assume the latitude information can be exploited for
%                %interpolation of the time: lat_surf original locations
%                %nadir along-track ??
%                time_surf_nonmov = time_surf; %time of the original surfacs non-moved or dispalced to be centered around TRP
%                time_surf = spline(lat_surf,time_surf_nonmov,latlon_sar_surf(:,1).');
%                
%                
%                % Compute the output information (satellite positions & T0/PRI at the new datation time of the dispalced surfaces around TRP)
%                %SATELLITE related info
%                x_vel_sat_surf = spline(time_sar_ku_pulse,x_vel_sat_sar_pulse,time_surf);
%                y_vel_sat_surf = spline(time_sar_ku_pulse,y_vel_sat_sar_pulse,time_surf);
%                z_vel_sat_surf = spline(time_sar_ku_pulse,z_vel_sat_sar_pulse,time_surf);
%                
%                alt_rate_sat_surf = spline(time_sar_ku_pulse,alt_rate_sar_sat_pulse,time_surf);
%                alt_sat_surf = spline(time_sar_ku_pulse,alt_sar_sat_pulse,time_surf);
%                
%                %xyz_sat_surf = lla2ecef([lat_surf.',lon_surf.',alt_sat_surf.'],flat_coeff_cst,semi_major_axis_cst);
%                x_sat_surf = spline(time_sar_ku_pulse,x_sar_sat_pulse,time_surf);
%                y_sat_surf = spline(time_sar_ku_pulse,y_sar_sat_pulse,time_surf);
%                z_sat_surf = spline(time_sar_ku_pulse,z_sar_sat_pulse,time_surf);
%                
%                
%                pitch_sat_surf = spline(time_sar_ku_pulse,pitch_sar_pulse,time_surf);
%                roll_sat_surf  = spline(time_sar_ku_pulse,roll_sar_pulse,time_surf);
%                yaw_sat_surf   = spline(time_sar_ku_pulse,yaw_sar_pulse,time_surf);
%                
%                %compute T0_surf and pri_surf for each surface
%                T0_surf        = spline(time_sar_ku_pulse,T0_sar_pre_dat_pulse,time_surf);
%                pri_surf       = spline(time_sar_ku_pulse,pri_sar_pre_dat_pulse,time_surf);
%                %------------------------------------------------------------
%                
%                %Define the lla of the moved surfaces                             
               lat_surf = latlon_sar_surf(:,1).';
               lon_surf = latlon_sar_surf(:,2).';
               %force all surfaces to have the same height as TRP
               %to analyze the IRF of the TRP (cut in the azimuth-along-track dimension)
               alt_surf = ones(1,length(lat_surf))*chd.alt_trp;
               %alt_surf = latlon_sar_surf(:,3).';
               %alt_surf(TRP_surf_idx) = chd.alt_trp;
               
                              %>> Geodetic to cartesian
               xyz_over_TRP = lla2ecef([lat_surf.', ...
                                        lon_surf.', ...
                                        alt_surf.'], ...
                                        cst.flat_coeff,cst.semi_major_axis);%
               x_surf = xyz_over_TRP(:,1).';
               y_surf = xyz_over_TRP(:,2).';
               z_surf = xyz_over_TRP(:,3).';
               
                             
%               wind delay
                win_delay_surf = (alt_sat_surf-alt_surf)*2/cst.c;
               
%                win_delay_surf =sqrt((x_sat_surf-x_surf).^2+...
%                                     (y_sat_surf-y_surf).^2+...
%                                     (z_sat_surf-z_surf).^2)*2/cst.c;
               

               %consider only surfaces with given distance around transponder
               %arclength distance
               [arclen,~]=distance(lat_surf, lon_surf, ...
                                    chd.lat_trp*ones(1,length(lat_surf)), ...
                                    chd.lon_trp*ones(1,length(lat_surf)), ...
                                    [cst.semi_major_axis sqrt(1-(cst.semi_minor_axis/cst.semi_major_axis).^2)]);
               

               idx_surf_int = find(abs(arclen)<=cnf.FFt.surf_dist_trp);
                                         
               TRP_surf_idx = find(lat_surf(idx_surf_int)==chd.lat_trp);
               
               % Window delay
               win_delay_surf = win_delay_surf(idx_surf_int);
               
               % Time
               time_surf = time_surf(idx_surf_int);
               
               %ECEF
               x_surf = x_surf(idx_surf_int);
               y_surf = y_surf(idx_surf_int);
               z_surf = z_surf(idx_surf_int);
               
               % LLA
               lat_surf = lat_surf(idx_surf_int);
               lon_surf = lon_surf(idx_surf_int);
               alt_surf = alt_surf(idx_surf_int);
               
               %SATELLITE related info
               x_vel_sat_surf = x_vel_sat_surf(idx_surf_int);
               y_vel_sat_surf = y_vel_sat_surf(idx_surf_int);
               z_vel_sat_surf = z_vel_sat_surf(idx_surf_int);
               
               x_sat_surf     = x_sat_surf(idx_surf_int);
               y_sat_surf     = y_sat_surf(idx_surf_int);
               z_sat_surf     = z_sat_surf(idx_surf_int);
               
               alt_rate_sat_surf = alt_rate_sat_surf(idx_surf_int);
               alt_sat_surf = alt_sat_surf(idx_surf_int);
               
               pitch_sat_surf = pitch_sat_surf(idx_surf_int);
               roll_sat_surf  = roll_sat_surf(idx_surf_int);            
               yaw_sat_surf   = yaw_sat_surf(idx_surf_int);  
               
%                T0_surf   = T0_surf(idx_surf_int);
%                pri_surf   = pri_surf(idx_surf_int);
               
           end
         
   end

end