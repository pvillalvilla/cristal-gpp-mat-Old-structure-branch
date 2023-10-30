%% HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
% ---------------------------------------------------------
% Objective: Compute the surface locations defined by the intersection 
% of the Doppler beams and the estimated surface positions along the 
% satellite track. 
% ----------------------------------------------------------
% Author:    Albert Garcia  / isardSAT
%            Eduard Makhoul / isardSAT
% ----------------------------------------------------------
% Version  record
% 1.0 2018/07/18 First version imported from Dedop rev 125
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [L1BS, i_surf, cnf] = surface_locations (L1A,L1BS,N_total_bursts, i_burst, i_surf, cnf, chd, cst)

%NOTE i_surf is refered to the next surface location. If i_surf = 4, we
%have already computed 3 locations and we are looking for the 4th.

%% 1.BURST POSITION ON GROUND done in reading routine     


if(i_burst==1)
% 1. Initialisation: the first surface location and all the other associated values are initialised:
    
    L1BS(i_surf).lat_surf    = L1A(i_burst).lat_sar_sat;
    L1BS(i_surf).lon_surf    = L1A(i_burst).lon_sar_sat;
    L1BS(i_surf).alt_surf    = L1A(i_burst).alt_sar_surf;
    L1BS(i_surf).x_surf      = L1A(i_burst).x_sar_surf;
    L1BS(i_surf).y_surf      = L1A(i_burst).y_sar_surf;
    L1BS(i_surf).z_surf      = L1A(i_burst).z_sar_surf;
    L1BS(i_surf).alt_rate_sat= L1A(i_burst).alt_rate_sar_sat;
    L1BS(i_surf).time_surf   = L1A(i_burst).time;
    L1BS(i_surf).lat_sat     = L1A(i_burst).lat_sar_sat;
    L1BS(i_surf).lon_sat     = L1A(i_burst).lon_sar_sat;
    L1BS(i_surf).alt_sat    	= L1A(i_burst).alt_sar_sat;
    L1BS(i_surf).x_sat       = L1A(i_burst).x_sar_sat ;
    L1BS(i_surf).y_sat       = L1A(i_burst).y_sar_sat ;
    L1BS(i_surf).z_sat       = L1A(i_burst).z_sar_sat ;
    L1BS(i_surf).x_vel_sat   = L1A(i_burst).x_vel_sat_sar;
    L1BS(i_surf).y_vel_sat   = L1A(i_burst).y_vel_sat_sar;
    L1BS(i_surf).z_vel_sat 	= L1A(i_burst).z_vel_sat_sar;
    L1BS(i_surf).pitch_surf 	= L1A(i_burst).pitch_sar;
    L1BS(i_surf).roll_surf  	= L1A(i_burst).roll_sar;
    L1BS(i_surf).yaw_surf    = L1A(i_burst).yaw_sar;
    L1BS(i_surf).win_delay_surf  = L1A(i_burst).win_delay;
    L1BS(i_surf).surf_counter    = i_surf;
    
    i_surf = i_surf+1;
else
    

%% 2. COMPUTE ANGULAR RESOLUTION   

vel_sat = [L1BS(i_surf-1).x_vel_sat, L1BS(i_surf-1).y_vel_sat, L1BS(i_surf-1).z_vel_sat];
ang_az_beam_res = asin (cst.c / (2 * chd.freq * norm(vel_sat) * (chd.burst_duration_sar(1))))/cnf.zp_fact_azimut;

% Tcycle = 4.66785e-02;
% Tcycle = (N_bursts_cycle_chd * bri_nom);
% (with another orbit, a Tcycle was enough but, as a new surface locations
% is not seen until 9 bursts, this Tcycle now means another thing.
% For safety, 2 "Tcycle" are applied.

% while (L1BS.time_surf(i_surf-1) <= (max (L1A.time) - Tcycle))
    
%% 3. Coarse intersection loop
    
    a = 0;
    i_curr = i_burst-1;
    
%     while a < ang_az_beam_res(i_surf-1) && (i_curr+i_loop) < L1A.N_total_bursts
        % Compute the angle a between the following vectors:
        v = [L1BS(i_surf - 1).x_surf, L1BS(i_surf - 1).y_surf, L1BS(i_surf - 1).z_surf] - [L1BS(i_surf - 1).x_sat, L1BS(i_surf - 1).y_sat, L1BS(i_surf - 1).z_sat];
        w = [L1A(i_burst).x_sar_surf, L1A(i_burst).y_sar_surf, L1A(i_burst).z_sar_surf] - [L1BS(i_surf - 1).x_sat, L1BS(i_surf - 1).y_sat, L1BS(i_surf - 1).z_sat];
        a = real(acos ((v*w') / (norm(v) * norm(w))));
%         i_loop = i_loop + 1;
%     end
if (a > ang_az_beam_res)
    
    L1BS(i_surf).surf_counter    = i_surf;

    N_bursts_interpol = 8;
    end_burst = length(L1A);
    start_burst = max(1,end_burst-N_bursts_interpol);
  
    method = 0;
    % 5. Fine intersection loop
    
    if method == 0 %CRYOSAT
        % A ratio between angles B1 and B2 is made
        % 'a' is B2
        B2 = a;
        % B1 has to be calculated (angle between v and w_aux)
        w_aux = [L1A(i_curr).x_sar_surf, L1A(i_curr).y_sar_surf, L1A(i_curr).z_sar_surf] - [L1BS(i_surf-1).x_sat, L1BS(i_surf-1).y_sat, L1BS(i_surf-1).z_sat];
        B1 = real(acos (dot(v,w_aux) / (norm(v) * norm(w_aux))));
        
        alfa = (ang_az_beam_res - B1) / (B2 - B1);
        
        
        
    elseif method == 1 %SENTINEL
        
%         % In order to precisely determine the intersection between the angular Doppler resolution direction and the surface, an interpolation over the surface is performed.
%         % Divide the interval where the interpolation has to be performed [time(i_curr), time(i_curr + 1)] in cnf.N_surf_samp_interpol equally spaced parts:
        for i_interp = 1:cnf.N_points_spline_surf
            time_interp (i_interp) = L1A(i_curr).time + (i_interp+1) * (L1A(i_curr+1).time  - L1A(i_curr).time) / cnf.N_points_spline_surf;
        end
        
%         % Then perform a cubic splines interpolation with a smoothing factor of cnf.smooth_fact_surf_pos over the surface taking cnf.N_points_spline_surf/2 samples
%         % from each side of the current position [x_surf_geoloc(i_curr), y_surf_geoloc(i_curr), z_surf_geoloc(i_curr)] if possible. If not, take the maximum locations available within
%         % the limits. The resulting interpolated surface is named [x_surf_interp, y_surf_interp, z_surf_interp].
%         
        start_surf = max (1, i_curr - cnf.N_points_spline_surf/2 + 1);
        end_surf = min (N_total_bursts, i_curr + cnf.N_points_spline_surf/2);
%         
%         %**** Tractament de n?mero de punts imparell ****%
%         if mod (end_surf - start_surf + 1, 2) ~= 0
%             if start_surf == 1
%                 end_surf = i_curr + cnf.N_points_spline_surf/2 + 1;
%             end
%             if end_surf == L1A.N_total_bursts
%                 start_surf = i_curr - cnf.N_points_spline_surf/2;
%             end
%         end
%         
        for i_interp = 1:cnf.N_points_spline_surf
            x_surf_interp(i_interp) = spline (L1A(start_surf:end_surf).time, L1A(start_surf:end_surf).x_sar_surf, time_interp(i_interp));
            y_surf_interp(i_interp) = spline (L1A(start_surf:end_surf).time, L1A(start_surf:end_surf).y_sar_surf, time_interp(i_interp));
            z_surf_interp(i_interp) = spline (L1A(start_surf:end_surf).time, L1A(start_surf:end_surf).z_sar_surf, time_interp(i_interp));
        end
%         
%         
%         % Now, compute the fine intersection following the same procedure as in the coarse intersection loop.
         i_loop_fine = 1;
         g = 0;
        while g < ang_az_beam_res(i_surf - 1)
        % Compute the angle g between the following vectors:
            v = [L1BS(i_surf - 1).x_surf, L1BS(i_surf - 1).y_surf, L1BS(i_surf - 1).z_surf] - [L1BS(i_surf - 1).x_sat, L1BS(i_surf - 1).y_sat, L1BS(i_surf - 1).z_sat];
            w = [x_surf_interp(i_loop_fine), y_surf_interp(i_loop_fine), z_surf_interp(i_loop_fine)] - [L1BS(i_surf - 1).x_sat, L1BS(i_surf - 1).y_sat, L1BS(i_surf - 1).z_sat];
            g = real(acos (dot(v,w) / (norm(v) * norm(w))));
            i_loop_fine = i_loop_fine + 1;
        end
        i_loop_fine = i_loop_fine - 1;
	
    end
    
%%    
    % 6. Determination of the new surface location
    if method == 0
        L1BS(i_surf).x_surf = L1A(i_curr).x_sar_surf + alfa * (L1A(i_curr+1).x_sar_surf - L1A(i_curr).x_sar_surf);
        L1BS(i_surf).y_surf = L1A(i_curr).y_sar_surf + alfa * (L1A(i_curr+1).y_sar_surf - L1A(i_curr).y_sar_surf);
        L1BS(i_surf).z_surf = L1A(i_curr).z_sar_surf + alfa * (L1A(i_curr+1).z_sar_surf - L1A(i_curr).z_sar_surf);
    elseif method == 1
        L1BS(i_surf).x_surf = x_surf_interp(i_loop_fine);
        L1BS(i_surf).y_surf = y_surf_interp(i_loop_fine);
        L1BS(i_surf).z_surf= z_surf_interp(i_loop_fine);
    end
    
    % ==> Then convert them to geodetic coordinates (see section 5.3.3) and get lat_surf(i_surf), lon_surf(i_surf) and alt_surf(i_surf).
%     [lat_surf(i_surf), lon_surf(i_surf), alt_surf(i_surf)] = gen_mec_con_06 (x_surf(i_surf), y_surf(i_surf), z_surf(i_surf), cst.semi_major_axis, cst.flat_coeff, cnf.accu_orbit_alt, cnf.accu_orbit_lat);
    lla = ecef2lla([L1BS(i_surf).x_surf,L1BS(i_surf).y_surf,L1BS(i_surf).z_surf],cst.flat_coeff,cst.semi_major_axis);
    L1BS(i_surf).lat_surf = lla(1);
    L1BS(i_surf).lon_surf = lla(2);
    L1BS(i_surf).alt_surf = lla(3);
    
%% 
    % 7. Determination of Orbit State and associated Window delay
    % The corresponding time to the new surface location is:
    if method == 0
        L1BS(i_surf).time_surf = L1A(i_curr).time + alfa * (L1A(i_curr+1).time - L1A(i_curr).time);
    elseif method == 1
        L1BS(i_surf).time_surf = time_interp(i_loop_fine);
    end
    
    
    % FER AMB EL GEOLOCATION (al codi en C++, no aqu?)!!! (8-point Lagrange interpolation)
    % ==> In order to get the satellite position and velocity at this time, perform an interpolation (the same type of interpolation used in step 5, now using cnf.N_points_spline_orbit
    % instead of cnf.N_points_spline_surf) for each parameter and get the value at time_surf(i_surf). Note that the parameters to be interpolated in this step are
    % [x_sar_sat, y_sar_sat, z_sar_sat] and [x_vel_sat_sar, y_vel_sat_sar, z_vel_sat_sar] respectively (using the corresponding smoothing factors). The result of
    % the interpolation will be [x_sat(i_surf), y_sat(i_surf), z_sat(i_surf)] and [x_vel_sat(i_surf), y_vel_sat(i_surf), z_vel_sat(i_surf)].
    
    
    % -------------------------------------------------------
%     % prova
%     x_surf2(i_surf) = spline (time(start_surf:end_surf), x_sar_surf(start_surf:end_surf), time_surf(i_surf));
%     y_surf2(i_surf) = spline (time(start_surf:end_surf), y_sar_surf(start_surf:end_surf), time_surf(i_surf));
%     z_surf2(i_surf) = spline (time(start_surf:end_surf), z_sar_surf(start_surf:end_surf), time_surf(i_surf));
%     [lat_surf2(i_surf), lon_surf2(i_surf), alt_surf2(i_surf)] = gen_mec_con_06 (x_surf2(i_surf), y_surf2(i_surf), z_surf2(i_surf), cst.semi_major_axis, cst.flat_coeff, cnf.accu_orbit_alt, cnf.accu_orbit_lat);
    %--------------------------------------------------------
    
    L1BS(i_surf).x_sat = spline ([L1A(start_burst:end_burst).time], [L1A(start_burst:end_burst).x_sar_sat], L1BS(i_surf).time_surf);
    L1BS(i_surf).y_sat = spline ([L1A(start_burst:end_burst).time], [L1A(start_burst:end_burst).y_sar_sat], L1BS(i_surf).time_surf);
    L1BS(i_surf).z_sat = spline ([L1A(start_burst:end_burst).time], [L1A(start_burst:end_burst).z_sar_sat], L1BS(i_surf).time_surf);
    
    lla = ecef2lla([L1BS(i_surf).x_sat,L1BS(i_surf).y_sat,L1BS(i_surf).z_sat],cst.flat_coeff,cst.semi_major_axis);
    L1BS(i_surf).lat_sat = lla(1);
    L1BS(i_surf).lon_sat = lla(2);
    L1BS(i_surf).alt_sat = lla(3);
    
    L1BS(i_surf).x_vel_sat      = spline ([L1A(start_burst:end_burst).time], [L1A(start_burst:end_burst).x_vel_sat_sar]   , L1BS(i_surf).time_surf);
    L1BS(i_surf).y_vel_sat      = spline ([L1A(start_burst:end_burst).time], [L1A(start_burst:end_burst).y_vel_sat_sar]   , L1BS(i_surf).time_surf);
    L1BS(i_surf).z_vel_sat      = spline ([L1A(start_burst:end_burst).time], [L1A(start_burst:end_burst).z_vel_sat_sar]   , L1BS(i_surf).time_surf);
    L1BS(i_surf).alt_rate_sat   = spline ([L1A(start_burst:end_burst).time], [L1A(start_burst:end_burst).alt_rate_sar_sat], L1BS(i_surf).time_surf);
    
    
    % ==> In addition, perform an interpolation of pitch_sar, roll_sar and yaw_sar in order to obtain the pitch_surf(i_surf), roll_surf(i_surf) and yaw_surf(i_surf).
    L1BS(i_surf).pitch_surf = spline ([L1A(start_burst:end_burst).time], [L1A(start_burst:end_burst).pitch_sar]  , L1BS(i_surf).time_surf);
    L1BS(i_surf).roll_surf  = spline ([L1A(start_burst:end_burst).time], [L1A(start_burst:end_burst).roll_sar]   , L1BS(i_surf).time_surf);
    L1BS(i_surf).yaw_surf   = spline ([L1A(start_burst:end_burst).time], [L1A(start_burst:end_burst).yaw_sar]    , L1BS(i_surf).time_surf);
    
    
    % On the other hand, the associated window delay is calculated as follows:
    L1BS(i_surf).win_delay_surf = (L1BS(i_surf).alt_sat - L1BS(i_surf).alt_surf) * 2 / cst.c;
    i_surf = i_surf + 1;
    




end


end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Si es canvia la posicio d'una surf loc, tb s'haurien de canviar les posicions corresponents del sat.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (cnf.trp_flag== 1&& i_surf > 2) %Transponder Flag
    

    
%     chd.lat_trp = [5.15026158104065033, 5.16332102359257661, 5.17635262542077754, 5.18935609932659219] * 10;
%     chd.lon_trp = [-1.08806662380434943, -1.08640334999590337, -1.08472972811216820, -1.08304567190662297] * 100;
%     chd.alt_trp = [93.2553967992012076, 139.666851119641251, 185.931829473920970, 232.051149092048945];
    
    if cnf.trp_flag_method == 1 %computing the surface locations to be replaced (choosing the nearest ones)
%         for i_pt = 1:length(chd.alt_trp)
           
        p = lla2ecef([chd.lat_trp.', chd.lon_trp.', chd.alt_trp.'], cst.flat_coeff, cst.semi_major_axis);
        x_trp = p(1);
        y_trp = p(2);
        z_trp = p(3);

		L1BS(i_surf-1).dist = norm([L1BS(i_surf-1).x_surf,L1BS(i_surf-1).y_surf,L1BS(i_surf-1).z_surf] - [x_trp,y_trp,z_trp]);
       
%             [~,index] = min(dist);
%             disp(index);
        if(L1BS(i_surf-1).dist-L1BS(i_surf-2).dist>0)
%             L1BS(i_surf-2).dist = 0; 
            L1BS(i_surf-2).surface_type_flag = 4;%use it as flag
            
            u = [L1BS(i_surf-2).x_surf,L1BS(i_surf-2).y_surf,L1BS(i_surf-2).z_surf] - [L1BS(i_surf-3).x_surf,L1BS(i_surf-3).y_surf,L1BS(i_surf-3).z_surf]; 
            u_norm = u / norm(u);
            w = [x_trp,y_trp,z_trp] - [L1BS(i_surf-3).x_surf,L1BS(i_surf-3).y_surf,L1BS(i_surf-3).z_surf];
            cos_a = (u_norm * w.')/norm(w);  % norm(u_norm) = 1
            v = (norm(w) * cos_a) * u_norm;
            % For normal processing
            L1BS(i_surf-2).x_surf = L1BS(i_surf-3).x_surf + v(1);
            L1BS(i_surf-2).y_surf = L1BS(i_surf-3).y_surf + v(2);
            L1BS(i_surf-2).z_surf = L1BS(i_surf-3).z_surf + v(3);
            L1BS(i_surf-2).win_delay_surf = norm([L1BS(i_surf-1).x_sat,L1BS(i_surf-1).y_sat,L1BS(i_surf-1).z_sat] - [x_trp,y_trp,z_trp])/cst.c*2;
            %Recomputing the window delay in order to have the wvf in
            %middle of window
            
            lla = ecef2lla([L1BS(i_surf-2).x_surf,L1BS(i_surf-2).y_surf,L1BS(i_surf-2).z_surf],cst.flat_coeff,cst.semi_major_axis);
            L1BS(i_surf-2).lat_surf = lla(1);
            L1BS(i_surf-2).lon_surf = lla(2);
            L1BS(i_surf-2).alt_surf = L1BS(i_surf-2).alt_sat-L1BS(i_surf-2).win_delay_surf*cst.c/2;

            p = lla2ecef([L1BS(i_surf-2).lat_surf,L1BS(i_surf-2).lon_surf,L1BS(i_surf-2).alt_surf],cst.flat_coeff,cst.semi_major_axis);
            L1BS(i_surf-2).x_surf = p(1);

            L1BS(i_surf-2).y_surf = p(2);

            L1BS(i_surf-2).z_surf = p(3);
            L1BS(i_surf-2).dist = norm([L1BS(i_surf-2).x_surf,L1BS(i_surf-2).y_surf,L1BS(i_surf-2).z_surf] - [x_trp,y_trp,z_trp]);

%             L1BS(i_surf-2).win_delay_surf = (L1BS(i_surf-2).alt_sat - L1BS(i_surf-2).alt_surf) * 2 / cst.c;
			L1BS(i_surf-2).time_surf = interp1 ([L1BS(i_surf-3).lat_surf L1BS(i_surf-1).lat_surf], [L1BS(i_surf-3).time_surf L1BS(i_surf-1).time_surf], L1BS(i_surf-2).lat_surf,'linear');
            cnf.trp_flag = -1;
            
        end
 
    elseif cnf.trp_flag_method == 2 %forcing the surface locations to be replaced (GPP Ph3 Scenario #1)
        pt_index = [30,90,150,210];
        
        for i_pt = 1:length(chd.alt_trp)
            p = lla2ecef([chd.lat_trp(i_pt).', chd.lon_trp(i_pt).', chd.alt_trp(i_pt).'], cst.flat_coeff, cst.semi_major_axis);
            x_trp(i_pt) = p(1);
            y_trp(i_pt) = p(2);
            z_trp(i_pt) = p(3);
            
%             lat_surf(pt_index(i_pt)) = chd.lat_trp(i_pt);
%             L1BS.lon_surf(pt_index(i_pt)) = chd.lon_trp(i_pt);
%             alt_surf(pt_index(i_pt)) = chd.alt_trp(i_pt);
%             x_surf(pt_index(i_pt)) = x_trp(i_pt);
%             y_surf(pt_index(i_pt)) = y_trp(i_pt);
%             z_surf(pt_index(i_pt)) = z_trp(i_pt);
%             win_delay_surf(pt_index(i_pt)) = (alt_sat(pt_index(i_pt)) - alt_surf(pt_index(i_pt))) * 2 / cst.c;
            
            u = [x_surf(pt_index(i_pt)),y_surf(pt_index(i_pt)),z_surf(pt_index(i_pt))] - [x_surf(pt_index(i_pt)-1),y_surf(pt_index(i_pt)-1),z_surf(pt_index(i_pt)-1)]; 
            u_norm = u / norm(u);
            w = [x_trp(i_pt),y_trp(i_pt),z_trp(i_pt)] - [x_surf(pt_index(i_pt)-1),y_surf(pt_index(i_pt)-1),z_surf(pt_index(i_pt)-1)];
            cos_a = (u_norm * w.')/norm(w);  % norm(u_norm) = 1
            v = (norm(w) * cos_a) * u_norm;
            
            x_surf(pt_index(i_pt)) = x_surf(pt_index(i_pt)-1) + v(1);
            y_surf(pt_index(i_pt)) = y_surf(pt_index(i_pt)-1) + v(2);
            z_surf(pt_index(i_pt)) = z_surf(pt_index(i_pt)-1) + v(3);
            lla = ecef2lla([x_surf(pt_index(i_pt)),y_surf(pt_index(i_pt)),z_surf(pt_index(i_pt))],cst.flat_coeff,cst.semi_major_axis);
            L1BS.lat_surf(pt_index(i_pt)) = lla(1);
            L1BS.lon_surf(pt_index(i_pt)) = lla(2);
            L1BS.alt_surf(pt_index(i_pt)) = chd.alt_trp(i_pt);
            p = lla2ecef([L1BS.lat_surf(index),L1BS.lon_surf(index),L1BS.alt_surf(index)],cst.flat_coeff,cst.semi_major_axis);
            x_surf(index) = p(1);
            y_surf(index) = p(2);
            z_surf(index) = p(3);
            win_delay_surf(index) = (alt_sat(index) - L1BS.alt_surf(index)) * 2 / cst.c;

        end
    end
    
    
% elseif (0) %Transponder Flag 21A
% %     [x_trp,y_trp,z_trp,~,~] = gen_mec_con_03 (chd.lat_trp,chd.lon_trp,chd.alt_trp,cst.semi_major_axis,cst.flat_coeff);
%     p = lla2ecef([chd.lat_trp, chd.lon_trp, chd.alt_trp], cst.flat_coeff, cst.semi_major_axis);
%     x_trp = p(1);
%     y_trp = p(2);
%     z_trp = p(3);
%     
% % % %     u = [x_surf2(index+5),y_surf2(index+5),z_surf2(index+5)] - [x_surf2(index-5),y_surf2(index-5),z_surf2(index-5)];
% % % %     u_norm = u / norm(u);
% % % %     w = [x_trp(1),y_trp(1),z_trp(1)] - [x_surf2(index-5),y_surf2(index-5),z_surf2(index-5)];
% % % %     cos_a = cos((u_norm * w.')/norm(w));  % norm(u_norm) = 1
% % % %     v = (norm(w) * cos_a) * u_norm;
% % % %     x_surf2(index) = x_surf2(index-5) + v(1);
% % % % 	y_surf2(index) = y_surf2(index-5) + v(2);
% % % %     z_surf2(index) = z_surf2(index-5) + v(3);
% % % %     
% % % %     [lat_surf2(index),lon_surf2(index),alt_surf2(index)] = gen_mec_con_06 (x_surf2(index),y_surf2(index),z_surf2(index),cst.semi_major_axis,cst.flat_coeff,cnf.accu_orbit_alt,cnf.accu_orbit_lat);
% % % %     
%     
%     index = 72;
%     u = [x_surf(index),y_surf(index),z_surf(index)] - [x_surf(index-1),y_surf(index-1),z_surf(index-1)]; %perqu? la surface 177 s'ha d'endarrerir.
%                                                                                                          %si fos al rev?s, el vector 'u' tb canviaria.
%     u_norm = u / norm(u);
%     w = [x_trp(1),y_trp(1),z_trp(1)] - [x_surf(index-1),y_surf(index-1),z_surf(index-1)];
%     cos_a = (u_norm * w.')/norm(w);  % norm(u_norm) = 1
%     v = (norm(w) * cos_a) * u_norm;
%     
%     
%     % For normal processing
%     x_surf(index) = x_surf(index-1) + v(1);
%     y_surf(index) = y_surf(index-1) + v(2);
%     z_surf(index) = z_surf(index-1) + v(3);
% %     [lat_surf(index),L1BS.lon_surf(index),alt_surf(index)] = gen_mec_con_06(x_surf(index),y_surf(index),z_surf(index),cst.semi_major_axis,cst.flat_coeff,cnf.accu_orbit_alt,cnf.accu_orbit_lat);
%     lla = ecef2lla([x_surf(index),y_surf(index),z_surf(index)],cst.flat_coeff,cst.semi_major_axis);
%     lat_surf(index) = lla(1);
%     L1BS.lon_surf(index) = lla(2);
%     alt_surf(index) = lla(3);
%     win_delay_surf(index) = (alt_sat(index) - alt_surf(index)) * 2 / cst.c;
%     
% %     % For interpolation
% %     for i_interp = 1:10
% %         v_aux(i_interp) = (norm(w)*cos_a) * (9.87+i_interp/1000)/10;
% %         x_aux(i_interp) = x_surf(index-1) + v_aux(i_interp) * u_norm(1);
% %         y_aux(i_interp) = y_surf(index-1) + v_aux(i_interp) * u_norm(2);
% %         z_aux(i_interp) = z_surf(index-1) + v_aux(i_interp) * u_norm(3);
% %         [lat_aux(i_interp),lon_aux(i_interp),alt_aux(i_interp)] = gen_mec_con_06(x_aux(i_interp),y_aux(i_interp),z_aux(i_interp),cst.semi_major_axis,cst.flat_coeff,cnf.accu_orbit_alt,cnf.accu_orbit_lat);
% %     end
% %     
% %     win_delay_aux(i_interp) = (alt_sat(index) - alt_aux(i_interp)) * 2 / cst.c;

end

% % For interpolation
% x_test = [x_surf(1:176),x_aux,x_surf(177:N_total_surf_loc)];
% y_test = [y_surf(1:176),y_aux,y_surf(177:N_total_surf_loc)];
% z_test = [z_surf(1:176),z_aux,z_surf(177:N_total_surf_loc)];
% lat_test = [lat_surf(1:176),lat_aux,lat_surf(177:N_total_surf_loc)];
% lon_test = [L1BS.lon_surf(1:176),lon_aux,L1BS.lon_surf(177:N_total_surf_loc)];
% alt_test = [alt_surf(1:176),alt_aux,alt_surf(177:N_total_surf_loc)];
% win_test = [win_delay_surf(1:176),win_delay_aux,win_delay_surf(177:N_total_surf_loc)];
% N_total_surf_loc = N_total_surf_loc + 10;


% % For deletion of the nearest surface locations (la 177 s'ha substitu?t i la 176 s'elimina aqu?)
% x_test = [x_surf(1:175),x_surf(177:N_total_surf_loc)];
% y_test = [y_surf(1:175),y_surf(177:N_total_surf_loc)];
% z_test = [z_surf(1:175),z_surf(177:N_total_surf_loc)];
% lat_test = [lat_surf(1:175),lat_surf(177:N_total_surf_loc)];
% lon_test = [L1BS.lon_surf(1:175),L1BS.lon_surf(177:N_total_surf_loc)];
% alt_test = [alt_surf(1:175),alt_surf(177:N_total_surf_loc)];
% win_test = [win_delay_surf(1:175),win_delay_surf(177:N_total_surf_loc)];
% N_total_surf_loc = N_total_surf_loc - 1;



% x_surf = x_test;
% y_surf = y_test;
% z_surf = z_test;
% lat_surf = lat_test;
% L1BS.lon_surf = lon_test;
% alt_surf = alt_test;
% win_delay_surf = win_test;



