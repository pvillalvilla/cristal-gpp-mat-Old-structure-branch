function [datation_bias_value, datation_bias_req_met, pos_max_along_focused] = datation_bias_SL(wfm_AC_interp, lat_surf,...
    lon_surf, x_vel_sat_sar_pulse, y_vel_sat_sar_pulse, z_vel_sat_sar_pulse, alt_sar_sat_pulse, ...
    cnf, chd, cst)
    % Compute distances along and slant range w.r.t TRP location
    % along
    [arclen,~]=distance(lat_surf,lon_surf,ones(1,length(lat_surf))*chd.lat_trp,ones(1,length(lat_surf))*chd.lon_trp,[cst.semi_major_axis sqrt(1-(cst.semi_minor_axis/cst.semi_major_axis).^2)]);
    [~,min_dist]=min(arclen);
    distances_over_arclengths=arclen;
    distances_over_arclengths(1:min_dist-1)=-1.0*distances_over_arclengths(1:min_dist-1);
    distances_over_arclengths_interp=interp(distances_over_arclengths,cnf.FFt.zp_al_TRP);
    
    % Along and across track cuts of the 2D PTR
    % along and across track cuts in the maximums
    [~,pos_max] = max(abs(wfm_AC_interp(:)));
    [pos_max_along_focused,~] = ind2sub(size(wfm_AC_interp),pos_max);
    
    %Along-track location error (position of peak in along-track);
    Along_track_error = abs(distances_over_arclengths_interp(pos_max_along_focused));
    %Datation error: translate along-track error location to time with
    %the ground velocity of the beam
    %compute mean values of velocities and orbital height
    vs=mean(sqrt(sum(([x_vel_sat_sar_pulse.',y_vel_sat_sar_pulse.',z_vel_sat_sar_pulse.']).^2,2)));
    H_orb = mean(alt_sar_sat_pulse);
    vg = vs/(1+H_orb/cst.earth_radius);
    
    datation_bias_value = Along_track_error/vg;
    datation_bias_req_met = datation_bias_value < chd.datation_bias; 

end