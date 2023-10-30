% compute_height_rate
% It computes an approximation of the height rate (in m/s)

function [height_rate,a] = compute_height_rate_generic(x_vel_sat_geoloc, y_vel_sat_geoloc, z_vel_sat_geoloc, ...
                                           x_sat_geoloc, y_sat_geoloc, z_sat_geoloc, ...
                                           x_surf_geoloc, y_surf_geoloc, z_surf_geoloc)


global pi_cst FFS_processing_active FFS_OSV_att_selec_pulse

N_points = length(x_vel_sat_geoloc);

a = zeros(1,N_points);
height_rate = zeros(1,N_points);

for i_point = 1:N_points
    if FFS_processing_active & FFS_OSV_att_selec_pulse
        progressbar([],[],[],[],[],[],[],[],[],i_point/N_points);
    end

    % nadir direction
    n = [x_surf_geoloc(i_point) - x_sat_geoloc(i_point), ...
         y_surf_geoloc(i_point) - y_sat_geoloc(i_point), ...
         z_surf_geoloc(i_point) - z_sat_geoloc(i_point)];
    n_norm = n/norm(n);
    
    % velocity vector
    v = [x_vel_sat_geoloc(i_point),y_vel_sat_geoloc(i_point),z_vel_sat_geoloc(i_point)];
    v_norm = v/norm(v);

    % Vector perpendicular a 'n' i 'v': 'w' (vector associat al pla A)
    w = cross(n,v);
    w_norm = w/norm(w);

    % Vector perpendicular a 'n' i 'w': 'm'
    m = cross(w_norm,n_norm);
    m_norm = m/norm(m);

    % angle between 'v' and 'm'
    a(i_point) = acos(dot(v_norm,m_norm));
    
    if (acos(dot(v_norm,n_norm)) < pi_cst / 2)

	a(i_point) = -a(i_point);
    end

    % height rate
    height_rate(i_point) = norm(v) * sin(a(i_point));
end
