% compute_height_rate
% It computes an approximation of the height rate (in m/s)
% v1.1 2019/06/10 Changed method to compute the angle from acos(v_norm*m_norm'); to real(-1.0*acos(v_norm*m_norm'));

function [height_rate,a] = compute_height_rate(N_total_burst,...
                                           x_vel_sat_geoloc, y_vel_sat_geoloc, z_vel_sat_geoloc, ...
                                           x_sat_geoloc, y_sat_geoloc, z_sat_geoloc, ...
                                           x_surf_geoloc, y_surf_geoloc, z_surf_geoloc, cst)

a = zeros(1,N_total_burst);
height_rate = zeros(1,N_total_burst);

for i_burst = 1:N_total_burst
    % nadir direction
    n = [x_surf_geoloc(i_burst) - x_sat_geoloc(i_burst), ...
         y_surf_geoloc(i_burst) - y_sat_geoloc(i_burst), ...
         z_surf_geoloc(i_burst) - z_sat_geoloc(i_burst)];
    n_norm = n/norm(n);
    
    % velocity vector
    v = [x_vel_sat_geoloc(i_burst),y_vel_sat_geoloc(i_burst),z_vel_sat_geoloc(i_burst)];
    v_norm = v/norm(v);

    % Vector perpendicular a 'n' i 'v': 'w' (vector associat al pla A)
    w = cross(n,v);
    w_norm = w/norm(w);

    % Vector perpendicular a 'n' i 'w': 'm'
    m = cross(w_norm,n_norm);
    m_norm = m/norm(m);

   % angle between 'v' and 'm'
    % issue getting complex angles
   % angle = acos(v_norm*m_norm')
  
   % angle between 'v' and 'm' 
   % angle between 'v' and 'm' 
   if acos(dot(v_norm,n_norm))<cst.pi/2 
       a(i_burst) = real(-1.0*acos(v_norm*m_norm')); 
   else
       a(i_burst) = real(1.0*acos(v_norm*m_norm')); 
   end
   
%    if atan2(dot(v_norm,n_norm))<cst.pi/2 
%        a(i_burst) = -1.0*atan2(norm(cross(v,m)),dot(v,m)); 
%    else
%        a(i_burst) = 1.0*atan2(norm(cross(v,m)),dot(v,m)); 
%    end
   
%    if acos(dot(v_norm,n_norm))<cst.pi/2 
%        a(i_burst) = -1.0*atan2(norm(cross(v,m)),dot(v,m)); 
%    else
%        a(i_burst) = 1.0*atan2(norm(cross(v,m)),dot(v,m)); 
%    end
%    a(i_burst) = 1.0*atan2(norm(cross(v,m)),dot(v,m));
    % height rate
    height_rate(i_burst) = norm(v) * sin(a(i_burst));
end
