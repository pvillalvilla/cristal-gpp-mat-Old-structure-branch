%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: Computation of the positions, velocities and attitude of
% satellite per pulse
%
% Calling: 
% INPUTs:
%   L1A_buffer: Buffer of L1A structures of all the bursts
%   cnf:        Configuration parameters structure
%   chd:        Characterization parameters structure
%   cst:        Constant parameters structure
%
% OUTPUTs:
%  x/y/z_sar_sat_pulse:     x/y/z positions of the satellite per pulse 
%  x/y/z_vel_sat_sar_pulse: x/y/z satellite velocity comp. per pulse 
%  lat_sar_sat_pulse:       latitude of satellite per pulse 
%  lon_sar_sat_pulse:       longitude of satellite per pulse
%  alt_sar_sat_pulse:       altitude of the satellite per pulse
%  alt_rate_sar_sat_pulse:  height rate of the satellite per pulse
%  roll_sar_pulse:          roll of the satellite per pulse
%  pitch_sar_pulse:         pitch of the satellite per pulse
%  yaw_sar_pulse:           yaw of the satellite per pulse
%
%
%
% COMMENTS: Current version uses a spline interpolator using the time and OSV info of all the bursts 
% called with the datation for all the pulses
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%            Albert Garcia / isardSAT
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_sar_sat_pulse, y_sar_sat_pulse, z_sar_sat_pulse,...
 x_vel_sat_sar_pulse,  y_vel_sat_sar_pulse,  z_vel_sat_sar_pulse,...
 lat_sar_sat_pulse, lon_sar_sat_pulse,...
 alt_sar_sat_pulse,alt_rate_sar_sat_pulse,...
 roll_sar_pulse, pitch_sar_pulse, yaw_sar_pulse]=...
 osv_attitude_pulse(L1A_buffer,time_pulse, cnf, chd, cst)

% % ADRIAN method based on S3 approach
% %% V2: 1. Work with offset time vector (It does not look it improves anything)
% %      2. Fit 5th order polynomial to spacecraft dynamics (location only?)
%ADRIAN test S3 approach, setting 1st time value to 0
order_p = 5;
% % Reset time vector
 time_vec_res = [L1A_buffer.time_rx_1st] - [L1A_buffer(1).time_rx_1st];
% %time_vec_res = time_pulse - time_pulse(1);
% % time_offset  = [L1A_buffer(1).time];
 x_sar_sat_polycoef = polyfit(time_vec_res, [L1A_buffer.x_sar_sat], order_p);
 y_sar_sat_polycoef = polyfit(time_vec_res, [L1A_buffer.y_sar_sat], order_p);
 z_sar_sat_polycoef = polyfit(time_vec_res, [L1A_buffer.z_sar_sat], order_p);
 x_sar_sat_int = polyval(x_sar_sat_polycoef, time_vec_res);
 y_sar_sat_int = polyval(y_sar_sat_polycoef, time_vec_res);
 z_sar_sat_int = polyval(z_sar_sat_polycoef, time_vec_res);
% %cartesian coordinates
 x_sar_sat_pulse = spline(time_vec_res, x_sar_sat_int, time_pulse);
 y_sar_sat_pulse = spline(time_vec_res, y_sar_sat_int, time_pulse);
 z_sar_sat_pulse = spline(time_vec_res, z_sar_sat_int, time_pulse);
% %geodetic coordinates
 lla = ecef2lla([x_sar_sat_pulse.', y_sar_sat_pulse.', z_sar_sat_pulse.'], 'WGS84');
% %lla = ecef2lla([x_sar_sat_pulse.', y_sar_sat_pulse.', z_sar_sat_pulse.'], cst.flat_coeff,cst.semi_major_axis);
lat_sar_sat_pulse = lla(:,1).';
lon_sar_sat_pulse = lla(:,2).';
alt_sar_sat_pulse = lla(:,3).';
% Alternatively?
% lla = ecef2lla([[L1A_buffer.x_sar_sat].', [L1A_buffer.y_sar_sat].', [L1A_buffer.z_sar_sat].'], 'WGS84');
% 
% lat_polycoef = polyfit(lla(:,1).', time_vec_res, order_p);
% lon_polycoef = polyfit(lla(:,2).', time_vec_res, order_p);
% alt_polycoef = polyfit(lla(:,3).', time_vec_res, order_p);
% 
% lat_int = polyval(lat_polycoef, time_vec_res);
% lon_int = polyval(lon_polycoef, time_vec_res);
% alt_int = polyval(alt_polycoef, time_vec_res);
% 
% lat_sar_sat_pulse = spline(time_vec_res, lat_int, time_pulse);
% lon_sar_sat_pulse = spline(time_vec_res, lon_int, time_pulse);
% alt_sar_sat_pulse = spline(time_vec_res, alt_int, time_pulse);
% 
% xyz_sar = lla2ecef([lat_sar_sat_pulse.', lon_sar_sat_pulse.', alt_sar_sat_pulse.'], 'WGS84');
% 
% x_sar_sat_pulse = xyz_sar(:,1).';
% y_sar_sat_pulse = xyz_sar(:,2).';
% z_sar_sat_pulse = xyz_sar(:,3).';
%% COMPUTATION OF VELOCITIES
x_vel_sar_sat_polycoef = polyfit(time_vec_res, [L1A_buffer.x_vel_sat_sar], order_p);
y_vel_sar_sat_polycoef = polyfit(time_vec_res, [L1A_buffer.y_vel_sat_sar], order_p);
z_vel_sar_sat_polycoef = polyfit(time_vec_res, [L1A_buffer.z_vel_sat_sar], order_p);
x_vel_sar_sat_int = polyval(x_vel_sar_sat_polycoef, time_vec_res);
y_vel_sar_sat_int = polyval(y_vel_sar_sat_polycoef, time_vec_res);
z_vel_sar_sat_int = polyval(z_vel_sar_sat_polycoef, time_vec_res);
x_vel_sat_sar_pulse = spline(time_vec_res, x_vel_sar_sat_int, time_pulse);
y_vel_sat_sar_pulse = spline(time_vec_res, y_vel_sar_sat_int, time_pulse);
z_vel_sat_sar_pulse = spline(time_vec_res, z_vel_sar_sat_int, time_pulse);
alt_rate_sar_sat_polycoef = polyfit(time_vec_res, [L1A_buffer.alt_rate_sar_sat], order_p);
alt_rate_sar_sat_int = polyval(alt_rate_sar_sat_polycoef, time_vec_res);
alt_rate_sar_sat_pulse = spline(time_vec_res, alt_rate_sar_sat_int, time_pulse);
%% COMPUTATION OF ATTITUDE
roll_sar_sat_polycoef  = polyfit(time_vec_res, [L1A_buffer.roll_sar], order_p);
pitch_sar_sat_polycoef = polyfit(time_vec_res, [L1A_buffer.pitch_sar], order_p);
yaw_sar_sat_polycoef   = polyfit(time_vec_res, [L1A_buffer.yaw_sar], order_p);
roll_sar_sat_int  = polyval(roll_sar_sat_polycoef, time_vec_res);
pitch_sar_sat_int = polyval(pitch_sar_sat_polycoef, time_vec_res);
yaw_sar_sat_int   = polyval(yaw_sar_sat_polycoef, time_vec_res);
roll_sar_pulse  = spline(time_vec_res, roll_sar_sat_int, time_pulse);
pitch_sar_pulse = spline(time_vec_res, pitch_sar_sat_int, time_pulse);
yaw_sar_pulse   = spline(time_vec_res, yaw_sar_sat_int, time_pulse);

% % OLD METHOD
% %% COMPUTATION OF POSITIONS
% %cartesian coordinates
% x_sar_sat_pulse = spline([L1A_buffer.time],[L1A_buffer.x_sar_sat],time_pulse);
% y_sar_sat_pulse = spline([L1A_buffer.time],[L1A_buffer.y_sar_sat],time_pulse);
% z_sar_sat_pulse = spline([L1A_buffer.time],[L1A_buffer.z_sar_sat],time_pulse);
% 
% %geodetic coordinates
% %lla = ecef2lla([x_sar_sat_pulse.',y_sar_sat_pulse.',z_sar_sat_pulse.'],cst.flat_coeff,cst.semi_major_axis);
% lla = ecef2lla([x_sar_sat_pulse.', y_sar_sat_pulse.', z_sar_sat_pulse.'], 'WGS84');
% lat_sar_sat_pulse = lla(:,1).';
% lon_sar_sat_pulse = lla(:,2).';
% alt_sar_sat_pulse = lla(:,3).';
% 
% 
% %% COMPUTATION OF VELOCITIES
% x_vel_sat_sar_pulse = spline([L1A_buffer.time],[L1A_buffer.x_vel_sat_sar],time_pulse);
% y_vel_sat_sar_pulse = spline([L1A_buffer.time],[L1A_buffer.y_vel_sat_sar],time_pulse);
% z_vel_sat_sar_pulse = spline([L1A_buffer.time],[L1A_buffer.z_vel_sat_sar],time_pulse);
% 
% alt_rate_sar_sat_pulse = spline([L1A_buffer.time],[L1A_buffer.alt_rate_sar_sat],time_pulse);
% 
% %% COMPUTATION OF ATTITUDE
% roll_sar_pulse  = spline([L1A_buffer.time],[L1A_buffer.roll_sar],time_pulse);
% pitch_sar_pulse = spline([L1A_buffer.time],[L1A_buffer.pitch_sar],time_pulse);
% yaw_sar_pulse   = spline([L1A_buffer.time],[L1A_buffer.yaw_sar],time_pulse);


end

