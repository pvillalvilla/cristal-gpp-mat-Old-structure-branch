%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v1a
%
% ---------------------------------------------------------
% Objective: Read ORBIT from L0
%
% Calling: 
% INPUTs:
%
%
% OUTPUTs:
%
%
% ----------------------------------------------------------
% Author:    Albert Garcia  / isardSAT
%            
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ORBIT] = read_orbit(files)
    
        ORBIT.com_position_vector   = ncread(files.filename_L0,'com_position_vector');
        ORBIT.com_velocity_vector   = ncread(files.filename_L0,'com_velocity_vector');
        ORBIT.alt_rate_sar_sat      = ncread(files.filename_L0,'com_altitude_rate');
        
        ORBIT.roll                  = ncread(files.filename_L0,'off_nadir_roll_angle');
        ORBIT.pitch                 = ncread(files.filename_L0,'off_nadir_pitch_angle');
        ORBIT.yaw                   = ncread(files.filename_L0,'off_nadir_yaw_angle');
        ORBIT.tracker_range         = ncread(files.filename_L0,'tracker_range_calibrated');
        ORBIT.time                  = ncread(files.filename_L0,'data_record_time');
        
end

