%% 
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% DeDop 
% This code implements the CODING & PACKING 
% algorithm for Level-1ASS product for input Fully Focused Chain
%
% ---------------------------------------------------------
% Objective: Pack variables and write the NETCDF
% 
% INPUTs : Workspace
% OUTPUTs: TM Structure as defined on isardSAT_JasonCS_DPM
%
% ----------------------------------------------------------
% Author:    Albert Garcia / isardSAT
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT ()
% 
% Versions
% 1.0 
% 1.1 Updated time conversion for data between 2010 and 2016 (CR2)
% 2.0 Transformed to a function. Writting one record
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NetCDF utils
function write_NetCDF_L1ASS_surf(files, time_surf, x_surf, y_surf, z_surf, cnf,chd,cst)

%% PACKING L1A

%----------A. Time variables ----------------------------------------------
var_id = netcdf.inqVarID(files.ncid_L1AS,'time_l1ass_echo');
netcdf.putVar(files.ncid_L1AS,var_id,time_surf);
%----------D. Surface Position variables ------------------------------
var_id = netcdf.inqVarID(files.ncid_L1AS,'x_pos_l1ass_echo');
netcdf.putVar(files.ncid_L1AS,var_id,x_surf);
var_id = netcdf.inqVarID(files.ncid_L1AS,'y_pos_l1ass_echo');
netcdf.putVar(files.ncid_L1AS,var_id,y_surf);
var_id = netcdf.inqVarID(files.ncid_L1AS,'z_pos_l1ass_echo');
netcdf.putVar(files.ncid_L1AS,var_id,z_surf);


%------------- N. CHD Parameters -----------------------------

var_id = netcdf.inqVarID(files.ncid_L1AS,'freq');
netcdf.putVar(files.ncid_L1AS,var_id,chd.freq);
var_id = netcdf.inqVarID(files.ncid_L1AS,'bw');
netcdf.putVar(files.ncid_L1AS,var_id,chd.bw);
var_id = netcdf.inqVarID(files.ncid_L1AS,'pulse_length');
netcdf.putVar(files.ncid_L1AS,var_id,chd.pulse_length);
var_id = netcdf.inqVarID(files.ncid_L1AS,'chirp_slope');
netcdf.putVar(files.ncid_L1AS,var_id,chd.chirp_slope);
var_id = netcdf.inqVarID(files.ncid_L1AS,'alt_freq');
netcdf.putVar(files.ncid_L1AS,var_id,chd.alt_freq);
var_id = netcdf.inqVarID(files.ncid_L1AS,'N_bursts_cycle');
netcdf.putVar(files.ncid_L1AS,var_id,chd.N_bursts_cycle_sar);
var_id = netcdf.inqVarID(files.ncid_L1AS,'N_samples_sar');
netcdf.putVar(files.ncid_L1AS,var_id,chd.N_samples_sar);
var_id = netcdf.inqVarID(files.ncid_L1AS,'N_pulses_burst');
netcdf.putVar(files.ncid_L1AS,var_id,chd.N_pulses_burst);
var_id = netcdf.inqVarID(files.ncid_L1AS,'antenna_beamwidth_along_track');
netcdf.putVar(files.ncid_L1AS,var_id,chd.antenna_beamwidth_along_track);
var_id = netcdf.inqVarID(files.ncid_L1AS,'antenna_beamwidth_across_track');
netcdf.putVar(files.ncid_L1AS,var_id,chd.antenna_beamwidth_across_track);




netcdf.close(files.ncid_L1AS);


