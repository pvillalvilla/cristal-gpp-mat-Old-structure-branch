%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: Read L0 record
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
%            Eduard Makhoul / isardSAT
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L1A,files] = read_L1A_record(files,L1A,i_burst,N_bursts, cnf, chd, cst)
    
        if i_burst==1
            files.fidL1A = netcdf.open(files.filename_L1A,'NC_NOWRITE'); %open file
        end
        [netCDF_L1A] = readanyNETCDF_record(files.fidL1A,i_burst,'nb');
        [L1A]        = adaptL1A(netCDF_L1A,files.filename_L1A,L1A,i_burst, cnf, chd, cst);
        if i_burst==N_bursts
            netcdf.close(files.fidL1A); %close the netcdf file
        end
    
 
    %fclose all;
end

