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

function [L1A,files] = read_adapt_L0_record(files,L1A,i_burst,N_bursts, cnf, chd, cst)
    
        if i_burst==1
            files.fidL0 = netcdf.open(files.filename_L0,'NC_NOWRITE'); %open file
        end
        [netCDF_L0] = readanyNETCDF_record(files.fidL0,i_burst,'nb');
        [L1A]        = adaptL0_to_L1A(netCDF_L0,files.filename_L0,L1A,i_burst, cnf, chd, cst);
        if i_burst==N_bursts
            netcdf.close(files.fidL0); %close the netcdf file
        end
    
 
    %fclose all;
end

