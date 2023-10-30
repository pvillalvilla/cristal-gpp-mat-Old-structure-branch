%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: Get number of records of a netCDF      
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
function N_records=get_num_record(filename, string)
ncid = netcdf.open(filename,'NC_NOWRITE'); %open file
[ndims,~,~,~] = netcdf.inq(ncid); %get global attributes
for i_dim=0:(ndims-1)
    [dimname, dimlen] = netcdf.inqDim(ncid,i_dim);
    if(strcmp(dimname,string))
        break;
    end
end
netcdf.close(ncid);
N_records = dimlen;