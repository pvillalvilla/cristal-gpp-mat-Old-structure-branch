%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    isardSAT S.L.   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% readanyNETCDF reads any file in NETCDF format
% INPUT example: 'C:\Users\Pablo\Desktop\RA2\PTR\NETDCFexample.nc'
% OUTPUT: a Matlab structure with all the data and attributes from the NETCDF file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alldata]=readanyNETCDF_V2(file)

ncid = netcdf.open(file,'NC_NOWRITE'); %open file
[~,nvars,~,~] = netcdf.inq(ncid); %get global attributes

for i=0:(nvars-1)
    [varname,~,~,natts] = netcdf.inqVar(ncid,i); %get variables name
    alldata.data.(varname) = netcdf.getVar(ncid,netcdf.inqVarID(ncid,varname)); %get variable data
    for j=0:(natts - 1)
        name = netcdf.inqAttName(ncid,netcdf.inqVarID(ncid,varname),j); %get attributes name
        name2put = name;
        if name(1) == '_'
            name2put(1) = [];
        end
        alldata.attributes.(varname).(name2put) = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,varname),(name)); %get attributes info
    end
end;

netcdf.close(ncid);

end

% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid)            Return information about netCDF file

% [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,varid)      Return information about variable
% data = netcdf.getVar(ncid,varid)                              Return data from netCDF variable
% varid = netcdf.inqVarID(ncid,varname)                         Return ID associated with variable name

% attname = netcdf.inqAttName(ncid,varid,attnum)                Return name of netCDF attribute
% [xtype,attlen] = netcdf.inqAtt(ncid,varid,attname)            Return information about netCDF attribute
% attrvalue = netcdf.getAtt(ncid,varid,attname)                 Return netCDF attribute