%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    isardSAT S.L.   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% readanyNETCDF reads any file in NETCDF format
% INPUT example: 'C:\Users\Pablo\Desktop\RA2\PTR\NETDCFexample.nc'
% OUTPUT: a Matlab structure with all the data and attributes from the NETCDF file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alldata]=readanyNETCDF_V2(file)

ncid = netcdf.open(file,'NC_NOWRITE'); %open file
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid); %get global attributes
[~,nvars,~,~] = netcdf.inq(ncid); %get global attributes

childGrps = netcdf.inqGrps(ncid); %retrieve groups of variables

% retrieve variables in root
for i=0:(nvars-1)
    [data(i+1).varname,data(i+1).xtype,data(i+1).dimids,data(i+1).natts] = netcdf.inqVar(ncid,i); %get variables name
    alldata.data.(data(i+1).varname) = netcdf.getVar(ncid,netcdf.inqVarID(ncid,data(i+1).varname)); %get variable data
    for j=0:(data(i+1).natts - 1)
        att.name = netcdf.inqAttName(ncid,netcdf.inqVarID(ncid,data(i+1).varname),j); %get attributes name
        att.name2put = att.name;
        if att.name(1) == '_'
           att.name2put(1) = []; 
        end
        alldata.attributes.(data(i+1).varname).(att.name2put) = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,data(i+1).varname),(att.name)); %get attributes info
    end
end
clear data

% retrieve variables in groups
for g=1:length(childGrps)
    [~,nvars,~,~] = netcdf.inq(childGrps(g)); %get group 'g' attributes
    groupName     = netcdf.inqGrpNameFull(childGrps(g)); %get group name
    groupName     = groupName(2:end);
    for i=0:(nvars-1)
        [data(i+1).varname,~,~,data(i+1).natts] = netcdf.inqVar(childGrps(g),i); %get variables name
        alldata.data.(groupName).(data(i+1).varname) = netcdf.getVar(childGrps(g),netcdf.inqVarID(childGrps(g),data(i+1).varname)); %get variable data
        for j=0:(data(i+1).natts - 1)
            att.name = netcdf.inqAttName(childGrps(g),netcdf.inqVarID(childGrps(g),data(i+1).varname),j); %get attributes name
            att.name2put = att.name;
            if att.name(1) == '_'
                att.name2put(1) = [];
            end
            alldata.attributes.(groupName).(data(i+1).varname).(att.name2put) = netcdf.getAtt(childGrps(g),netcdf.inqVarID(childGrps(g),data(i+1).varname),(att.name)); %get attributes info
        end
    end
end

% finfo = ncinfo(file);
% alldata.attributes.global = finfo.Attributes;

netcdf.close(ncid);

end