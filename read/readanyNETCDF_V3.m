%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    isardSAT S.L.   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% readanyNETCDF reads any file in NETCDF format
% INPUT example: 'C:\Users\Pablo\Desktop\RA2\PTR\NETDCFexample.nc'
% OUTPUT: a Matlab structure with all the data and attributes from the NETCDF file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alldata]=readanyNETCDF_V3(file)

ncid = netcdf.open(file,'NC_NOWRITE'); %open file
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid); %get global attributes
[~,nvars,~,~] = netcdf.inq(ncid); %get global attributes

childGrps = netcdf.inqGrps(ncid); %retrieve groups of variables

% retrieve variables in root
for i=0:(nvars-1)
    [data(i+1).varname,data(i+1).xtype,data(i+1).dimids,data(i+1).natts] = netcdf.inqVar(ncid,i); %get variables name
    alldata.data.(data(i+1).varname) = ncread(file,  (data(i+1).varname));

end
clear data

% retrieve variables in groups
for g=1:length(childGrps)
    [~,nvars,~,~] = netcdf.inq(childGrps(g)); %get group 'g' attributes
    groupName     = netcdf.inqGrpNameFull(childGrps(g)); %get group name
    groupName     = groupName(2:end);
    for i=0:(nvars-1)
        [data(i+1).varname,~,~,data(i+1).natts] = netcdf.inqVar(childGrps(g),i); %get variables name
%          alldata.data.(groupName).(data(i+1).varname) = netcdf.getVar(childGrps(g),netcdf.inqVarID(childGrps(g),data(i+1).varname)); %get variable data
        alldata.(groupName).(data(i+1).varname) = ncread(file, ['/' groupName '/' (data(i+1).varname)]);

    end
    childGrps2 = netcdf.inqGrps(childGrps(g)); % mirar si hi ha grups dins dels groups
    if length(childGrps2) >= 1 
        for n=1:length(childGrps2)
            [~,nvars,~,~] = netcdf.inq(childGrps2(n)); %get group 'g' attributes
            groupNameFull     = netcdf.inqGrpNameFull(childGrps2(n)); %get group name
            subgroupName = split(groupNameFull(2:end),'/');
%             groupNameFu     = groupNameFull(2:end);
            for i=0:(nvars-1)
                [data(i+1).varname,~,~,data(i+1).natts] = netcdf.inqVar(childGrps2(n),i); %get variables name
        %          alldata.data.(groupName).(data(i+1).varname) = netcdf.getVar(childGrps(g),netcdf.inqVarID(childGrps(g),data(i+1).varname)); %get variable data
                alldata.(subgroupName{1}).(subgroupName{2}).(data(i+1).varname) = ncread(file, [groupNameFull '/' (data(i+1).varname)]);
            end
        end
    end
end

% finfo = ncinfo(file);
% alldata.attributes.global = finfo.Attributes;

netcdf.close(ncid);
 end