function    [chd] = check_dimensions(filesBulk,chd, level)

switch level
    case 'L0'
        filesBulk.fid = netcdf.open(filesBulk.filename_L0,'NC_NOWRITE'); %open file
      case 'L1A' 
          filesBulk.fid = netcdf.open(filesBulk.filename_L1A,'NC_NOWRITE'); %open file
    
end
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(filesBulk.fid); %get global attributes


for i_dim=0:(ndims-1)
    [dimname, dimlen(i_dim+1)] = netcdf.inqDim(filesBulk.fid,i_dim);
    if(strcmp(dimname,'np'))
        chd.N_pulses_burst     = dimlen(i_dim+1);
    end
    if(strcmp(dimname,'ns'))
        chd.N_samples_sar     = dimlen(i_dim+1);
    end
end
netcdf.close(filesBulk.fid); %close the netcdf file

end
