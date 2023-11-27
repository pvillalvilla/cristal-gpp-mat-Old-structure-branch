tpt1 = 'G:/My Drive/S3NGT/codes/S3NGT_TPT1/inputs//S3NGT_SIRS___KUHR__1A_20230702T213953_20230702T213955_0001.nc';
tpt2 = 'G:\My Drive\S3NGT\Aresys\TPT2_grid\PIC_SIRS___KUHR__1A_20280101T234415_20280101T234425_0001.nc';


% Open the NetCDF file in write mode
ncid_copy = netcdf.open(tpt2, 'NC_NOWRITE');

% Read existing data
id_0 = netcdf.inqVarID(ncid_copy, 'data_record_time');
data0 = netcdf.getVar(ncid_copy, id_0);

id_1 = netcdf.inqVarID(ncid_copy, 'com_position_vector');
data1 = netcdf.getVar(ncid_copy, id_1);

id_2 = netcdf.inqVarID(ncid_copy, 'com_velocity_vector');
data2 = netcdf.getVar(ncid_copy, id_2);

id_3 = netcdf.inqVarID(ncid_copy, 'com_altitude_rate');
data3 = netcdf.getVar(ncid_copy, id_3);

% Open the NetCDF file in write mode
ncid = netcdf.open(tpt1, 'NC_WRITE');

% Read existing data
varID0 = netcdf.inqVarID(ncid, 'data_record_time');
var0 = netcdf.getVar(ncid, varID0);
varID1 = netcdf.inqVarID(ncid, 'com_position_vector');
var1 = netcdf.getVar(ncid, varID1);
varID2 = netcdf.inqVarID(ncid, 'com_velocity_vector');
var2 = netcdf.getVar(ncid, varID2);
varID3 = netcdf.inqVarID(ncid, 'com_altitude_rate');
var3 = netcdf.getVar(ncid, varID3);

% Write modified data back to the NetCDF file
netcdf.putVar(ncid, varID1, data2(:,1:376));


% Close the NetCDF file
netcdf.close(ncid);
netcdf.close(ncid_copy);
