tpt1 = 'G:/My Drive/S3NGT/codes/S3NGT_TPT1/inputs/S3NGT_SIRS___KUHR__1A_20230702T213953_20230702T213955_0001.nc';
tpt2 = 'G:\My Drive\S3NGT\codes\S3NGT_TPT2\inputs\S3NGT_SIRS___KUHR__1A_20280101T234415_20280101T234425_0001.nc';


% Open the NetCDF file in write mode
ncid_tpt2 = netcdf.open(tpt2, 'NC_NOWRITE');

% Read existing data
id_0 = netcdf.inqVarID(ncid_tpt2, 'data_record_time');
data0 = netcdf.getVar(ncid_tpt2, id_0);

id_1 = netcdf.inqVarID(ncid_tpt2, 'com_position_vector');
data1 = netcdf.getVar(ncid_tpt2, id_1);

id_2 = netcdf.inqVarID(ncid_tpt2, 'com_velocity_vector');
data2 = netcdf.getVar(ncid_tpt2, id_2);

id_3 = netcdf.inqVarID(ncid_tpt2, 'com_altitude_rate');
data3 = netcdf.getVar(ncid_tpt2, id_3);

% Open the NetCDF file in write mode
ncid_tpt1 = netcdf.open(tpt1, 'NC_NOWRITE');

% Read existing data
varID0 = netcdf.inqVarID(ncid_tpt1, 'data_record_time');
var0 = netcdf.getVar(ncid_tpt1, varID0);
varID1 = netcdf.inqVarID(ncid_tpt1, 'com_position_vector');
var1 = netcdf.getVar(ncid_tpt1, varID1);
varID2 = netcdf.inqVarID(ncid_tpt1, 'com_velocity_vector');
var2 = netcdf.getVar(ncid_tpt1, varID2);
varID3 = netcdf.inqVarID(ncid_tpt1, 'com_altitude_rate');
var3 = netcdf.getVar(ncid_tpt1, varID3);

% Write modified data back to the NetCDF file
netcdf.putVar(ncid_tpt1, varID1, data2(:,1:376));


% Close the NetCDF file
netcdf.close(ncid_tpt1);
netcdf.close(ncid_tpt2);
