% Read CESM FOSI netcdf

clear all
close all

fpath='/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/';

%% All inputs
ncdisp([fpath 'cyclic_test_forcing.nc'])

%%
% Size:       22x366
% Dimensions: location,time
tp_long_name   = 'T_pelagic';
tp_units       = 'degC';
tb_long_name   = 'T_bottom';
tb_units       = 'degC';
det_long_name  = 'POC flux';
det_units      = 'g/m^2/d';
det_b          = 0.7;
zoo_long_name  = 'Zooplankton biomass';
zoo_units      = 'g/m^2';
mort_long_name = 'Zooplankton quadratic mortality';
mort_units     = 'g/m^2/d';
time_units = 'days';

%%
ncid = netcdf.open([fpath 'cyclic_test_forcing.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:5
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
for i = 7:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% bathym
ncdisp([fpath 'bathym_test_forcing.nc'])

%%
% Size:       22x366
% Dimensions: location,time
FillValue   = NaN;
bathy_long_name  = 'depth';
bathy_units      = 'm';

%%
ncid = netcdf.open([fpath 'bathym_test_forcing.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
save([fpath 'cyclic_test_forcing.mat']);

%% struct for code
ESM.Tp = T_pelagic(:,1:365);
ESM.Tb = T_bottom(:,1:365);
ESM.Zm = zooC(:,1:365);
ESM.det = poc_flux_bottom(:,1:365);
ESM.dZm = zoo_mort(:,1:365);
GRD.Z = bathymetry;
GRD.area = nan*ones(size(bathymetry));

save([fpath 'Data_cyclic_test_forcing.mat'],'ESM');
save([fpath 'Grid_test_forcing.mat'],'GRD');




