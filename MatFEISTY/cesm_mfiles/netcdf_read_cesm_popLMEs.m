% Read CESM FOSI netcdf

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';

%% All inputs
ncdisp([fpath 'LME-mask-POP_gx1v6.nc'])

%%
ncid = netcdf.open([fpath 'LME-mask-POP_gx1v6.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%%
save([fpath 'LME-mask-POP_gx1v6.mat']);





