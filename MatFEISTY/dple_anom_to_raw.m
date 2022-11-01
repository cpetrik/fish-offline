% Function to read in DPLE member
function [TEMP_150m,TEMP_bottom,POC_FLUX_IN_bottom,LzooC_150m,Lzoo_loss_150m] ...
    = dple_anom_to_raw(fpath,im,yr,ni,nj,L,clim_Tp,clim_Tb,clim_POC,clim_zooC,clim_loss)

% READ IN ENSEMBLE MEMBER
% Pelagic temperature
ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_TEMP_150m.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1, [0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
netcdf.close(ncid);

% Bottom temperature
ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_TEMP_bottom.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
[ni,nj] = size(TLONG);
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1, [0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
netcdf.close(ncid);

%Bottom detritus
ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_POC_FLUX_IN_bottom.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
[ni,nj] = size(TLONG);
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1, [0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
netcdf.close(ncid);

%spC
ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_spC_150m.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
[ni,nj] = size(TLONG);
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1, [0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
netcdf.close(ncid);

%diatC
ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_diatC_150m.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1,[0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
netcdf.close(ncid);

%diazC
ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_diazC_150m.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    eval([ varname ' = netcdf.getVar(ncid,n-1,[0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
netcdf.close(ncid);

% Zooplankton
ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_zooC_150m.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    %eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs(1)-1],[360 180 length(z200) length(runs)]);
    eval([ varname ' = netcdf.getVar(ncid,n-1,[0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
netcdf.close(ncid);

% ZooLoss
ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_zoo_loss_150m.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
for n = nvars
    varname = netcdf.inqVar(ncid, n-1);
    %eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs(1)-1],[360 180 length(z200) length(runs)]);
    eval([ varname ' = netcdf.getVar(ncid,n-1,[0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
netcdf.close(ncid);

% Doubles and nans
TEMP_150m = double(TEMP_150m);
TEMP_bottom = double(TEMP_bottom);
POC_FLUX_IN_bottom = double(POC_FLUX_IN_bottom);
zooC_150m = double(zooC_150m);
zoo_loss_150m = double(zoo_loss_150m);
diatC_150m = double(diatC_150m);
diazC_150m = double(diazC_150m);
spC_150m = double(spC_150m);

POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;
zooC_150m(zooC_150m >= 9.9e+36) = nan;
zoo_loss_150m(zoo_loss_150m >= 9.9e+36) = nan;


% CALC LARGE FRACTION OF ZOOP FROM ALL PHYTO
fracL = diatC_150m ./ (diatC_150m + spC_150m + diazC_150m);
LzooC_150m = fracL .* zooC_150m;
Lzoo_loss_150m = fracL .* zoo_loss_150m;

clear diatC_150m spC_150m zooC_150m zoo_loss_150m diazC_150m


% ADD TO FOSI CLIMATOL

% ADD DPLE DRIFT-CORR ANOMALIES TO FOSI VALUES
%fosi_Tp fosi_Tb fosi_POC fosi_zooC fosi_loss
TEMP_150m = TEMP_150m + clim_Tp;
TEMP_bottom = TEMP_bottom + clim_Tb;
POC_FLUX_IN_bottom = POC_FLUX_IN_bottom + clim_POC;
LzooC_150m = LzooC_150m + clim_zooC;
Lzoo_loss_150m = Lzoo_loss_150m + clim_loss;

POC_FLUX_IN_bottom(POC_FLUX_IN_bottom<0) = 0.0;
LzooC_150m(LzooC_150m<0) = 0.0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0.0;

%histogram(TEMP_150m)

end
