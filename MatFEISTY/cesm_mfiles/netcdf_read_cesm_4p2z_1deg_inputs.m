% Read CESM 1 degree companion run for JRA55-forced high res
% Early cycle


clear all
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/4P2Z/';
spath = '/Volumes/petrik-lab/Feisty/test/';

%% All inputs
ncdisp([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001.marbl-scope-FIESTY-forcing.nc'])

%%
% Size:       300x370x408
% Dimensions: nlon,nlat,time
FillValue   = NaN;
missing_value = 9.969209968386869e+36;

% TEMP_150m
TEMP_150m_name  = 'Potential Temperature';
TEMP_150m_units = 'degC';

% mesozooC_150m
mesozooC_150m_units = 'mmol/m^3 cm';
mesozooC_150m_name  = 'Mesozooplankton Biomass Vertical Integral';

% mesozoo_loss_150m
mesozoo_loss_150m_units = 'mmol/m^3 cm/s';
mesozoo_loss_150m_name  = 'Mesozooplankton Loss Vertical Integral';

% TEMP_bottom
TEMP_bottom_name  = 'Potential Temperature';
TEMP_bottom_units = 'degC';

% POC_FLUX_IN_bottom
POC_FLUX_IN_bottom_name  = 'POC Flux into Cell';
POC_FLUX_IN_bottom_units = 'mmol/m^3 cm/s'; 

%%
ncid = netcdf.open([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001.marbl-scope-FIESTY-forcing.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
%61 years, 1957? to 2018?
yr = (time/365)+ 1712;

%%
save([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001.marbl-scope-FIESTY-forcing.mat']);
save([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001.grid.mat'],'TLAT','TLONG','TAREA');






