% Read CESM High Res netcdf Cal Curr only
% all means and depth integrals are over the top 100m 
% (rather than the top 150m). This is because the zoo_loss variable was only 
% output for this hi-res model run as a depth integral over the top 100m, 
% so the other variables are the same to be consistent.

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM_HR/';

%% All inputs
ncdisp([fpath 'g.e22.G1850ECO_JRA_HR.TL319_t13.004.FIESTY-forcing_hiresJRA_CAcurr.nc'])

%%
% Size:       300x370x408
% Dimensions: nlon,nlat,time
FillValue   = NaN;
missing_value = 9.969209968386869e+36;

% TEMP_100m
TEMP_150m_name  = 'Potential Temperature';
TEMP_150m_units = 'degC';

% zooC_100m
zooC_100m_units = 'mmol/m^3 cm';
zooC_100m_name  = 'thickness of layer k';

% zoo_loss_100m
zoo_loss_100m_units = 'mmol/m^3 cm/s';
zoo_loss_100m_name  = 'thickness of layer k';

% diazC_100m
diazC_100m_units = 'mmol/m^3 cm';
diazC_100m_name  = 'thickness of layer k';

% diatC_100m
diatC_100m_units = 'mmol/m^3 cm';
diatC_100m_name  = 'thickness of layer k';

% coccoC_100m
coccoC_100m_units = 'mmol/m^3 cm';
coccoC_100m_name  = 'thickness of layer k';

% spC_100m
spC_100m_units = 'mmol/m^3 cm';
spC_100m_name  = 'thickness of layer k';

% TEMP_bottom
TEMP_bottom_name  = 'Potential Temperature';
TEMP_bottom_units = 'degC';

% POC_FLUX_IN_bottom
POC_FLUX_IN_bottom_name  = 'POC Flux into Cell';
POC_FLUX_IN_bottom_units = 'mmol/m^3 cm/s'; 

%%
ncid = netcdf.open([fpath 'g.e22.G1850ECO_JRA_HR.TL319_t13.004.FIESTY-forcing_hiresJRA_CAcurr.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
%34 years, 1958 to 1991
time = 1:408;
yr = time/12 + 1958;
%yr = 1958+(1/12):(1/12):1991;

%%
save([fpath 'g.e22.G1850ECO_JRA_HR.TL319_t13.004.FIESTY-forcing_hiresJRA_CAcurr.mat']);





