% Read CESM FOSI netcdf

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';

%% All inputs
ncdisp([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.nc'])

%%
% Size:       320x384x816
% Dimensions: nlon,nlat,time
FillValue   = NaN;
missing_value = 9.969209968386869e+36;

% TEMP_150m
TEMP_150m_name  = 'Potential Temperature';
TEMP_150m_units = 'degC';

% zooC_150m
zooC_150m_units = 'mmol/m^3 cm';
zooC_150m_name  = 'thickness of layer k';

% zoo_loss_150m
zoo_loss_150m_units = 'mmol/m^3/s cm';
zoo_loss_150m_name  = 'thickness of layer k';

% diatC_150m
diatC_150m_units = 'mmol/m^3 cm';
diatC_150m_name  = 'thickness of layer k';

% spC_150m
spC_150m_units = 'mmol/m^3 cm';
spC_150m_name  = 'thickness of layer k';

% TEMP_bottom
TEMP_bottom_name  = 'Potential Temperature';
TEMP_bottom_units = 'degC';

% POC_FLUX_IN_bottom
POC_FLUX_IN_bottom_name  = 'POC Flux into Cell';
POC_FLUX_IN_bottom_units = 'mmol/m^3 cm/s'; 

%%
ncid = netcdf.open([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yr = (time-time(1)+1)/365;

%%
save([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat']);





