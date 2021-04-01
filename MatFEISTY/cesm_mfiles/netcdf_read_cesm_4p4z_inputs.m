% Read CESM 4P4Z netcdf

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/4P4Z/';

%% All inputs
ncdisp([fpath 'g.e22a06.G1850ECOIAF_JRA_PHYS_DEV.TL319_g17.4p4z.004.FIESTY-forcing.nc'])

%%
% source: 'CCSM POP2, the CCSM Ocean Component'
% Size:       320x384x732
% Dimensions: nlon,nlat,time
time_units = 'days since 0000-01-01 00:00:00';
FillValue   = NaN;
missing_value = 9.969209968386869e+36;

% TEMP_150m
TEMP_150m_name  = 'Potential Temperature';
TEMP_150m_units = 'degC';

% zoo3C_150m
zoo3C_150m_units = 'mmol/m^3 cm';
zoo3C_150m_name  = 'thickness of layer k';

% zoo4C_150m
zoo4C_150m_units = 'mmol/m^3 cm';
zoo4C_150m_name  = 'thickness of layer k';

% zoo_loss_150m
zoo3_loss_150m_units = 'mmol/m^3/s cm';
zoo3_loss_150m_name  = 'thickness of layer k';

% zoo_loss_150m
zoo4_loss_150m_units = 'mmol/m^3/s cm';
zoo4_loss_150m_name  = 'thickness of layer k';

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
ncid = netcdf.open([fpath 'g.e22a06.G1850ECOIAF_JRA_PHYS_DEV.TL319_g17.4p4z.004.FIESTY-forcing.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
netcdf.close(ncid);

%% Time
%yr = (time-time(1)+1)/365;
%yr = (time-time(1))/365;
yr = (time)/365;

%%
save([fpath 'g.e22a06.G1850ECOIAF_JRA_PHYS_DEV.TL319_g17.4p4z.004.FIESTY-forcing.mat']);





