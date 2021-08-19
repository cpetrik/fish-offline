% Read CESM FOSI netcdf

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';

%% All inputs
ncdisp([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.phytoNPP_FIESTY-forcing.nc'])

%%
% Size:       320x384x816
% Dimensions: nlon,nlat,time
FillValue   = NaN;
missing_value = 9.969209968386869e+36;

% photoC_diaz_150m
photoC_diaz_150m_units = 'mmol/m^3/s cm';
photoC_diaz_150m_name  = 'thickness of layer k';

% photoC_diat_150m
photoC_diat_150m_units = 'mmol/m^3/s cm';
photoC_diat_150m_name  = 'thickness of layer k';

% photoC_sp_150m
photoC_sp_150m_units = 'mmol/m^3/s cm';
photoC_sp_150m_name  = 'thickness of layer k';

% time            
time_units = 'days since 0000-01-01 00:00:00';
calendar   = 'noleap';

%%
ncid = netcdf.open([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.phytoNPP_FIESTY-forcing.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yr = (time/365) + 1699;

%%
save([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.phytoNPP_FIESTY-forcing.mat']);





