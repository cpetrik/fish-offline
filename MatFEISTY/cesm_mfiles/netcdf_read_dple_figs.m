% Read CESM 1 degree companion run for JRA55-forced high res
% Later cycle


clear 
close all

fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/DPLE_ms/';

%% Fig 1
fold1 = 'Fig.1/';
file1 = 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.TEMP_mean_150m.19580101-20211231.nc';
ncdisp([fpath file1])

%% Fig 2
file2 = 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.TEMP_BOTTOM_2.19580101-20211231.nc';
ncdisp([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.TEMP_BOTTOM_2.19580101-20211231.nc'])

%% Fig 3
file3 = 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.pocToFloor_2.19580101-20211231.nc';
ncdisp([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.pocToFloor_2.19580101-20211231.nc'])

%% Fig 4
file4 = 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.mesozooC_zint_150m.19580101-20211231.nc';
ncdisp([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.mesozooC_zint_150m.19580101-20211231.nc'])

%% Figs 5-8
file5 = 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.19580101-20211231.nc';
ncdisp([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.19580101-20211231.nc'])

%% Figs S5
fileS5 = 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.19580101-20211231.nc';
ncdisp([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.19580101-20211231.nc'])

%% Figs S7
fileS7 = 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.19580101-20211231.nc';
ncdisp([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.19580101-20211231.nc'])

%% Figs S14-15
fileS14 = 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.19580101-20211231.nc';
ncdisp([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.19580101-20211231.nc'])

%% Figs S19
fileS19 = 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.19580101-20211231.nc';
ncdisp([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.19580101-20211231.nc'])

%%
dinfo = dir([fpath fold1]);
full_filenames = fullfile({dinfo.folder}, {dinfo.name});



%%
% Size:       320x384x23360
% Dimensions: nlon,nlat,time
FillValue   = NaN;
missing_value = 9.969209968386869e+36;

% TEMP_150m var=18
TEMP_150m_name  = 'Potential Temperature 0-150m Vertical Mean';
TEMP_150m_units = 'degC';

% mesozooC_150m var=36
mesozooC_150m_units = 'mmol/m^3 cm';
mesozooC_150m_name  = 'Mesozooplankton Carbon 0-150m Vertical Integral';

% mesozoo_loss_150m var=36
mesozoo_loss_150m_units = 'mmol/m^3 cm/s';
mesozoo_loss_150m_name  = 'Mesozooplankton Loss Vertical Integral, 0-150m';

% TEMP_bottom var=18
TEMP_bottom_name  = 'Potential temperature Value at Sea Floor';
TEMP_bottom_units = 'degC';

% POC_FLUX_IN_bottom var=41
POC_FLUX_IN_bottom_name  = 'POC Flux Hitting Sea Floor';
POC_FLUX_IN_bottom_units = 'nmol/cm^2/s';      %'mmol/m^3 cm/s'; -->'nmol/cm^2/s'

%% TP
ncid = netcdf.open([fpath file1],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 18 %1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% TB
ncid = netcdf.open([fpath file2],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 18 %1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Det
ncid = netcdf.open([fpath file3],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 41 %1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Zmeso
ncid = netcdf.open([fpath file4],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 36 %1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Zloss
ncid = netcdf.open([fpath file5],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 36 %1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yrs = 1958:2021;
nyrs = 23360/365;

%% Just save 1st year to test
ptemp = TEMP_mean_150m(:,:,1:365);
btemp = TEMP_BOTTOM_2(:,:,1:365);
det = pocToFloor_2(:,:,1:365);
zmeso = mesozooC_zint_150m(:,:,1:365);
zloss = mesozoo_loss_zint_150m(:,:,1:365);

%%
save([spath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001branch.pop.h.ecosys.nday1.FEISTY_inputs.1958.mat'],...
    'ptemp','btemp','det','zmeso','zloss','TEMP_150m_name','TEMP_150m_units','mesozooC_150m_units','mesozooC_150m_name',...
    'mesozoo_loss_150m_units','mesozoo_loss_150m_name','TEMP_bottom_name','TEMP_bottom_units','POC_FLUX_IN_bottom_name',...
    'POC_FLUX_IN_bottom_units','FillValue');





