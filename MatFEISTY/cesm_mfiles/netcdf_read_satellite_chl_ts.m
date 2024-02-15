% Create one time series of Seawifs & Modis chl
% Park used, SeaWiFS: SEP1997-DEC2002 and MODIS: JAN2003-DEC2017

clear
close all

cpath = '/Volumes/petrik-lab/Feisty/Obs_data/Chl/';

%% Seawifs ---------------------------------------------------------
ncdisp([cpath 'SEASTAR_SEAWIFS_GAC.19970101_19971231.L3m.YR.CHL.chlor_a.9km.nc'])

%%
% chlor_a
% Size:       4320x2160
% Dimensions: lon,lat
% Datatype:   single
% Attributes:
seawifs_long_name     = 'Chlorophyll Concentration, OCI Algorithm';
seawifs_units         = 'mg m^-3';
seawifs_standard_name = 'mass_concentration_of_chlorophyll_in_sea_water';
seawifs_FillValue    = -32767;
% valid_min     = 0.001;
% valid_max     = 100;
% reference     = 'Hu, C., Lee Z., and Franz, B.A. (2012). Chlorophyll-a algorithms for oligotrophic oceans: A novel approach based on three-band reflectance difference, J. Geophys. Res., 117, C01011, doi:10.1029/2011JC007395.';
% display_scale = 'log';
% display_min   = 0.01;
% display_max   = 20;

%%
tstart = 19970101:10000:20090101;
tend = 19971231:10000:20091231;

%each file has 1 yr 
yrs = 1997:2009;
nyrs = length(tstart);

%% regrid to 1deg grid
lat1 = -89.5:89.5;
lon1 = -179.5:179.5;
[lat_g,lon_g] = meshgrid(lat1,lon1);

[ni,nj] = size(lon_g);
chl = nan(ni,nj,nyrs);

%%
for t=1:length(tstart)

    clear chlor_a chlos ncid

    %%
    ncid = netcdf.open([cpath 'SEASTAR_SEAWIFS_GAC.',num2str(tstart(t)),...
        '_',num2str(tend(t)),'.L3m.YR.CHL.chlor_a.9km.nc'],'NC_NOWRITE');

    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

    % Get all other vars 1st
    for n = 1:(nvars)
        varname = netcdf.inqVar(ncid, n-1);
        eval([ varname ' = netcdf.getVar(ncid,n-1);']);
        eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
    end
    netcdf.close(ncid);

    chlor_a = double(chlor_a);
    chlor_a(chlor_a < 0) = NaN;
    
    %% Use interp2 for data that are already on a regular grid
    [LAT,LON] = meshgrid(lat,lon);
    chlos = interp2(LAT,LON,chlor_a,lat_g,lon_g);

    %% concatenate
    chl(:,:,t) = chlos;

end

%%
save([cpath 'SEASTAR_SEAWIFS_GAC.L3m.YR.CHL.chlor_a_1deg_1997_2009.mat'],'chl',...
    'seawifs_long_name','seawifs_units','seawifs_standard_name',...
    'seawifs_FillValue','lat_g','lon_g','yrs');


%% MODIS ---------------------------------------------------------
clear
close all

cpath = '/Volumes/petrik-lab/Feisty/Obs_data/Chl/';

ncdisp([cpath 'AQUA_MODIS.20020101_20021231.L3m.YR.CHL.chlor_a.9km.nc'])

%%
% chlor_a
% Size:       4320x2160
% Dimensions: lon,lat
% Datatype:   single
% Attributes:
modis_long_name     = 'Chlorophyll Concentration, OCI Algorithm';
modis_units         = 'mg m^-3';
modis_standard_name = 'mass_concentration_of_chlorophyll_in_sea_water';
modis_FillValue    = -32767;
% valid_min     = 0.001;
% valid_max     = 100;
% reference     = 'Hu, C., Lee Z., and Franz, B.A. (2012). Chlorophyll-a algorithms for oligotrophic oceans: A novel approach based on three-band reflectance difference, J. Geophys. Res., 117, C01011, doi:10.1029/2011JC007395.';
% display_scale = 'log';
% display_min   = 0.01;
% display_max   = 20;

%%
tstart = 20020101:10000:20220101;
tend = 20021231:10000:20221231;

%each file has 1 yr 
yrs = 2002:2022;
nyrs = length(tstart);

%% regrid to 1deg grid
lat1 = -89.5:89.5;
lon1 = -179.5:179.5;
[lat_g,lon_g] = meshgrid(lat1,lon1);

[ni,nj] = size(lon_g);
chl = nan(ni,nj,nyrs);

%%
for t=1:length(tstart)

    clear chlor_a chlos ncid

    %%
    ncid = netcdf.open([cpath 'AQUA_MODIS.',num2str(tstart(t)),...
        '_',num2str(tend(t)),'.L3m.YR.CHL.chlor_a.9km.nc'],'NC_NOWRITE');

    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

    % Get all other vars 1st
    for n = 1:(nvars)
        varname = netcdf.inqVar(ncid, n-1);
        eval([ varname ' = netcdf.getVar(ncid,n-1);']);
        eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
    end
    netcdf.close(ncid);

    chlor_a = double(chlor_a);
    chlor_a(chlor_a < 0) = NaN;
    
    %% Use interp2 for data that are already on a regular grid
    [LAT,LON] = meshgrid(lat,lon);
    chlos = interp2(LAT,LON,chlor_a,lat_g,lon_g);

    %% concatenate
    chl(:,:,t) = chlos;

end

%%
myrs = yrs;
save([cpath 'AQUA_MODIS.L3m.YR.CHL.chlor_a_1deg_2002_2022.mat'],'chl',...
    'modis_long_name','modis_units','modis_standard_name',...
    'modis_FillValue','lat_g','lon_g','myrs');


%% Seawifs & MODIS ---------------------------------------------------------
%SeaWiFS: SEP1997-DEC2002 and MODIS: JAN2003-DEC2017

clear chlor_a chlos chl

load([cpath 'SEASTAR_SEAWIFS_GAC.L3m.YR.CHL.chlor_a_1deg_1997_2009.mat']);
syrs = yrs;
schl = chl;
clear yr chl

load([cpath 'AQUA_MODIS.L3m.YR.CHL.chlor_a_1deg_2002_2022.mat'])
mchl = chl;
clear chl

yrs = syrs(1):myrs(end);
sid = find(syrs<2002);

[ni,nj,nt] = size(schl);
chl = nan*ones(ni,nj,length(yrs));
chl(:,:,1:length(sid)) = schl(:,:,sid);
chl(:,:,length(sid)+1:end) = mchl;

save([cpath 'SEAWIFS_1997_2001_MODIS_2002_2022.YR.CHL.chlor_a_1deg.mat'],...
    'seawifs_long_name','seawifs_units','seawifs_standard_name',...
    'modis_long_name','modis_units','modis_standard_name',...
    'lat_g','lon_g','yrs','chl','sid');
