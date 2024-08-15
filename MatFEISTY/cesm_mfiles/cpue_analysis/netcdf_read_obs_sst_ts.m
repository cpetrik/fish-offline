% Read obs SST
% 1971-2000 from OISST
% downloaded 2/16/21

clear
close all

fpath='/Volumes/petrik-lab/Feisty/Obs_data/OISST/';

%% sst
ncdisp([fpath 'sst.day.mean.1982.nc'])

%%
% dimensions:
% lat = 720 ;
% lon = 1440 ;
% time = UNLIMITED ; // (365 currently)
% nbnds = 2 ;

% lat:actual_range = -89.875f, 89.875f ;
% lon:actual_range = 0.125f, 359.875f ;
time_units = "days since 1800-01-01 00:00:0.0" ;
%
% sst(time, lat, lon) ;
sst_units = "degC" ;
sst_valid_range   = [-3  45];
sst_missing_value = -9.96921e+36;
sst_dataset = "NOAA High-resolution Blended Analysis" ;
% sst:actual_range = -1.8, 34.48 ;
sst_long_name = "Daily Long Term Mean Sea Surface Temperature" ;


%%
%each file has 1 yr
yrs = 1982:2020;
nyrs = length(yrs);

%% regrid to 1deg grid
lat1 = -89.5:89.5;
lon1 = -179.5:179.5;
[lat_g,lon_g] = meshgrid(lat1,lon1);

[ni,nj] = size(lon_g);
tos = nan(ni,nj,nyrs);

%%
for t=1:nyrs

    clear sst temp ncid msst msst2

    %%
    ncid = netcdf.open([fpath 'sst.day.mean.',num2str(yrs(t)),...
        '.nc'],'NC_NOWRITE');

    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

    % Get all vars
    if t==1
        for n = 1:(nvars)
            varname = netcdf.inqVar(ncid, n-1);
            eval([ varname ' = netcdf.getVar(ncid,n-1);']);
            eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
        end

        % shift sst lon by 180
        lons = lon;
        lons(lons>180) = lons(lons>180)-360;
        sid = find((lon>180));

        lons2(1:720) = lons(sid);
        lons2(721:1440) = lons(1:720);

        [LAT,LON] = meshgrid(double(lat),double(lons2));

    else
        for n = (nvars)
            varname = netcdf.inqVar(ncid, n-1);
            eval([ varname ' = netcdf.getVar(ncid,n-1);']);
            eval([ varname '(' varname ' == 1.0e+20) = NaN;']);
        end
    end
    netcdf.close(ncid);

    sst = double(sst);
    sst(sst < -9) = NaN;

    %annual mean
    msst = nanmean(sst,3);

    %% Use interp2 for data that are already on a regular grid
    msst2(1:720,:) = msst(sid,:);
    msst2(721:1440,:) = msst(1:720,:);

    temp = interp2(LAT,LON,msst2,lat_g,lon_g);

    %% concatenate
    tos(:,:,t) = temp;

end

%%
save([fpath 'oisst.annual.mean.1deg_1982_2020.mat'],'tos',...
    'sst_units','sst_valid_range','sst_dataset','sst_long_name',...
    'lat_g','lon_g','yrs');

%%
clatlim=[-90 90];
clonlim=[-180 180];
figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,temp)

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(LAT,LON,squeeze(tos(:,:,21)))


