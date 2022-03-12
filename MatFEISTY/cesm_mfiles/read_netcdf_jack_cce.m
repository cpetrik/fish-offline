% Read CESM DPLE netcdf for CCE from Jack

clear all
close all

fpath='/Users/cpetrik/Dropbox/TAMU/iHESP/NSF_convergence/modeling/Jack/';

%% All inputs
ncdisp([fpath 'CCLME_HMXL_FINAL.nc'])

%%
ncid = netcdf.open([fpath 'CCLME_HMXL_FINAL.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end

netcdf.close(ncid);

HMXL(HMXL>1e30) = nan;
HMXL =double(HMXL);

%%
test = squeeze(HMXL(:,:,1,1));

[ni,nj,nt,nm] = size(HMXL);
t2=reshape(HMXL,ni*nj,nt,nm);
mts = squeeze(nanmean(t2,1));

%%
[lat,lon] = meshgrid(TLAT,TLONG); 
%TLAT:18.9-46.9; TLONG:223.8-251.9 = off
%TLAT:21.9-48.4; TLONG:227.3-252.8 = should be

plotminlat=18; %Set these bounds for your data
plotmaxlat=50;
plotminlon=220;
plotmaxlon=255;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;     

%%
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat,lon,test)
% pcolor(lat,lon,test)
% shading flat
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
%caxis([0 1.1]);
set(gcf,'renderer','painters')


figure(2)
plot(1:122,mts)

%%
figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat2,lon2,AllF(:,:,1))
% pcolor(lat,lon,test)
% shading flat
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
colorbar
%caxis([0 1.1]);
set(gcf,'renderer','painters')