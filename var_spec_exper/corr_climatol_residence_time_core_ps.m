% Corr residence time with spectral slope
% Residence time calc from Climatol
% Spectral slope from CORE-forced
% On different grids and forcing was not identical

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

%% Climatol
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
ppath = [pp cfile '/Climatol/'];
% Grid
load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_gridspec.mat');
% residence times
load([fpath 'Residence_time_means_Climatol_' harv '_' cfile '.mat'])

%% CORE
spath='/Volumes/MIP/GCM_DATA/CORE-forced/';
% power spectra
load([fpath 'powerspec_feisty_core_1950_2007_ln.mat'])
% Grid
load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

%% Put on each grid first
%Climatol
[ni,nj]=size(lon);
Csf=NaN*ones(ni,nj);
Csp=NaN*ones(ni,nj);
Csd=NaN*ones(ni,nj);
Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);
Csf(ID)=log10(sf_mres2);
Csp(ID)=log10(sp_mres2);
Csd(ID)=log10(sd_mres2);
Cmf(ID)=log10(mf_mres2);
Cmp(ID)=log10(mp_mres2);
Cmd(ID)=log10(md_mres2);
Clp(ID)=log10(lp_mres2);
Cld(ID)=log10(ld_mres2);
% Csf(ID)=1./(sf_mres2);
% Csp(ID)=1./(sp_mres2);
% Csd(ID)=1./(sd_mres2);
% Cmf(ID)=1./(mf_mres2);
% Cmp(ID)=1./(mp_mres2);
% Cmd(ID)=1./(md_mres2);
% Clp(ID)=1./(lp_mres2);
% Cld(ID)=1./(ld_mres2);

%Hist (CORE)
Hsf=sf_ps;
Hsp=sp_ps;
Hsd=sd_ps;
Hmf=mf_ps;
Hmp=mp_ps;
Hmd=md_ps;
Hlp=lp_ps;
Hld=ld_ps;


%% Interpolate to same grid
%lat        [-89.5 89.5]
%geolat_t   [-81.5 89.4879]
%lon        [0.5 359.5]
%geolon_t   [-279.9803 79.9803]

% Need to fix both longitudes
test = lon-360;
id=find(test<-180);
test(id)=test(id)+360;
lon = test;

test2=geolon_t;
id=find(test2<-180);
test2(id)=test2(id)+360;
geolon_t = test2;

geolat = double(geolat_t');
geolon = double(geolon_t');

lats = -89.5:89.5;
lons = -179.5:179.5;
[glon,glat] = meshgrid(lons,lats);

clat = lat;
clon = lon;

hSF = griddata(geolat,geolon,Hsf',glat,glon);
hSP = griddata(geolat,geolon,Hsp',glat,glon);
hSD = griddata(geolat,geolon,Hsd',glat,glon);
hMF = griddata(geolat,geolon,Hmf',glat,glon);
hMP = griddata(geolat,geolon,Hmp',glat,glon);
hMD = griddata(geolat,geolon,Hmd',glat,glon);
hLP = griddata(geolat,geolon,Hlp',glat,glon);
hLD = griddata(geolat,geolon,Hld',glat,glon);

cSF = griddata(clat,clon,Csf,glat,glon);
cSP = griddata(clat,clon,Csp,glat,glon);
cSD = griddata(clat,clon,Csd,glat,glon);
cMF = griddata(clat,clon,Cmf,glat,glon);
cMP = griddata(clat,clon,Cmp,glat,glon);
cMD = griddata(clat,clon,Cmd,glat,glon);
cLP = griddata(clat,clon,Clp,glat,glon);
cLD = griddata(clat,clon,Cld,glat,glon);

%% Correlations
%r
r(1) = corr(hSF(:),cSF(:),'rows','complete');
r(2) = corr(hSP(:),cSP(:),'rows','complete');
r(3) = corr(hSD(:),cSD(:),'rows','complete');
r(4) = corr(hMF(:),cMF(:),'rows','complete');
r(5) = corr(hMP(:),cMP(:),'rows','complete');
r(6) = corr(hMD(:),cMD(:),'rows','complete');
r(7) = corr(hLP(:),cLP(:),'rows','complete');
r(8) = corr(hLD(:),cLD(:),'rows','complete');

%% Figure of scatter plots
figure(1)
subplot(3,3,1)
%scatter(cSF(:),hSF(:),20,ptemp(:),'filled'); hold on;
%cmocean('thermal');
%colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
scatter(cSF(:),hSF(:),10,'k','filled'); hold on;
text(0.5,-1.5,['r = ' sprintf('%2.2f',r(1))],'color','r')
axis([-3 3 -6 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('SF')

subplot(3,3,2)
scatter(cSP(:),hSP(:),10,'k','filled'); hold on;
text(-20,-1.5,['r = ' sprintf('%2.2f',r(2))],'color','r')
axis([-50 3 -6 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('SP')

subplot(3,3,3)
scatter(cSD(:),hSD(:),10,'k','filled'); hold on;
text(0.5,-1.5,['r = ' sprintf('%2.2f',r(3))],'color','r')
axis([-4 3 -6 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('SD')

subplot(3,3,4)
scatter(cMF(:),hMF(:),10,'k','filled'); hold on;
text(0,-1.5,['r = ' sprintf('%2.2f',r(4))],'color','r')
axis([-3 3 -6 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('MF')

subplot(3,3,5)
scatter(cMP(:),hMP(:),10,'k','filled'); hold on;
text(-20,-1.5,['r = ' sprintf('%2.2f',r(5))],'color','r')
axis([-50 3 -6 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('MP')

subplot(3,3,6)
scatter(cMD(:),hMD(:),10,'k','filled'); hold on;
text(2,-1.5,['r = ' sprintf('%2.2f',r(6))],'color','r')
axis([0.5 3 -7 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('MD')

subplot(3,3,8)
scatter(cLP(:),hLP(:),10,'k','filled'); hold on;
text(-10,-1.5,['r = ' sprintf('%2.2f',r(7))],'color','r')
axis([-30 3 -6 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('LP')

subplot(3,3,9)
scatter(cLD(:),hLD(:),10,'k','filled'); hold on;
text(2,-1.5,['r = ' sprintf('%2.2f',r(8))],'color','r')
axis([1.5 3.5 -7 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('LD')
stamp('')
print('-dpng',[ppath 'Corr_Climatol_log10residence_time_CORE_reddening.png'])

%%
figure(2)
subplot(3,3,1)
%scatter(cSF(:),hSF(:),20,ptemp(:),'filled'); hold on;
%cmocean('thermal');
%colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
scatter(cSF(:),hSF(:),10,'k','filled'); hold on;
%text(0.5,-1.5,['r = ' sprintf('%2.2f',r(1))],'color','r')
axis([0 300 -6 -1])
xlabel('1/ Residence time')
ylabel('Spectral slope')
title('SF')

subplot(3,3,2)
scatter(cSP(:),hSP(:),10,'k','filled'); hold on;
%text(-20,-1.5,['r = ' sprintf('%2.2f',r(2))],'color','r')
axis([0 4e3 -6 -1])
xlabel('1/ Residence time')
ylabel('Spectral slope')
title('SP')

subplot(3,3,3)
scatter(cSD(:),hSD(:),10,'k','filled'); hold on;
%text(0.5,-1.5,['r = ' sprintf('%2.2f',r(3))],'color','r')
axis([0 4e3 -6 -1])
xlabel('1/ Residence time')
ylabel('Spectral slope')
title('SD')

subplot(3,3,4)
scatter(cMF(:),hMF(:),10,'k','filled'); hold on;
%text(0,-1.5,['r = ' sprintf('%2.2f',r(4))],'color','r')
axis([0 0.1 -6 -1])
xlabel('1/ Residence time')
ylabel('Spectral slope')
title('MF')

subplot(3,3,5)
scatter(cMP(:),hMP(:),10,'k','filled'); hold on;
%text(-20,-1.5,['r = ' sprintf('%2.2f',r(5))],'color','r')
axis([0 1 -6 -1])
xlabel('1/ Residence time')
ylabel('Spectral slope')
title('MP')

subplot(3,3,6)
scatter(cMD(:),hMD(:),10,'k','filled'); hold on;
%text(2,-1.5,['r = ' sprintf('%2.2f',r(6))],'color','r')
axis([0 0.15 -7 -1])
xlabel('1/ Residence time')
ylabel('Spectral slope')
title('MD')

subplot(3,3,8)
scatter(cLP(:),hLP(:),10,'k','filled'); hold on;
%text(-10,-1.5,['r = ' sprintf('%2.2f',r(7))],'color','r')
axis([0 4e3 -6 -1])
xlabel('1/ Residence time')
ylabel('Spectral slope')
title('LP')

subplot(3,3,9)
scatter(cLD(:),hLD(:),10,'k','filled'); hold on;
text(2,-1.5,['r = ' sprintf('%2.2f',r(8))],'color','r')
axis([0 0.015 -7 -1])
xlabel('1/ Residence time')
ylabel('Spectral slope')
title('LD')
stamp('')
%print('-dpng',[ppath 'Corr_Climatol_inv_residence_time_CORE_reddening.png'])
