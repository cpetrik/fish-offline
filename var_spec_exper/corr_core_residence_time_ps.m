% Corr residence time with spectral slope
% Residence time calc from CORE
% Spectral slope from CORE-forced
% On different grids and forcing was not identical

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/',...
    cfile,'/CORE/'];

%% CORE
spath='/Volumes/MIP/GCM_DATA/CORE-forced/';
% power spectra
load([fpath 'powerspec_feisty_core_1950_2007_ln.mat'])
% residence times
load([fpath 'Residence_time_means_Core_' harv '_' cfile '.mat'])

% Grid
load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

%% Put on each grid first
[ni,nj]=size(geolon_t);
ID = GRD.ID;

%CORE residence times
cSF=NaN*ones(ni,nj);
cSP=NaN*ones(ni,nj);
cSD=NaN*ones(ni,nj);
cMF=NaN*ones(ni,nj);
cMP=NaN*ones(ni,nj);
cMD=NaN*ones(ni,nj);
cLP=NaN*ones(ni,nj);
cLD=NaN*ones(ni,nj);
cSF(ID)=log10(sf_mres2);
cSP(ID)=log10(sp_mres2);
cSD(ID)=log10(sd_mres2);
cMF(ID)=log10(mf_mres2);
cMP(ID)=log10(mp_mres2);
cMD(ID)=log10(md_mres2);
cLP(ID)=log10(lp_mres2);
cLD(ID)=log10(ld_mres2);
% Csf(ID)=1./(sf_mres2);
% Csp(ID)=1./(sp_mres2);
% Csd(ID)=1./(sd_mres2);
% Cmf(ID)=1./(mf_mres2);
% Cmp(ID)=1./(mp_mres2);
% Cmd(ID)=1./(md_mres2);
% Clp(ID)=1./(lp_mres2);
% Cld(ID)=1./(ld_mres2);

%CORE powerspec slopes
hSF=sf_ps;
hSP=sp_ps;
hSD=sd_ps;
hMF=mf_ps;
hMP=mp_ps;
hMD=md_ps;
hLP=lp_ps;
hLD=ld_ps;

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
text(3,-1.5,['r = ' sprintf('%2.2f',r(1))],'color','r')
axis([0 4 -6 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('SF')

subplot(3,3,2)
scatter(cSP(:),hSP(:),10,'k','filled'); hold on;
text(3,-1.5,['r = ' sprintf('%2.2f',r(2))],'color','r')
axis([0 4 -6 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('SP')

subplot(3,3,3)
scatter(cSD(:),hSD(:),10,'k','filled'); hold on;
text(3,-1.5,['r = ' sprintf('%2.2f',r(3))],'color','r')
axis([0 4 -6 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('SD')

subplot(3,3,4)
scatter(cMF(:),hMF(:),10,'k','filled'); hold on;
text(2.5,-1.5,['r = ' sprintf('%2.2f',r(4))],'color','r')
axis([1 3 -6 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('MF')

subplot(3,3,5)
scatter(cMP(:),hMP(:),10,'k','filled'); hold on;
text(2.5,-1.5,['r = ' sprintf('%2.2f',r(5))],'color','r')
axis([1 3 -6 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('MP')

subplot(3,3,6)
scatter(cMD(:),hMD(:),10,'k','filled'); hold on;
text(2.5,-1.5,['r = ' sprintf('%2.2f',r(6))],'color','r')
axis([1 3 -7 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('MD')

subplot(3,3,8)
scatter(cLP(:),hLP(:),10,'k','filled'); hold on;
text(2.8,-1.5,['r = ' sprintf('%2.2f',r(7))],'color','r')
axis([2.5 3 -6 -1])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('LP')

subplot(3,3,9)
scatter(cLD(:),hLD(:),10,'k','filled'); hold on;
text(2.8,-2.5,['r = ' sprintf('%2.2f',r(8))],'color','r')
axis([2.5 3 -7 -2])
xlabel('log_1_0 Residence time')
ylabel('Spectral slope')
title('LD')
stamp('')
print('-dpng',[ppath 'Corr_log10residence_time_CORE_reddening.png'])

