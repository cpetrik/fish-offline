% Map clusters of driver-fish mult linear regressions
% For all 63 LMEs
% Includes SST and satellite chl also
% CESM Tp, Tb, Det, ZmLoss
% Fish biomass, nu, gamma

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/' cfile '/corrs/'];

mod = 'v15_All_fish03';

%spath ='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';
load([spath,'LME_biom_nu_gamma_catch_cpue_A_corr_coeffs_cluster_matched.mat'])

%%  Colormap

% colorblind friendly
load('paul_tol_cmaps.mat')

%try muted and add greys
mcol = muted ./ 255;
%add greys 
mcol(10,:) = zmeso(10,:);
% mcol(11,:) = zmeso(9,:);
% mcol(12,:) = zmeso(8,:);

%% Map
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

tlme = double(lme_mask);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
clatlim=[plotminlat plotmaxlat];
clonlim=[plotminlon plotmaxlon];

load coastlines;

%% on grid after matching
Bcorr  = nan(ni,nj);
Pcorr  = nan(ni,nj);
Gcorr  = nan(ni,nj);
Ccorr  = nan(ni,nj);
Ecorr  = nan(ni,nj);

lid = biomA(:,1);
%use col = 2 after finding matching # to text
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Bcorr(id) = ClusterB(i,2);
    Pcorr(id) = ClusterP(i,2);
    Gcorr(id) = ClusterG(i,2);
    Ccorr(id) = ClusterC(i,2);
    Ecorr(id) = ClusterE(i,2);
end


%% Matching
%cme
f3 = figure('Units','inches','Position',[1 3 7.5 5]);

% biom
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Bclus)
colormap(mcol(1:5,:))
clim([1 5])
colorbar('southoutside','Ticks',1:5,'TickLabels',1:5,'Direction','reverse')
title('Biom')

%prod
subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pclus)
colormap(mcol(1:5,:))
clim([1 5])
colorbar('southoutside','Ticks',1:5,'TickLabels',1:5,'Direction','reverse')
title('Prod')

%Gam
subplot('Position',[0.65 0.575 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Gclus)
colormap(mcol(1:5,:))
clim([1 5])
colorbar('southoutside','Ticks',1:5,'TickLabels',1:5,'Direction','reverse')
title('Gamma')

% catch
subplot('Position',[0.01 0.10 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Cclus)
colormap(mcol(1:5,:))
clim([1 5])
colorbar('southoutside','Ticks',1:5,'TickLabels',1:5,'Direction','reverse')
title('Catch')

%cpue
subplot('Position',[0.33 0.10 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Eclus)
colormap(mcol(1:5,:))
clim([1 5])
colorbar('southoutside','Ticks',1:5,'TickLabels',1:5,'Direction','reverse')
title('CPUE')

%print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_gamma_catch_cpue_A_corr_coeffs_cluster_Match.png'])

%% All same 12 colors and text ---------------------------

f1 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Bcorr)
colormap(mcol)
clim([1 10])
title('Biomass')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pcorr)
colormap(mcol)
clim([1 10])
title('Prod')

subplot('Position',[0.65 0.575 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Gcorr)
colormap(mcol)
clim([1 10])
title('Recruitment')

subplot('Position',[0.01 0.10 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Ccorr)
colormap(mcol)
clim([1 10])
title('Catch')

subplot('Position',[0.33 0.10 0.32 0.4]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Ecorr)
colormap(mcol)
clim([1 10])
title('CPUE')
colorbar('Position',[0.675 0.125 0.03 0.35],'Ticks',1:12,'TickLabels',alltex,...
    'Direction','reverse')
%print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_gamma_catch_cpue_A_corr_coeffs_cluster.png'])

