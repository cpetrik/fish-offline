% Map R2 of driver-fish mult linear regressions
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/' cfile '/corrs/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

%spath ='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';
load([spath,'LME_nu_mlr_coeffs_ALLdiv2SD_alllags3_R2_07_cluster.mat'])

%%  Colormap
cmR=cbrewer('seq','Reds',20,'PCHIP');

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

%%
%R2 on grid
Fbcorr  = nan(ni,nj);
Pbcorr  = nan(ni,nj);
Dbcorr  = nan(ni,nj);
Abcorr  = nan(ni,nj);
Bbcorr  = nan(ni,nj);

Fccorr  = nan(ni,nj);
Pccorr  = nan(ni,nj);
Dccorr  = nan(ni,nj);
Accorr  = nan(ni,nj);

Fncorr  = nan(ni,nj);
Pncorr  = nan(ni,nj);
Dncorr  = nan(ni,nj);
Ancorr  = nan(ni,nj);

lid = Fbiom(:,1);

for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Fbcorr(id) = Fbiom(i,15);
    Pbcorr(id) = Pbiom(i,15);
    Dbcorr(id) = Dbiom(i,15);
    Abcorr(id) = Abiom(i,15);
    Bbcorr(id) = Bbiom(i,15);

    Fccorr(id) = Fcatch(i,15);
    Pccorr(id) = Pcatch(i,15);
    Dccorr(id) = Dcatch(i,15);
    Accorr(id) = Acatch(i,15);
    
    Fncorr(id) = Fnu(i,15);
    Pncorr(id) = Pnu(i,15);
    Dncorr(id) = Dnu(i,15);
    Ancorr(id) = Anu(i,15);
    
end

%% Biomass
f1 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Fbcorr)
colormap(cmR)
caxis([0 1])
title('Forage fishes')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pbcorr)
colormap(cmR)
caxis([0 1])
title('Large pelagics')

subplot('Position',[0.33 0.10 0.32 0.4]) %A
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Abcorr)
colormap(cmR)
caxis([0 1])
title('All fishes')
c=colorbar('Position',[0.675 0.125 0.03 0.35]);
c.Label.String = 'R^2 biomass-drivers';

subplot('Position',[0.01 0.10 0.32 0.4]) %D
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Dbcorr)
colormap(cmR)
caxis([0 1])
title('Demersals')

subplot('Position',[0.65 0.575 0.32 0.4]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Bbcorr)
colormap(cmR)
caxis([0 1])
title('Benthos')
print('-dpng',[ppath 'Map_LMEs_driver_mlr_alllags3_noint_biom_R2_fntypes.png'])


%% Nu
f2 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Fncorr)
colormap(cmR)
caxis([0 1])
title('Forage fishes')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pncorr)
colormap(cmR)
caxis([0 1])
title('Large pelagics')

 %All
subplot('Position',[0.33 0.10 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Ancorr)
colormap(cmR)
caxis([0 1])
title('All fishes')
c=colorbar('Position',[0.675 0.125 0.03 0.35]);
c.Label.String = 'R^2 nu-drivers';

subplot('Position',[0.01 0.10 0.32 0.4]) %D
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Dncorr)
colormap(cmR)
caxis([0 1])
title('Demersals')

% subplot('Position',[0.65 0.575 0.32 0.4]) %B
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
% hold on
% surfm(TLAT,TLONG,Bncorr)
% colormap(cmR)
% caxis([0 1])
% title('Benthos')
print('-dpng',[ppath 'Map_LMEs_driver_mlr_alllags3_noint_nu_R2_fntypes.png'])


%% Obs Catch
f3 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Fccorr)
colormap(cmR)
caxis([0 1])
title('Forage fishes')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pccorr)
colormap(cmR)
caxis([0 1])
title('Large pelagics')

subplot('Position',[0.33 0.10 0.32 0.4]) %A
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Accorr)
colormap(cmR)
caxis([0 1])
title('All fishes')
c=colorbar('Position',[0.675 0.125 0.03 0.35]);
c.Label.String = 'R^2 catch-drivers';

subplot('Position',[0.01 0.10 0.32 0.4]) %D
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Dccorr)
colormap(cmR)
caxis([0 1])
title('Demersals')

% subplot('Position',[0.65 0.575 0.32 0.4]) %B
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
% hold on
% surfm(TLAT,TLONG,Bccorr)
% colormap(cmR)
% caxis([0 1])
%title('Benthos')
print('-dpng',[ppath 'Map_LMEs_driver_mlr_alllags3_noint_FishMIPcatch_R2_fntypes.png'])

