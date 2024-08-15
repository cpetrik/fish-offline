% Plot max corr of driver-fish corrs
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
%spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',cfile,'/corrs/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';

% Biomass
load([spath,'LMEs_corr_driver_maxcorrs.mat'],'LAtab','LFtab','LPtab','LDtab','lid')

BAtab = LAtab;
BFtab = LFtab;
BPtab = LPtab;
BDtab = LDtab;

clear LAtab LFtab LPtab LDtab

%% Nu
load([spath,'LMEs_nu_corr_driver_maxcorrs.mat'],'LAtab','LFtab','LPtab','LDtab')

PAtab = LAtab;
PFtab = LFtab;
PPtab = LPtab;
PDtab = LDtab;

clear LAtab LFtab LPtab LDtab

%%  ---------------------------------------------------------
cnam = {'coef','p','lag','idriver','driver'};
ctex = {'TP','TB','Det','ZmLoss'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

load('paul_tol_cmaps.mat')

%% Map
%cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
cpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';

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

%% Correlation value map
cmR = cbrewer('seq','Reds',9,'PCHIP');
cmRB = cbrewer('div','RdBu',21,'PCHIP');
cmRB = flipud(cmRB);

% put on grid only significant ones
BARcorr  = nan(ni,nj);
BFRcorr  = nan(ni,nj);
BPRcorr  = nan(ni,nj);
BDRcorr  = nan(ni,nj);
PARcorr  = nan(ni,nj);
PFRcorr  = nan(ni,nj);
PPRcorr  = nan(ni,nj);
PDRcorr  = nan(ni,nj);

for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);

    if (BAtab(i,2) <= 0.05)
        BARcorr(id) = BAtab(i,1);
    end
    if (BFtab(i,2) <= 0.05)
        BFRcorr(id) = BFtab(i,1);
    end
    if (BPtab(i,2) <= 0.05)
        BPRcorr(id) = BPtab(i,1);
    end
    if (BDtab(i,2) <= 0.05)
        BDRcorr(id) = BDtab(i,1);
    end

    if (PAtab(i,2) <= 0.05)
        PARcorr(id) = PAtab(i,1);
    end
    if (PFtab(i,2) <= 0.05)
        PFRcorr(id) = PFtab(i,1);
    end
    if (PPtab(i,2) <= 0.05)
        PPRcorr(id) = PPtab(i,1);
    end
    if (PDtab(i,2) <= 0.05)
        PDRcorr(id) = PDtab(i,1);
    end

end

%% 4plot - Biomass R2 by fn type
f2 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.55 0.45 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,BFRcorr.^2)
hold on
colormap(cmR)
clim([0 1])
title('Forage fishes')
colorbar('Position',[0.25 0.5 0.5 0.025],'Orientation','horizontal')

subplot('Position',[0.5 0.55 0.45 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,BPRcorr.^2)
hold on
colormap(cmR)
clim([0 1])
title('Large pelagics')

subplot('Position',[0.5 0.02 0.45 0.4]) %A
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,BARcorr.^2)
hold on
colormap(cmR)
clim([0 1])
title('All fishes')

subplot('Position',[0.01 0.02 0.45 0.4]) %D
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,BDRcorr.^2)
hold on
colormap(cmR)
clim([0 1])
title('Demersals')

print('-dpng',[ppath 'Map_LMEs_driver_biom_fntypes_corr_R2.png'])

%% 4plot - Production R2 by fn type
f4 = figure('Units','inches','Position',[1 4 7.5 5]);
subplot('Position',[0.01 0.55 0.45 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,PFRcorr.^2)
hold on
colormap(cmR)
clim([0 1])
title('Forage fishes')
colorbar('Position',[0.25 0.5 0.5 0.025],'Orientation','horizontal')

subplot('Position',[0.5 0.55 0.45 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,PPRcorr.^2)
hold on
colormap(cmR)
clim([0 1])
title('Large pelagics')

subplot('Position',[0.5 0.02 0.45 0.4]) %A
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,PARcorr.^2)
hold on
colormap(cmR)
clim([0 1])
title('All fishes')

subplot('Position',[0.01 0.02 0.45 0.4]) %D
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,PDRcorr.^2)
hold on
colormap(cmR)
clim([0 1])
title('Demersals')

print('-dpng',[ppath 'Map_LMEs_driver_nu_fntypes_corr_R2.png'])

%% All R2 of biomass & prod side by side ---------------------------

f1 = figure('Units','inches','Position',[1 5 7.5 2.75]);
subplot('Position',[0.01 0.05 0.425 0.8])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,BARcorr.^2)
colormap(cmR)
clim([0 1])
title('Biomass')

subplot('Position',[0.45 0.05 0.425 0.8])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,PARcorr.^2)
colormap(cmR)
clim([0 1])
title('Production')
colorbar('Position',[0.88 0.1 0.03 0.7])
print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_All_corr_R2.png'])


