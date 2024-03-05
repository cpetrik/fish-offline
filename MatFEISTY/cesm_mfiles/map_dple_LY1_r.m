% Plot r value of DPLE pred pot LY1
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/'];

mod = 'v15_All_fish03';


load('paul_tol_cmaps.mat')
mcol(1,:) = drainbow(12,:)/255; % orange
mcol(2,:) = drainbow(4,:)/255; %dk blue
mcol(3,:) = drainbow(15,:)/255; %grey 
mcol(4,:) = drainbow(6,:)/255; %lt blue

%%
Ftab = nan*ones(66,1);
Ptab = nan*ones(66,1);
Dtab = nan*ones(66,1);

Ftab(3) = 0.47;
Ftab(13) = 0.78;
Ftab(25) = 0.67;
Ftab(27) = 0.78;
Ftab(29) = 0.75;

Ptab(13) = 0.66;
Ptab(25) = 0.60;
Ptab(49) = 0.53;

Dtab(1) = 0.49;
Dtab(8) = 0.58;
Dtab(18) = 0.70;
Dtab(21) = 0.73;
Dtab(39) = 0.49;

%% Map
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
Fcorr  = nan(ni,nj);
Pcorr  = nan(ni,nj);
Dcorr  = nan(ni,nj);
for i=1:66
    id = find(tlme==i);
    Fcorr(id) = Ftab(i);
    Pcorr(id) = Ptab(i);
    Dcorr(id) = Dtab(i);
end

%%
% f2 = figure('Units','inches','Position',[1 3 7.5 2.75]);
% subplot('Position',[0.01 0.05 0.45 0.8])
figure(1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Fcorr)
hold on
colormap(cmR)
clim([0 0.9])
colorbar('SouthOutside')
title('Forage Fish')
print('-dpng',[ppath 'Map_LMEs_DPLE_LY1_corr_F.png'])

%subplot('Position',[0.47 0.05 0.45 0.8])
figure(2)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Pcorr)
hold on
colormap(cmR)
clim([0 0.9])
colorbar('SouthOutside')
title('Large Pelagics')
print('-dpng',[ppath 'Map_LMEs_DPLE_LY1_corr_P.png'])

figure(3)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Dcorr)
hold on
colormap(cmR)
clim([0 0.9])
colorbar('SouthOutside')
title('Demersals')
print('-dpng',[ppath 'Map_LMEs_DPLE_LY1_corr_D.png'])


