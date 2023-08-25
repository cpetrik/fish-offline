% Map clusters of driver-fish mult linear regressions
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
load([spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_alllags_reduc.mat'])

%% Cluster descrips
ftex = {'+TP',...             = 1
    '+ZmLoss',...             = 2
    '+Zmeso',...              = 3
    '+Zmeso, -ZmLoss',...     = 4
    '+TP*Zmeso, -TP*ZmLoss'};% = 5
    
ptex = {'+TP',...               = 1
    '+Zmeso, -ZmLoss',...       = 2
    '-Zmeso, +ZmLoss',...       = 3
    '-TP, +Zmeso',...           = 4
    '-TP'};                %    = 5
    
dtex = {'+Det',...                           = 1
    '-Zmeso, +ZmLoss',...                    = 2
    '+Zmeso, -ZmLoss',...                    = 3
    '-Det?',...                              = 4
    '-Zm, +ZmL, -TP*Zm, +TP*ZmL'}; %        = 5

atex = {'+TP',...                 = 1
    '-Zmeso, +ZmLoss, +TP',...    = 2
    '+Zmeso',...                  = 3
    '+Zmeso, -ZmLoss',...         = 4
    '+Zmeso, -ZmLoss',...         = 5
    '-Zmeso, +ZmLoss, -TP'}; %    = 6

%% Add a text descript vector that matches to clusters
CatF = cell(length(ClusterF),1);
CatF(ClusterF==1) = ftex(1);
CatF(ClusterF==2) = ftex(2);
CatF(ClusterF==3) = ftex(3);
CatF(ClusterF==4) = ftex(4);
CatF(ClusterF==5) = ftex(5);

CatP = cell(length(ClusterP),1);
CatP(ClusterP==1) = ptex(1);
CatP(ClusterP==2) = ptex(2);
CatP(ClusterP==3) = ptex(3);
CatP(ClusterP==4) = ptex(4);
CatP(ClusterP==5) = ptex(5);

CatD = cell(length(ClusterD),1);
CatD(ClusterD==1) = dtex(1);
CatD(ClusterD==2) = dtex(2);
CatD(ClusterD==3) = dtex(3);
CatD(ClusterD==4) = dtex(4);
CatD(ClusterD==5) = dtex(5);

CatA = cell(length(ClusterA),1);
CatA(ClusterA==1) = atex(1);
CatA(ClusterA==2) = atex(2);
CatA(ClusterA==3) = atex(3);
CatA(ClusterA==4) = atex(4);
CatA(ClusterA==5) = atex(5);
CatA(ClusterA==6) = atex(6);

CatF = string(CatF);
CatP = string(CatP);
CatD = string(CatD);
CatA = string(CatA);


%% Create one master colormap for all categories

% Biomass clusters (use same colors if same)
% alltex = {'+Det',...                                = 1
%     '+Zmeso, -ZmLoss',...                            = 2
%     '-Zmeso, +ZmLoss',...                            = 3
%     '-TP*Zmeso, +TP*ZmLoss',...                      = 4
%     '+Det, -Zmeso, +ZmLoss',...                     = 5
%     '-Det, -Zmeso, +ZmLoss',...                     = 6
%     '-Zmeso, +ZmLoss, -TP*Zmeso, +TP*ZmLoss',...    = 7
%     '-Zmeso, +ZmLoss, +TP*Zmeso, -TP*ZmLoss',...    = 8
%     '+Zmeso, +ZmLoss',...                           = 9
%     'unk1',...                                      = 10
%     'unk2',...                                      = 11
%     'unk3'};                                       %= 12

alltex = {'+Zmeso',...
    '+Zmeso, -ZmLoss',...           %biom=2 [0.5333    0.8000    0.9333]
    '-Zmeso, +ZmLoss',...           %biom=3 [0.2667    0.6667    0.6000]
    '-Zmeso, +ZmLoss, +TP',...
    '-Zmeso, +ZmLoss, -TP',...
    '-Zm, +ZmL, -TP*Zm, +TP*ZmL',...%biom=7 [0.8000    0.4000    0.4667]
    '+ZmLoss',...
    '+TP*Zmeso, -TP*ZmLoss',...
    '+TP',...
    '-TP',...
    '-TP, +Zmeso',...
    '+Det',...                      %biom=1 [0.2000    0.1333    0.5333]
    '-Det?'};
    

ClusterF(:,2) = nan;
ClusterF((CatF==alltex(1)),2) = 1;
ClusterF((CatF==alltex(2)),2) = 2;
ClusterF((CatF==alltex(3)),2) = 3;
ClusterF((CatF==alltex(4)),2) = 4;
ClusterF((CatF==alltex(5)),2) = 5;
ClusterF((CatF==alltex(6)),2) = 6;
ClusterF((CatF==alltex(7)),2) = 7;
ClusterF((CatF==alltex(8)),2) = 8;
ClusterF((CatF==alltex(9)),2) = 9;
ClusterF((CatF==alltex(10)),2) = 10;
ClusterF((CatF==alltex(11)),2) = 11;
ClusterF((CatF==alltex(12)),2) = 12;
ClusterF((CatF==alltex(13)),2) = 13;

ClusterP(:,2) = nan;
ClusterP((CatP==alltex(1)),2) = 1;
ClusterP((CatP==alltex(2)),2) = 2;
ClusterP((CatP==alltex(3)),2) = 3;
ClusterP((CatP==alltex(4)),2) = 4;
ClusterP((CatP==alltex(5)),2) = 5;
ClusterP((CatP==alltex(6)),2) = 6;
ClusterP((CatP==alltex(7)),2) = 7;
ClusterP((CatP==alltex(8)),2) = 8;
ClusterP((CatP==alltex(9)),2) = 9;
ClusterP((CatP==alltex(10)),2) = 10;
ClusterP((CatP==alltex(11)),2) = 11;
ClusterP((CatP==alltex(12)),2) = 12;
ClusterP((CatP==alltex(13)),2) = 13;

ClusterD(:,2) = nan;
ClusterD((CatD==alltex(1)),2) = 1;
ClusterD((CatD==alltex(2)),2) = 2;
ClusterD((CatD==alltex(3)),2) = 3;
ClusterD((CatD==alltex(4)),2) = 4;
ClusterD((CatD==alltex(5)),2) = 5;
ClusterD((CatD==alltex(6)),2) = 6;
ClusterD((CatD==alltex(7)),2) = 7;
ClusterD((CatD==alltex(8)),2) = 8;
ClusterD((CatD==alltex(9)),2) = 9;
ClusterD((CatD==alltex(10)),2) = 10;
ClusterD((CatD==alltex(11)),2) = 11;
ClusterD((CatD==alltex(12)),2) = 12;
ClusterD((CatD==alltex(13)),2) = 13;

ClusterA(:,2) = nan;
ClusterA((CatA==alltex(1)),2) = 1;
ClusterA((CatA==alltex(2)),2) = 2;
ClusterA((CatA==alltex(3)),2) = 3;
ClusterA((CatA==alltex(4)),2) = 4;
ClusterA((CatA==alltex(5)),2) = 5;
ClusterA((CatA==alltex(6)),2) = 6;
ClusterA((CatA==alltex(7)),2) = 7;
ClusterA((CatA==alltex(8)),2) = 8;
ClusterA((CatA==alltex(9)),2) = 9;
ClusterA((CatA==alltex(10)),2) = 10;
ClusterA((CatA==alltex(11)),2) = 11;
ClusterA((CatA==alltex(12)),2) = 12;
ClusterA((CatA==alltex(13)),2) = 13;


%%  Colormap

% colorblind friendly
load('paul_tol_cmaps.mat')

%try night w/ added biomass matching colors
mcol = night ./ 255;
cbio = muted ./ 255;

mcol2(1,:) = mcol(7,:); %5?
mcol2(2,:) = [0.5333    0.8000    0.9333];
mcol2(3,:) = [0.2667    0.6667    0.6000];
mcol2(6,:) = [0.8000    0.4000    0.4667];
mcol2(12,:) = [0.2000    0.1333    0.5333];

mcol2(4,:) = mcol(3,:);
mcol2(5,:) = mcol(1,:);
mcol2(7:11,:) = mcol(8:2:16,:);
mcol2(13,:) = drainbow(4,:)./ 255;


%% Shorter text for colorbar labels
% alltex2 = {'+Det',...                                = 1
%     '+Zmeso, -ZmLoss',...                            = 2
%     '-Zmeso, +ZmLoss',...                            = 3
%     '-TP*Z, +TP*ZL',...                      = 4
%     '+Det, -Zmeso, +ZmLoss',...                     = 5
%     '-Det, -Zmeso, +ZmLoss',...                     = 6
%     '-Zmeso, +ZmLoss, -TP*Z, +TP*ZL',...    = 7
%     '-Zmeso, +ZmLoss, +TP*Z, -TP*ZL',...    = 8
%     '+Zmeso, +ZmLoss',...                           = 9
%     'unk1',...                                      = 10
%     'unk2',...                                      = 11
%     'unk3'};                                       %= 12

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
% on grid
Fcorr  = nan(ni,nj);
Pcorr  = nan(ni,nj);
Dcorr  = nan(ni,nj);
Acorr  = nan(ni,nj);

for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Fcorr(id) = ClusterF(i,2);
    Pcorr(id) = ClusterP(i,2);
    Dcorr(id) = ClusterD(i,2);
    Acorr(id) = ClusterA(i,2);
end

%% All
figure(10)
subplot('Position',[0.01 0.13 0.725 0.725])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Acorr)
colormap(mcol2)
caxis([1 13])
colorbar('Position',[0.71 0.25 0.03 0.5],'Ticks',1:13,'TickLabels',alltex,...
    'Direction','reverse') 
% for i=1:length(atex)
%     text(3.4,(-1.4+(i-1)*0.6),sprintf(atex{i}))
% end
title('All fishes')
print('-dpng',[ppath 'Map_LMEs_driver_mlr_nu_cluster_v2_A.png'])

%%
f1 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Fcorr)
colormap(mcol2)
caxis([1 13])
title('Forage fishes')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pcorr)
colormap(mcol2)
caxis([1 13])
title('Large pelagics')

subplot('Position',[0.33 0.10 0.32 0.4]) %A
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Acorr)
colormap(mcol2)
caxis([1 13])
title('All fishes')

subplot('Position',[0.01 0.10 0.32 0.4]) %D
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Dcorr)
colormap(mcol2)
caxis([1 13])
title('Demersals')

colorbar('Position',[0.675 0.125 0.03 0.35],'Ticks',1:13,'TickLabels',alltex,...
    'Direction','reverse')

print('-dpng',[ppath 'Map_LMEs_driver_mlr_nu_cluster_v2_fntypes.png'])

