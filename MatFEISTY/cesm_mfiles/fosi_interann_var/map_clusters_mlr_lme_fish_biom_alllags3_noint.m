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
load([spath,'LME_nu_mlr_coeffs_ALLdiv2SD_alllags3_R2_07_cluster.mat'])

%% vector of clusters to use old code
ClusterF = Fbiom(:,22);
ClusterP = Pbiom(:,22);
ClusterD = Dbiom(:,22);
ClusterB = Bbiom(:,22);
ClusterA = Abiom(:,22);

%% Cluster descrips
% alltex = {'+TB, -TP',...    = 1
%     '-TP',...               = 2
%     '-TP, +ZL',...          = 3
%     '+ZL',...               = 4
%     '-Det, +ZL',...         = 5
%     '+Det, +ZL',...         = 6
%     '+Det, -ZL',...         = 7
%     '+Det',...              = 8
%     '+Det &',...            = 9
%     '-Det',...              = 10
%     '+TB, -TP, +ZL',...     = 11
%     '-TB, +TP'};          % = 12

ftex = {'+ZL',...       = 1
    '-TP',...           = 2 
    'other',...         = 3
    '-Det',...          = 4
    '+ZL, +ZL',...      = 5
    '+TB',...           = 6 
    '+TB, -TP'};    %   = 7
    
ptex = {'-TP, +ZL',...    = 1
    '-Det, -Det, +ZL',... = 2
    '+ZL',...             = 3 
    '+TB, -TP, +ZL',...   = 4
    '-Det, +ZL',...       = 5
    '+Det, -ZL',...       = 6
    'other'};           % = 7
    
dtex = {'+Det, +Det',...      = 1 
    '+Det, +ZL',...           = 2
    '+Det & others',...       = 3
    '+Det, +Det, +ZL'};    %  = 4

atex = {'+Det, +Det',... = 1
    '-TP',...            = 2
    '+ZL, +ZL',...       = 3
    '-Det, +ZL',...      = 4
    '+Det & others',...  = 5
    '+ZL'}; %            = 6
     

btex = {'+Det0, +Det1, -Det2?',...  = 1
    '+Det, +Det',...                = 2
    '+Det0, -Det2',...              = 3
    '+TB, -TP',...                  = 4
    '-TB, +TP'};      %             = 5
          

%% Add a text descript vector that matches to clusters
% F = 7
CatF = cell(length(ClusterF),1);
CatF(ClusterF==1) = ftex(1);
CatF(ClusterF==2) = ftex(2);
CatF(ClusterF==3) = ftex(3);
CatF(ClusterF==4) = ftex(4);
CatF(ClusterF==5) = ftex(5);
CatF(ClusterF==6) = ftex(6);
CatF(ClusterF==7) = ftex(7);
CatF(isnan(ClusterF)) = {''};

% P = 7
CatP = cell(length(ClusterP),1);
CatP(ClusterP==1) = ptex(1);
CatP(ClusterP==2) = ptex(2);
CatP(ClusterP==3) = ptex(3);
CatP(ClusterP==4) = ptex(4);
CatP(ClusterP==5) = ptex(5);
CatP(ClusterP==6) = ptex(6);
CatP(ClusterP==7) = ptex(7);
CatP(isnan(ClusterP)) = {''};

% D = 4
CatD = cell(length(ClusterD),1);
CatD(ClusterD==1) = dtex(1);
CatD(ClusterD==2) = dtex(2);
CatD(ClusterD==3) = dtex(3);
CatD(ClusterD==4) = dtex(4);
CatD(isnan(ClusterD)) = {''};

% A = 6
CatA = cell(length(ClusterA),1);
CatA(ClusterA==1) = atex(1);
CatA(ClusterA==2) = atex(2);
CatA(ClusterA==3) = atex(3);
CatA(ClusterA==4) = atex(4);
CatA(ClusterA==5) = atex(5);
CatA(ClusterA==6) = atex(6);
CatA(isnan(ClusterA)) = {''};

% B = 5
CatB = cell(length(ClusterB),1);
CatB(ClusterB==1) = btex(1);
CatB(ClusterB==2) = btex(2);
CatB(ClusterB==3) = btex(3);
CatB(ClusterB==4) = btex(4);
CatB(ClusterB==5) = btex(5);
CatB(isnan(ClusterB)) = {''};

CatF = string(CatF);
CatP = string(CatP);
CatD = string(CatD);
CatA = string(CatA);
CatB = string(CatB);


%% Create one master colormap for all categories
alltex = {'+TB, -TP',...    = 1
    '-TP',...               = 2
    '-TP, +ZL',...          = 3
    '+ZL',...               = 4
    '-Det, +ZL',...         = 5
    '+Det, +ZL',...         = 6
    '+Det, -ZL',...         = 7
    '+Det',...              = 8
    '+Det &',...            = 9
    '-Det',...              = 10
    '+TB, -TP, +ZL',...     = 11
    '-TB, +TP'};          % = 12


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

ClusterB(:,2) = nan;
ClusterB((CatB==alltex(1)),2) = 1;
ClusterB((CatB==alltex(2)),2) = 2;
ClusterB((CatB==alltex(3)),2) = 3;
ClusterB((CatB==alltex(4)),2) = 4;
ClusterB((CatB==alltex(5)),2) = 5;
ClusterB((CatB==alltex(6)),2) = 6;
ClusterB((CatB==alltex(7)),2) = 7;
ClusterB((CatB==alltex(8)),2) = 8;
ClusterB((CatB==alltex(9)),2) = 9;
ClusterB((CatB==alltex(10)),2) = 10;
ClusterB((CatB==alltex(11)),2) = 11;
ClusterB((CatB==alltex(12)),2) = 12;


%%  Colormap

% mcol = [
%     34/255 136/255 51/255;...   %green
%     170/255 51/255 119/255;...  %purple
%     238/255 102/255 119/255;... %red
%     0/255 68/255 136/255;...    %blue
%     51/255 187/255 238/255;...  %cyan
%     153/255 153/255 51/255;...  %olive
%     0 0 0;...                   %black
%     0.50 0.50 0.50;...          % grey
%     ];


% colorblind friendly
load('paul_tol_cmaps.mat')

%try night w/16 colors
% mcol = night(1:12,:) ./ 255;
% %reorder to put light colors as unknowns
% mcol2 = mcol(1:6,:);
% mcol2(7:9,:) = mcol(10:12,:);
% mcol2(10:12,:) = mcol(7:9,:);

%try muted and add greys
mcol = muted ./ 255;
%add greys as unknowns
mcol(10,:) = zmeso(10,:);
mcol(11,:) = zmeso(9,:);
mcol(12,:) = zmeso(8,:);

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
Bcorr  = nan(ni,nj);

lid = Fbiom(:,1);
%usse col = 2 after finding matching # to text
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Fcorr(id) = ClusterF(i,1);
    Pcorr(id) = ClusterP(i,1);
    Dcorr(id) = ClusterD(i,1);
    Acorr(id) = ClusterA(i,1);
    Bcorr(id) = ClusterB(i,1);
end

%% Matching
% F
f1 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Fcorr)
colormap(mcol(1:7,:))
clim([1 7])
colorbar('southoutside','Ticks',1:7,'TickLabels',ftex,'Direction','reverse')
title('Forage fishes')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pcorr)
colormap(mcol(1:7,:))
clim([1 7])
colorbar('southoutside','Ticks',1:7,'TickLabels',ptex,'Direction','reverse')
title('Large pelagics')

subplot('Position',[0.65 0.575 0.32 0.4]) %A
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Acorr)
colormap(mcol(1:7,:))
clim([1 7])
colorbar('southoutside','Ticks',1:6,'TickLabels',atex,'Direction','reverse')
title('All fishes')

subplot('Position',[0.01 0.10 0.32 0.4]) %D
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Dcorr)
colormap(mcol(1:7,:))
clim([1 7])
colorbar('southoutside','Ticks',1:4,'TickLabels',dtex,'Direction','reverse')
title('Demersals')

subplot('Position',[0.33 0.10 0.32 0.4]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Bcorr)
colormap(mcol(1:7,:))
clim([1 7])
colorbar('southoutside','Ticks',1:5,'TickLabels',btex,'Direction','reverse')
title('Benthos')
print('-dpng',[ppath 'Map_LMEs_driver_mlr_biom_alllags3_noint_cluster_fntypes.png'])


%% All same 12 colors and text ---------------------------
% All
figure(10)
subplot('Position',[0.01 0.13 0.725 0.725])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Acorr)
colormap(mcol)
caxis([1 12])
colorbar('Position',[0.71 0.25 0.03 0.5],'Ticks',1:12,'TickLabels',alltex2,...
    'Direction','reverse') 
% for i=1:length(atex)
%     text(3.4,(-1.4+(i-1)*0.6),sprintf(atex{i}))
% end
title('All fishes')
print('-dpng',[ppath 'Map_LMEs_driver_mlr_biom_cluster_v2_A.png'])

%%
f1 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Fcorr)
colormap(mcol)
caxis([1 12])
title('Forage fishes')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pcorr)
colormap(mcol)
caxis([1 12])
title('Large pelagics')

subplot('Position',[0.65 0.575 0.32 0.4]) %A
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Acorr)
colormap(mcol)
caxis([1 12])
title('All fishes')

subplot('Position',[0.01 0.10 0.32 0.4]) %D
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Dcorr)
colormap(mcol)
caxis([1 12])
title('Demersals')

subplot('Position',[0.33 0.10 0.32 0.4]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Bcorr)
colormap(mcol)
caxis([1 12])
colorbar('Position',[0.675 0.125 0.03 0.35],'Ticks',1:12,'TickLabels',alltex2,...
    'Direction','reverse')
title('Benthos')
print('-dpng',[ppath 'Map_LMEs_driver_mlr_biom_cluster_v2_fntypes.png'])