% Map clusters of driver-fish mult linear regressions
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
%spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/' cfile '/corrs/'];

spath ='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';
load([spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_alllags_reduc.mat'])

%% Cluster descrips
% ftex = {'+ZmLoss',...             = 2  
%     '+TP',...                     = 1
%     '+Zmeso',...                  = 3
%     '+TP*Zmeso, -TP*ZmLoss',...   = 5
%     '+Zmeso, -ZmLoss'};           = 4
ftex = {'+TP',...             = 1
    '+ZmLoss',...             = 2
    '+Zmeso',...              = 3
    '+Zmeso, -ZmLoss',...     = 4
    '+TP*Zmeso, -TP*ZmLoss'};% = 5
    
% ptex = {'-TP, +Zmeso',...         = 4
%     '+TP',...                     = 1
%     '-TP',...                     = 5
%     '-Zmeso, +ZmLoss, (+TP)',...  = 3
%     '+Zmeso, -ZmLoss'};           = 2
ptex = {'+TP',...               = 1
    '+Zmeso, -ZmLoss',...       = 2
    '-Zmeso, +ZmLoss, (+TP)',...= 3
    '-TP, +Zmeso',...           = 4
    '-TP'};                %    = 5
    
% dtex = {'+Zmeso, -ZmLoss',...                 = 3
%     '+Det',...                                = 1
%     '-Zmeso, +ZmLoss',...                     = 2
%     '-Det?',...                               = 4
%     '-Zmeso, +ZmLoss, -TP*Zmeso, +TP*ZmLoss'};= 5
dtex = {'+Det',...                           = 1
    '-Zmeso, +ZmLoss',...                    = 2
    '+Zmeso, -ZmLoss',...                    = 3
    '-Det?',...                              = 4
    '-Zm, +ZmL, -TP*Zm, +TP*ZmL'}; %        = 5

% atex = {'+TP',...                 = 1
%     '-Zmeso, +ZmLoss, +TP',...    = 2
%     '-Zmeso, +ZmLoss, -TP',...    = 6
%     '+Zmeso',...                  = 3
%     '+Zmeso, -ZmLoss',...         = 4
%     '++Zmeso, --ZmLoss'}; %       = 5
atex = {'+TP',...                 = 1
    '-Zmeso, +ZmLoss, +TP',...    = 2
    '+Zmeso',...                  = 3
    '+Zmeso, -ZmLoss',...         = 4
    '++Zmeso, --ZmLoss',...       = 5
    '-Zmeso, +ZmLoss, -TP'}; %    = 6


%%  ---------------------------------------------------------
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

mcol = [
    34/255 136/255 51/255;...   %green
    170/255 51/255 119/255;...  %purple
    238/255 102/255 119/255;... %red
    0/255 68/255 136/255;...    %blue
    51/255 187/255 238/255;...  %cyan
    153/255 153/255 51/255;...  %olive
    0 0 0;...                   %black
    0.50 0.50 0.50;...          % grey
    ];


%colorblind friendly
% cb=[34/255 136/255 51/255;...   %green
%     153/255 153/255 51/255;...  %olive
%     51/255 187/255 238/255;...  %cyan
%     0/255 68/255 136/255;...    %blue
%     238/255 102/255 119/255;... %red
%     170/255 51/255 119/255;...  %purple
%     0 0 0;...                   %black
%     0.25 0.25 0.25;...          %dk grey
%     0.50 0.50 0.50;...          % grey
%     0.75 0.75 0.75];            %lt grey

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
    Fcorr(id) = ClusterF(i);
    Pcorr(id) = ClusterP(i);
    Dcorr(id) = ClusterD(i);
    Acorr(id) = ClusterA(i); 
end

%% All
figure(1)
subplot('Position',[0.01 0.13 0.725 0.725])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.85 0.85 0.85]);
hold on
surfm(TLAT,TLONG,Acorr)
colormap(turbo(6))
colorbar('Position',[0.74 0.25 0.03 0.5],'Ticks',1:6,'TickLabels',[]) %,'TickLabels',atex)
for i=1:length(atex)
    text(3.4,(-1.4+(i-1)*0.6),sprintf(atex{i}))
end
% colormap(mcol)
% colorbar('Ticks',1:8,'TickLabels',ctex)
title('All fishes')
print('-dpng',[ppath 'Map_LMEs_driver_mlr_nu_cluster_A.png'])

%% F
figure(2)
subplot('Position',[0.01 0.13 0.725 0.725])%F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.85 0.85 0.85]);
hold on
surfm(TLAT,TLONG,Fcorr)
colormap(copper(5))
colorbar('Position',[0.74 0.25 0.03 0.5],'Ticks',1:5,'TickLabels',[]) %,'TickLabels',atex)
for i=1:length(ftex)
    text(3.4,(-1.4+(i-1)*0.725),sprintf(ftex{i}))
end
title('Forage fishes')
print('-dpng',[ppath 'Map_LMEs_driver_mlr_nu_cluster_F.png'])

%% P
figure(3)
subplot('Position',[0.01 0.13 0.725 0.725]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.85 0.85 0.85]);
hold on
surfm(TLAT,TLONG,Pcorr)
colormap(parula(5))
colorbar('Position',[0.74 0.25 0.03 0.5],'Ticks',1:4,'TickLabels',[]) %,'TickLabels',atex)
for i=1:length(ptex)
    text(3.4,(-1.4+(i-1)*0.725),sprintf(ptex{i}))
end
title('Large pelagics')
print('-dpng',[ppath 'Map_LMEs_driver_mlr_nu_cluster_P.png'])

%% D
figure(4)
subplot('Position',[0.01 0.13 0.725 0.725]) %D
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.85 0.85 0.85]);
hold on
surfm(TLAT,TLONG,Dcorr)
colormap(parula(5))
colorbar('Position',[0.74 0.25 0.03 0.5],'Ticks',1:5,'TickLabels',[]) %,'TickLabels',atex)
for i=1:length(dtex)
    text(3.4,(-1.4+(i-1)*0.725),sprintf(dtex{i}))
end
title('Demersals')
print('-dpng',[ppath 'Map_LMEs_driver_mlr_nu_cluster_D.png'])


