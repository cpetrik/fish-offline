% Map clusters of driver-fish mult linear regressions
% For all 63 LMEs
% CESM Tp, Tb, Det, ZmLoss
% Fish biomass, nu, CPUE

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/' cfile '/corrs/'];

mod = 'v15_All_fish03';

%spath ='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';
load([spath,'LME_biom_nu_cpue_A_sig_corr_coeffs_cluster.mat'])

%% vector of clusters to use old code
ClusterB = biomA(:,6);
ClusterP = prodA(:,6);
ClusterE = cpueA(:,6);

%% Cluster descrips
% alltex = 

btex = {'-T, +prey',...   %1
    '+prey'};             %2
    
ptex = {'+T',...        %1
        '-T, +prey',... %2
        '+prey',...     %3
        '~'};           %4
    
    
    
%cpue
etex = {'-ZL',...   %1
    '-T',...        %2
    '-T, +prey',... %3
    '~',...         %4
    '+T',...        %5
    '+prey'};       %6
    
             
   
%% Add a text descript vector that matches to clusters
% 
%biomass
CatB = cell(length(ClusterB),1);
CatB(ClusterB==1) = btex(1);
CatB(ClusterB==2) = btex(2);
CatB(isnan(ClusterB)) = {''};

%nu, prod
CatP = cell(length(ClusterP),1);
CatP(ClusterP==1) = ptex(1);
CatP(ClusterP==2) = ptex(2);
CatP(ClusterP==3) = ptex(3);
CatP(ClusterP==4) = ptex(4);
CatP(isnan(ClusterP)) = {''};

%cpue
CatE = cell(length(ClusterE),1);
CatE(ClusterE==1) = etex(1);
CatE(ClusterE==2) = etex(2);
CatE(ClusterE==3) = etex(3);
CatE(ClusterE==4) = etex(4);
CatE(ClusterE==5) = etex(5);
CatE(ClusterE==6) = etex(6);
CatE(isnan(ClusterE)) = {''};

CatB = string(CatB);
CatP = string(CatP);
CatE = string(CatE);

% Put text in table
cCatB = char(CatB);
cCatP = char(CatP);
cCatE = char(CatE);

Atab = table(cCatB,cCatP,cCatE,'VariableNames',...
    {'Biomass','Production','CPUE'});

writetable(Atab,[spath,'LME_biom_nu_cpue_A_corr_coeffs_cluster_matched.csv'])


%% Create one master colormap for all categories
alltex = {'+T',...  %1
'-T',...            %2
'-T, +prey',...     %3
'+prey',...         %4
'-ZL'};             %5

ClusterB(:,2) = nan;
ClusterB((CatB==alltex(1)),2) = 1;
ClusterB((CatB==alltex(2)),2) = 2;
ClusterB((CatB==alltex(3)),2) = 3;
ClusterB((CatB==alltex(4)),2) = 4;
ClusterB((CatB==alltex(5)),2) = 5;


ClusterP(:,2) = nan;
ClusterP((CatP==alltex(1)),2) = 1;
ClusterP((CatP==alltex(2)),2) = 2;
ClusterP((CatP==alltex(3)),2) = 3;
ClusterP((CatP==alltex(4)),2) = 4;
ClusterP((CatP==alltex(5)),2) = 5;

ClusterE(:,2) = nan;
ClusterE((CatE==alltex(1)),2) = 1;
ClusterE((CatE==alltex(2)),2) = 2;
ClusterE((CatE==alltex(3)),2) = 3;
ClusterE((CatE==alltex(4)),2) = 4;
ClusterE((CatE==alltex(5)),2) = 5;

%% save
save([spath,'LME_biom_nu_gamma_catch_cpue_A_corr_coeffs_cluster_matched.mat'])

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

%%
% on grid
Bclus  = nan(ni,nj);
Pclus  = nan(ni,nj);
Eclus  = nan(ni,nj);

lid = biomA(:,1);
%use col = 2 after finding matching # to text
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Eclus(id) = ClusterE(i,1);
    Pclus(id) = ClusterP(i,1);
    Bclus(id) = ClusterB(i,1);
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
colormap(mcol(1:6,:))
clim([1 6])
colorbar('southoutside','Ticks',1:5,'TickLabels',1:5,'Direction','reverse')
title('Biom')

%prod
subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pclus)
colormap(mcol(1:6,:))
clim([1 6])
colorbar('southoutside','Ticks',1:5,'TickLabels',1:5,'Direction','reverse')
title('Prod')

%Gam
%subplot('Position',[0.65 0.575 0.32 0.4])


% catch
%subplot('Position',[0.01 0.10 0.32 0.4]) 


%cpue
subplot('Position',[0.33 0.10 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Eclus)
colormap(mcol(1:6,:))
clim([1 6])
colorbar('southoutside','Ticks',1:6,'TickLabels',1:6,'Direction','reverse')
title('CPUE')

print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_cpue_A_corr_coeffs_cluster_Match.png'])

%% on grid after matching
Bcorr  = nan(ni,nj);
Pcorr  = nan(ni,nj);
Ecorr  = nan(ni,nj);

lid = biomA(:,1);
%use col = 2 after finding matching # to text
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Bcorr(id) = ClusterB(i,2);
    Pcorr(id) = ClusterP(i,2);
    Ecorr(id) = ClusterE(i,2);
end

%%
mcol2 = [
    238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    0/255 68/255 136/255;...    %blue
    51/255 187/255 238/255;...  %cyan
    34/255 136/255 51/255];  %green

%% All same 5 colors and text ---------------------------

f1 = figure('Units','inches','Position',[1 3 5.5 5]);

subplot('Position',[0.01 0.51 0.45 0.45]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Bcorr)
%colormap(mcol(1:5,:))
colormap(mcol2)
clim([1 5])
title('Biomass')

subplot('Position',[0.5 0.51 0.45 0.45]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pcorr)
%colormap(mcol(1:5,:))
colormap(mcol2)
clim([1 5])
title('Prod')

%subplot('Position',[0.65 0.575 0.32 0.4]) 


%subplot('Position',[0.01 0.10 0.32 0.4])


subplot('Position',[0.25 0.05 0.45 0.45]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Ecorr)
%colormap(mcol(1:5,:))
colormap(mcol2)
clim([1 5])
title('CPUE')
colorbar('Position',[0.725 0.1 0.03 0.35],'Ticks',1:12,'TickLabels',alltex,...
    'Direction','reverse')
print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_cpue_A_corr_coeffs_cluster.png'])

%% Agree +T --> change to Biom and/or Prod agree with CPUE
%'+T',...  %1
id1 = find(Bcorr==1);
id2 = find(Pcorr==1);
id3 = find(Ecorr==1);
tid = intersect(id2,id3);

Bcorr2 = nan*ones(ni,nj);
Pcorr2 = nan*ones(ni,nj);
Pcorr2(tid) = 1;
Ecorr2 = nan*ones(ni,nj);
Ecorr2(tid) = 1;

f5 = figure('Units','inches','Position',[1 3 5.5 5]);

subplot('Position',[0.01 0.51 0.45 0.45]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Bcorr2)
%colormap(mcol(1:5,:))
colormap(mcol2)
clim([1 5])
title('Biomass')

subplot('Position',[0.5 0.51 0.45 0.45]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pcorr2)
%colormap(mcol(1:5,:))
colormap(mcol2)
clim([1 5])
title('Prod')

subplot('Position',[0.25 0.05 0.45 0.45]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Ecorr2)
%colormap(mcol(1:5,:))
colormap(mcol2)
clim([1 5])
title('CPUE')
colorbar('Position',[0.725 0.1 0.03 0.35],'Ticks',1:12,'TickLabels',alltex,...
    'Direction','reverse')
print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_cpue_A_corr_coeffs_posT.png'])

%% Agree -T, +prey --> change to Biom and/or Prod agree with CPUE
%'-T, +prey',...     %3
id1 = find(Bcorr==3);
id2 = find(Pcorr==3);
id3 = find(Ecorr==3);
bid = intersect(id1,id3);
mid = intersect(id2,id3);

Bcorr3 = nan*ones(ni,nj);
Bcorr3(bid) = 3;
Pcorr3 = nan*ones(ni,nj);
Pcorr3(mid) = 3;
Ecorr3 = nan*ones(ni,nj);
Ecorr3(bid) = 3;
Ecorr3(mid) = 3;

f2 = figure('Units','inches','Position',[1 3 5.5 5]);

subplot('Position',[0.01 0.51 0.45 0.45]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Bcorr3)
%colormap(mcol(1:5,:))
colormap(mcol2)
clim([1 5])
title('Biomass')

subplot('Position',[0.5 0.51 0.45 0.45]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pcorr3)
%colormap(mcol(1:5,:))
colormap(mcol2)
clim([1 5])
title('Prod')

subplot('Position',[0.25 0.05 0.45 0.45]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Ecorr3)
%colormap(mcol(1:5,:))
colormap(mcol2)
clim([1 5])
title('CPUE')
colorbar('Position',[0.725 0.1 0.03 0.35],'Ticks',1:12,'TickLabels',alltex,...
    'Direction','reverse')
print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_cpue_A_corr_coeffs_negTposP.png'])

%% Agree +prey --> change to Biom and/or Prod agree with CPUE
%'+prey',...         %4
id1 = find(Bcorr==4);
id2 = find(Pcorr==4);
id3 = find(Ecorr==4);
pid = intersect(id1,id3);
zid = intersect(id2,id3);


Bcorr4 = nan*ones(ni,nj);
Bcorr4(pid) = 4;
Pcorr4 = nan*ones(ni,nj);
Pcorr4(zid) = 4;
Ecorr4 = nan*ones(ni,nj);
Ecorr4(pid) = 4;
Ecorr4(zid) = 4;

f3 = figure('Units','inches','Position',[1 3 5.5 5]);

subplot('Position',[0.01 0.51 0.45 0.45]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Bcorr4)
%colormap(mcol(1:5,:))
colormap(mcol2)
clim([1 5])
title('Biomass')

subplot('Position',[0.5 0.51 0.45 0.45]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pcorr4)
%colormap(mcol(1:5,:))
colormap(mcol2)
clim([1 5])
title('Prod')

subplot('Position',[0.25 0.05 0.45 0.45]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Ecorr4)
%colormap(mcol(1:5,:))
colormap(mcol2)
clim([1 5])
title('CPUE')
colorbar('Position',[0.725 0.1 0.03 0.35],'Ticks',1:12,'TickLabels',alltex,...
    'Direction','reverse')
print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_cpue_A_corr_coeffs_posP.png'])

