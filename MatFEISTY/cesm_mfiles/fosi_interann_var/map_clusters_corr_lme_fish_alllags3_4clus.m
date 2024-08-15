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
load([spath,'LME_biom_nu_gamma_catch_cpue_A_corr_coeffs_cluster.mat'])

%% vector of clusters to use old code
ClusterB = biomA(:,6);
ClusterP = prodA(:,6);
ClusterG = gammaA(:,6);
ClusterC = catchA(:,6);
ClusterE = cpueA(:,6);

%% Cluster descrips
% alltex = 

btex = {'-T',...          %1
    '-T, +prey',...       %2
    '+T, +prey',...       %3
    '+prey'};             %4
    
ptex = {'+T, -prey',...  %1
    '-T, +prey',...     %2
    '+prey',...         %3
    '+T',...            %4
    '+ZL'};             %5
    
gtex = {'+T, -Det, +ZL',... %1
    '+T, -prey',...         %2
    '-T, +prey',...         %3
    '+T, +prey'};           %4

%catch
ctex = {'+T, -prey',... %1
    '-T, +prey',...    %2
    'weak',...          %3
    '+T, +prey'};       %4

%cpue
etex = {'-T, +prey',...  %1
    '-T, -Det',...      %2    
    '-ZL',...           %3
    '+T, +prey'};       %4
             
   
%% Add a text descript vector that matches to clusters
% 
%biomass
CatB = cell(length(ClusterB),1);
CatB(ClusterB==1) = btex(1);
CatB(ClusterB==2) = btex(2);
CatB(ClusterB==3) = btex(3);
CatB(ClusterB==4) = btex(4);
CatB(isnan(ClusterB)) = {''};

%nu, prod
CatP = cell(length(ClusterP),1);
CatP(ClusterP==1) = ptex(1);
CatP(ClusterP==2) = ptex(2);
CatP(ClusterP==3) = ptex(3);
CatP(ClusterP==4) = ptex(4);
CatP(ClusterP==5) = ptex(5);
CatP(isnan(ClusterP)) = {''};

%gamma, recruit flux
CatG = cell(length(ClusterG),1);
CatG(ClusterG==1) = gtex(1);
CatG(ClusterG==2) = gtex(2);
CatG(ClusterG==3) = gtex(3);
CatG(ClusterG==4) = gtex(4);
CatG(isnan(ClusterG)) = {''};

%catch
CatC = cell(length(ClusterC),1);
CatC(ClusterC==1) = ctex(1);
CatC(ClusterC==2) = ctex(2);
CatC(ClusterC==3) = ctex(3);
CatC(ClusterC==4) = ctex(4);
CatC(isnan(ClusterC)) = {''};

%cpue
CatE = cell(length(ClusterE),1);
CatE(ClusterE==1) = etex(1);
CatE(ClusterE==2) = etex(2);
CatE(ClusterE==3) = etex(3);
CatE(ClusterE==4) = etex(4);
CatE(isnan(ClusterE)) = {''};

CatB = string(CatB);
CatP = string(CatP);
CatG = string(CatG);
CatC = string(CatC);
CatE = string(CatE);

% Put text in table
cCatB = char(CatB);
cCatP = char(CatP);
cCatG = char(CatG);
cCatC = char(CatC);
cCatE = char(CatE);

Atab = table(cCatB,cCatP,cCatG,cCatC,cCatE,'VariableNames',...
    {'Biomass','Production','Recruitment','Catch','CPUE'});

writetable(Atab,[spath,'LME_biom_nu_gamma_catch_cpue_A_corr_coeffs_cluster_matched.csv'])


%% Create one master colormap for all categories
alltex = {'-T',...  %1
'-T, +prey',...     %2
'-T, -Det',...      %3
'+T',...            %4
'+T, +prey',...     %5
'+T, -prey',...     %6
'+T, -Det, +ZL',... %7
'+prey',...         %8
'+ZL',...           %9
'-ZL'};             %10

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

ClusterG(:,2) = nan;
ClusterG((CatG==alltex(1)),2) = 1;
ClusterG((CatG==alltex(2)),2) = 2;
ClusterG((CatG==alltex(3)),2) = 3;
ClusterG((CatG==alltex(4)),2) = 4;
ClusterG((CatG==alltex(5)),2) = 5;
ClusterG((CatG==alltex(6)),2) = 6;
ClusterG((CatG==alltex(7)),2) = 7;
ClusterG((CatG==alltex(8)),2) = 8;
ClusterG((CatG==alltex(9)),2) = 9;
ClusterG((CatG==alltex(10)),2) = 10;

ClusterC(:,2) = nan;
ClusterC((CatC==alltex(1)),2) = 1;
ClusterC((CatC==alltex(2)),2) = 2;
ClusterC((CatC==alltex(3)),2) = 3;
ClusterC((CatC==alltex(4)),2) = 4;
ClusterC((CatC==alltex(5)),2) = 5;
ClusterC((CatC==alltex(6)),2) = 6;
ClusterC((CatC==alltex(7)),2) = 7;
ClusterC((CatC==alltex(8)),2) = 8;
ClusterC((CatC==alltex(9)),2) = 9;
ClusterC((CatC==alltex(10)),2) = 10;

ClusterE(:,2) = nan;
ClusterE((CatE==alltex(1)),2) = 1;
ClusterE((CatE==alltex(2)),2) = 2;
ClusterE((CatE==alltex(3)),2) = 3;
ClusterE((CatE==alltex(4)),2) = 4;
ClusterE((CatE==alltex(5)),2) = 5;
ClusterE((CatE==alltex(6)),2) = 6;
ClusterE((CatE==alltex(7)),2) = 7;
ClusterE((CatE==alltex(8)),2) = 8;
ClusterE((CatE==alltex(9)),2) = 9;
ClusterE((CatE==alltex(10)),2) = 10;

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
Gclus  = nan(ni,nj);
Cclus  = nan(ni,nj);
Eclus  = nan(ni,nj);

lid = biomA(:,1);
%use col = 2 after finding matching # to text
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Eclus(id) = ClusterE(i,1);
    Pclus(id) = ClusterP(i,1);
    Gclus(id) = ClusterG(i,1);
    Cclus(id) = ClusterC(i,1);
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

print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_gamma_catch_cpue_A_corr_coeffs_cluster_Match.png'])

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
print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_gamma_catch_cpue_A_corr_coeffs_cluster.png'])

