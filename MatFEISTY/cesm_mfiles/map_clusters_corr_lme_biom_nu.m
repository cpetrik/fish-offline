% Map clusters of driver-fish mult linear regressions
% For all 63 LMEs
% CESM Tp, Tb, Det, ZmLoss
% Fish biomass, nu

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

CatB = string(CatB);
CatP = string(CatP);

% Put text in table
cCatB = char(CatB);
cCatP = char(CatP);

Atab = table(cCatB,cCatP,'VariableNames',...
    {'Biomass','Production'});

writetable(Atab,[spath,'LME_biom_nu_A_corr_coeffs_cluster_matched.csv'])


%% Create one master colormap for all categories
alltex = {'+T',...  %1
'-T, +prey',...     %2
'+prey'};           %3

ClusterB(:,2) = nan;
ClusterB((CatB==alltex(1)),2) = 1;
ClusterB((CatB==alltex(2)),2) = 2;
ClusterB((CatB==alltex(3)),2) = 3;


ClusterP(:,2) = nan;
ClusterP((CatP==alltex(1)),2) = 1;
ClusterP((CatP==alltex(2)),2) = 2;
ClusterP((CatP==alltex(3)),2) = 3;


%% save
save([spath,'LME_biom_nu_gamma_A_corr_coeffs_cluster_matched.mat'])

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
%mcol = muted ./ 255;

%try night
mcol(1,:) = night(14,:)./255;
mcol(2,:) = night(1,:)./255;
mcol(3,:) = night(5,:)./255;

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

lid = biomA(:,1);
%use col = 2 after finding matching # to text
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
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
colormap(mcol(1:4,:))
clim([1 4])
colorbar('southoutside','Ticks',1:5,'TickLabels',1:5,'Direction','reverse')
title('Biom')

%prod
subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pclus)
colormap(mcol(1:4,:))
clim([1 4])
colorbar('southoutside','Ticks',1:5,'TickLabels',1:5,'Direction','reverse')
title('Prod')

print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_A_corr_coeffs_cluster_Match.png'])

%% on grid after matching
Bcorr  = nan(ni,nj);
Pcorr  = nan(ni,nj);

lid = biomA(:,1);
%use col = 2 after finding matching # to text
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Bcorr(id) = ClusterB(i,2);
    Pcorr(id) = ClusterP(i,2);
end


%% All same 5 colors and text ---------------------------

f1 = figure('Units','inches','Position',[1 3 7.5 2.75]);
subplot('Position',[0.01 0.05 0.425 0.8])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Bcorr)
colormap(mcol)
clim([1 3])
title('Biomass')

subplot('Position',[0.45 0.05 0.425 0.8])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,Pcorr)
colormap(mcol)
clim([1 3])
title('Production')
colorbar('Position',[0.88 0.1 0.03 0.7],'Ticks',1:3,'TickLabels',alltex,...
    'Direction','reverse')
print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_A_corr_coeffs_cluster.png'])

