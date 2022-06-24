% Map correlations of FEISTY FOSI
% w/climate indices

clear all
close all

apath = '/Users/cpetrik/Dropbox/Princeton/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'Climate_anomalies_annual_means.mat']);

tanom = canom;
clear canom

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
clatlim=[plotminlat plotmaxlat];
clonlim=[plotminlon plotmaxlon];

load coastlines;

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
dpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

sims = {'v15_All_fish03_';'v15_climatol_';'v15_varFood_';'v15_varTemp_'};
mod = sims{1};

%'Annual_Means_FOSI_v15_All_fish03_Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100'
load([dpath 'Annual_Means_FOSI_',mod,cfile '.mat']);

%% Group types and sizes
allS = sf_abio+sp_abio+sd_abio;
allM = mf_abio+mp_abio+md_abio;
allL = lp_abio+ld_abio;
allF = sf_abio+mf_abio;
allP = sp_abio+mp_abio+lp_abio;
allD = sd_abio+md_abio+ld_abio;
allV = allF + allP + allD;
%allB = b_abio;

clear sf_abio sp_abio sd_abio mf_abio mp_abio md_abio lp_abio ld_abio

%% Calc linear trend and remove

%% Calc fish anomalies
%mean biomass
mba = b_abio - nanmean(b_abio,2);
mFa = allF - nanmean(allF,2);
mPa = allP - nanmean(allP,2);
mDa = allD - nanmean(allD,2);
mSa = allS - nanmean(allS,2);
mMa = allM - nanmean(allM,2);
mLa = allL - nanmean(allL,2);
mVa = allV - nanmean(allV,2);

%%
clear allS allM allL allF allP allD allV

%% rep climate anoms
canom = repmat(manom(1,:),length(ID),1);
canom(:,:,2) = repmat(manom(2,:),length(ID),1);
canom(:,:,3) = repmat(manom(3,:),length(ID),1);
canom(:,:,4) = repmat(manom(4,:),length(ID),1);
canom(:,:,5) = repmat(manom(5,:),length(ID),1);
canom(:,:,6) = repmat(manom(6,:),length(ID),1);
canom(:,:,7) = repmat(manom(7,:),length(ID),1);
canom(:,:,8) = repmat(manom(8,:),length(ID),1);
canom(:,:,9) = repmat(manom(9,:),length(ID),1);
canom(:,:,10) = repmat(manom(10,:),length(ID),1);
canom(:,:,11) = repmat(manom(11,:),length(ID),1);

%% correlations by grid cell
yr = 1:3; %lags
rS = nan*ones(length(ID),11,length(yr));
rM = rS;
rL = rS;
rF = rS;
rP = rS;
rD = rS;
rV = rS;
rB = rS;

for j=1:11
    
    for k=1:3
        t = yr(k);
        %% Corr at diff lags
        rS1 = diag(corr((mSa(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS2 = diag(corr((mSa(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS3 = diag(corr((mSa(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS(:,j,k) = [rS1;rS2;rS3];
        clear rS1 rS2 rS3
        
        rM1 = diag(corr((mMa(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rM2 = diag(corr((mMa(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rM3 = diag(corr((mMa(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rM(:,j,k) = [rM1;rM2;rM3];
        clear rM1 rM2 rM3
        
        rL1 = diag(corr((mLa(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rL2 = diag(corr((mLa(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rL3 = diag(corr((mLa(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rL(:,j,k) = [rL1;rL2;rL3];
        clear rL1 rL2 rL3
        
        rF1 = diag(corr((mFa(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rF2 = diag(corr((mFa(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rF3 = diag(corr((mFa(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rF(:,j,k) = [rF1;rF2;rF3];
        clear rF1 rF2 rF3
        
        rP1 = diag(corr((mPa(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rP2 = diag(corr((mPa(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rP3 = diag(corr((mPa(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rP(:,j,k) = [rP1;rP2;rP3];
        clear rP1 rP2 rP3
        
        rD1 = diag(corr((mDa(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rD2 = diag(corr((mDa(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rD3 = diag(corr((mDa(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rD(:,j,k) = [rD1;rD2;rD3];
        clear rD1 rD2 rD3
        
        rV1 = diag(corr((mVa(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rV2 = diag(corr((mVa(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rV3 = diag(corr((mVa(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rV(:,j,k) = [rV1;rV2;rV3];
        clear rV1 rV2 rV3
        
        rB1 = diag(corr((mba(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rB2 = diag(corr((mba(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rB3 = diag(corr((mba(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rB(:,j,k) = [rB1;rB2;rB3];
        clear rB1 rB2 rB3
        
    end
end

save([dpath 'grid_corrs_fish_climate_FOSI_fished_',mod,'.mat'],...
    'rS','rM','rL','rF','rP','rD','rV','rB','tanom','yr');

%% put together and on grid
for j=1:11
    for k=1:3
        
        Scorr = nan(ni,nj);
        Mcorr = nan(ni,nj);
        Lcorr = nan(ni,nj);
        Fcorr = nan(ni,nj);
        Pcorr = nan(ni,nj);
        Dcorr = nan(ni,nj);
        Vcorr = nan(ni,nj);
        Bcorr = nan(ni,nj);
        
        Scorr(ID) = rS(:,j,k);
        Mcorr(ID) = rM(:,j,k);
        Lcorr(ID) = rL(:,j,k);
        Fcorr(ID) = rF(:,j,k);
        Pcorr(ID) = rP(:,j,k);
        Dcorr(ID) = rD(:,j,k);
        Vcorr(ID) = rV(:,j,k);
        Bcorr(ID) = rB(:,j,k);
        
        %% Fix seam
        [lat_s,lon_s,Scorr2] = cyclic_map_seam_cesm(TLAT,TLONG,Scorr);
        [~,~,Mcorr2] = cyclic_map_seam(TLAT,TLONG,Mcorr);
        [~,~,Lcorr2] = cyclic_map_seam(TLAT,TLONG,Lcorr);
        [~,~,Fcorr2] = cyclic_map_seam(TLAT,TLONG,Fcorr);
        [~,~,Pcorr2] = cyclic_map_seam(TLAT,TLONG,Pcorr);
        [~,~,Dcorr2] = cyclic_map_seam(TLAT,TLONG,Dcorr);
        [~,~,Bcorr2] = cyclic_map_seam(TLAT,TLONG,Bcorr);
        [~,~,Vcorr2] = cyclic_map_seam(TLAT,TLONG,Vcorr);
        
        %% Plots
        close all
        f5 = figure('Units','inches','Position',[1 3 16 8]);
        
        %1,1 - forage
        subplot('Position',[0.01 0.675 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Fcorr2)
        cmocean('balance')
        caxis([-1 1])
        %         text(-2.5,2.25,'a','FontWeight','bold','FontSize',14)
        %         text(0,2.2,'Historic corr2','HorizontalAlignment','center','FontWeight','bold')
        text(0,1.675,'Forage fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %2,1 - all verts
        subplot('Position',[0.16 0.375 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Vcorr2)
        cmocean('balance')
        caxis([-1 1])
        text(0,1.675,'All fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %3,1 - small
        subplot('Position',[0.01 0.075 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Scorr2)
        cmocean('balance')
        caxis([-1 1])
        text(0,1.675,'Small fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %1,2 - lg pel
        subplot('Position',[0.32 0.675 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Pcorr2)
        cmocean('balance')
        caxis([-1 1])
        %         text(-2.5,2.25,'b','FontWeight','bold','FontSize',14)
        text(0,1.85,tanom{j},'HorizontalAlignment','center','FontWeight','bold')
        text(0,1.675,'Large pelagics','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %2,2 -
        %         subplot('Position',[0.32 0.375 0.3 0.275])
        %         axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        %             'Grid','off','FLineWidth',1)
        %         surfm(lat_s,lon_s,pdiff_nc2)
        %         cmocean('balance')
        %         caxis([-1 1])
        %         text(0,1.675,'CNRM','HorizontalAlignment','center')
        %         h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        %
        %3,2 - med
        subplot('Position',[0.32 0.075 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Mcorr2)
        cmocean('balance')
        caxis([-1 1])
        text(0,1.675,'Medium fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %1,3 - dems
        subplot('Position',[0.63 0.675 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Dcorr2)
        cmocean('balance')
        caxis([-1 1])
        %         text(-2.5,2.25,'c','FontWeight','bold','FontSize',14)
        %         text(0,2.2,'% \Delta mesoz','HorizontalAlignment','center','FontWeight','bold')
        text(0,1.675,'Demersals','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %2,3 - benthos
        subplot('Position',[0.475 0.375 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Bcorr2)
        cmocean('balance')
        caxis([-1 1])
        text(0,1.675,'Benthic inverts','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %3,3 - large
        subplot('Position',[0.63 0.075 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Lcorr2)
        cmocean('balance')
        caxis([-1 1])
        text(0,1.675,'Large fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        colorbar('Position',[0.325 0.035 0.3 0.02],'orientation','horizontal');
        
        print('-dpng',[ppath 'Map_climate_corr_FOSI_',mod,tanom{j},'_','lag',num2str(k),'_types.png'])
        
    end
end

%% add yr of lag for bigger
for j=1:11
    for k=1:2
        
        Scorr = nan(ni,nj);
        Mcorr = nan(ni,nj);
        Lcorr = nan(ni,nj);
        Fcorr = nan(ni,nj);
        Pcorr = nan(ni,nj);
        Dcorr = nan(ni,nj);
        Vcorr = nan(ni,nj);
        Bcorr = nan(ni,nj);
        
        Scorr(ID) = rS(:,j,k);
        Mcorr(ID) = rM(:,j,k+1);
        Lcorr(ID) = rL(:,j,k+1);
        Fcorr(ID) = rF(:,j,k);
        Pcorr(ID) = rP(:,j,k+1);
        Dcorr(ID) = rD(:,j,k+1);
        Vcorr(ID) = rV(:,j,k+1);
        Bcorr(ID) = rB(:,j,k);
        
        %% Fix seam
        [lat_s,lon_s,Scorr2] = cyclic_map_seam_cesm(TLAT,TLONG,Scorr);
        [~,~,Mcorr2] = cyclic_map_seam(TLAT,TLONG,Mcorr);
        [~,~,Lcorr2] = cyclic_map_seam(TLAT,TLONG,Lcorr);
        [~,~,Fcorr2] = cyclic_map_seam(TLAT,TLONG,Fcorr);
        [~,~,Pcorr2] = cyclic_map_seam(TLAT,TLONG,Pcorr);
        [~,~,Dcorr2] = cyclic_map_seam(TLAT,TLONG,Dcorr);
        [~,~,Bcorr2] = cyclic_map_seam(TLAT,TLONG,Bcorr);
        [~,~,Vcorr2] = cyclic_map_seam(TLAT,TLONG,Vcorr);
        
        %% Plots
        close all
        f5 = figure('Units','inches','Position',[1 3 16 8]);
        
        %1,1 - forage
        subplot('Position',[0.01 0.675 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Fcorr2)
        cmocean('balance')
        caxis([-1 1])
        %         text(-2.5,2.25,'a','FontWeight','bold','FontSize',14)
        %         text(0,2.2,'Historic corr','HorizontalAlignment','center','FontWeight','bold')
        text(0,1.675,'Forage fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %2,1 - all verts
        subplot('Position',[0.16 0.375 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Vcorr2)
        cmocean('balance')
        caxis([-1 1])
        text(0,1.675,'All fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %3,1 - small
        subplot('Position',[0.01 0.075 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Scorr2)
        cmocean('balance')
        caxis([-1 1])
        text(0,1.675,'Small fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %1,2 - lg pel
        subplot('Position',[0.32 0.675 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Pcorr2)
        cmocean('balance')
        caxis([-1 1])
        %         text(-2.5,2.25,'b','FontWeight','bold','FontSize',14)
        text(0,1.85,tanom{j},'HorizontalAlignment','center','FontWeight','bold')
        text(0,1.675,'Large pelagics','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %2,2 -
        %         subplot('Position',[0.32 0.375 0.3 0.275])
        %         axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        %             'Grid','off','FLineWidth',1)
        %         surfm(lat_s,lon_s,pdiff_nc2)
        %         cmocean('balance')
        %         caxis([-1 1])
        %         text(0,1.675,'CNRM','HorizontalAlignment','center')
        %         h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        %
        %3,2 - med
        subplot('Position',[0.32 0.075 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Mcorr2)
        cmocean('balance')
        caxis([-1 1])
        text(0,1.675,'Medium fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %1,3 - dems
        subplot('Position',[0.63 0.675 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Dcorr2)
        cmocean('balance')
        caxis([-1 1])
        %         text(-2.5,2.25,'c','FontWeight','bold','FontSize',14)
        %         text(0,2.2,'% \Delta mesoz','HorizontalAlignment','center','FontWeight','bold')
        text(0,1.675,'Demersals','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %2,3 - benthos
        subplot('Position',[0.475 0.375 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Bcorr2)
        cmocean('balance')
        caxis([-1 1])
        text(0,1.675,'Benthic inverts','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %3,3 - large
        subplot('Position',[0.63 0.075 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Lcorr2)
        cmocean('balance')
        caxis([-1 1])
        text(0,1.675,'Large fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        colorbar('Position',[0.325 0.035 0.3 0.02],'orientation','horizontal');
        
        print('-dpng',[ppath 'Map_climate_corr_FOSI_',mod,tanom{j},'_','lag',num2str(k),'-',num2str(k+1),'_types.png'])
        
    end
end

%% focus on N Amer & PDO

plotminlat=0; %Set these bounds for your data
plotmaxlat=80;
plotminlon=-190;
plotmaxlon=-20;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

for j=10
    k=1;
    
    Scorr = nan(ni,nj);
    Mcorr = nan(ni,nj);
    Lcorr = nan(ni,nj);
    Fcorr = nan(ni,nj);
    Pcorr = nan(ni,nj);
    Dcorr = nan(ni,nj);
    Vcorr = nan(ni,nj);
    Bcorr = nan(ni,nj);
    
    Scorr(ID) = rS(:,j,k);
    Mcorr(ID) = rM(:,j,k+1);
    Lcorr(ID) = rL(:,j,k+2);
    Fcorr(ID) = rF(:,j,k);
    Pcorr(ID) = rP(:,j,k+1);
    Dcorr(ID) = rD(:,j,k+1);
    Vcorr(ID) = rV(:,j,k+1);
    Bcorr(ID) = rB(:,j,k);
    
    %% Fix seam
    [lat_s,lon_s,Scorr2] = cyclic_map_seam_cesm(TLAT,TLONG,Scorr);
    [~,~,Mcorr2] = cyclic_map_seam(TLAT,TLONG,Mcorr);
    [~,~,Lcorr2] = cyclic_map_seam(TLAT,TLONG,Lcorr);
    [~,~,Fcorr2] = cyclic_map_seam(TLAT,TLONG,Fcorr);
    [~,~,Pcorr2] = cyclic_map_seam(TLAT,TLONG,Pcorr);
    [~,~,Dcorr2] = cyclic_map_seam(TLAT,TLONG,Dcorr);
    [~,~,Bcorr2] = cyclic_map_seam(TLAT,TLONG,Bcorr);
    [~,~,Vcorr2] = cyclic_map_seam(TLAT,TLONG,Vcorr);
    
    %% Plots
    close all
    f5 = figure('Units','inches','Position',[1 3 16 8]);
    
    %1,1 - forage
    subplot('Position',[0.01 0.675 0.3 0.275])
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Fcorr2)
    cmocean('balance')
    caxis([-1 1])
    %         text(-2.5,2.25,'a','FontWeight','bold','FontSize',14)
    %         text(0,2.2,'Historic corr','HorizontalAlignment','center','FontWeight','bold')
    text(0,1.925,'Forage fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2,1 - all verts
    subplot('Position',[0.16 0.375 0.3 0.275])
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Vcorr2)
    cmocean('balance')
    caxis([-1 1])
    text(0,1.925,'All fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3,1 - small
    subplot('Position',[0.01 0.075 0.3 0.275])
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr2)
    cmocean('balance')
    caxis([-1 1])
    text(0,1.925,'Small fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %1,2 - lg pel
    subplot('Position',[0.32 0.675 0.3 0.275])
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Pcorr2)
    cmocean('balance')
    caxis([-1 1])
    %         text(-2.5,2.25,'b','FontWeight','bold','FontSize',14)
    text(0,2.15,tanom{j},'HorizontalAlignment','center','FontWeight','bold')
    text(0,1.925,'Large pelagics','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2,2 -
    %         subplot('Position',[0.32 0.375 0.3 0.275])
    %         axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    %             'Grid','off','FLineWidth',1)
    %         surfm(lat_s,lon_s,pdiff_nc2)
    %         cmocean('balance')
    %         caxis([-1 1])
    %         text(0,1.925,'CNRM','HorizontalAlignment','center')
    %         h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    %
    %3,2 - med
    subplot('Position',[0.32 0.075 0.3 0.275])
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr2)
    cmocean('balance')
    caxis([-1 1])
    text(0,1.925,'Medium fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %1,3 - dems
    subplot('Position',[0.63 0.675 0.3 0.275])
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Dcorr2)
    cmocean('balance')
    caxis([-1 1])
    %         text(-2.5,2.25,'c','FontWeight','bold','FontSize',14)
    %         text(0,2.2,'% \Delta mesoz','HorizontalAlignment','center','FontWeight','bold')
    text(0,1.925,'Demersals','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2,3 - benthos
    subplot('Position',[0.475 0.375 0.3 0.275])
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Bcorr2)
    cmocean('balance')
    caxis([-1 1])
    text(0,1.925,'Benthic inverts','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3,3 - large
    subplot('Position',[0.63 0.075 0.3 0.275])
    axesm ('Miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr2)
    cmocean('balance')
    caxis([-1 1])
    text(0,1.925,'Large fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    colorbar('Position',[0.325 0.035 0.3 0.02],'orientation','horizontal');
    
    print('-dpng',[ppath 'Map_climate_corr_FOSI_NAmer_',mod,tanom{j},'_','lag',num2str(k),'-',num2str(k+1),'_types.png'])
    
    
end


