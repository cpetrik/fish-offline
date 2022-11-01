% Map linear regressions of FEISTY FOSI inputs & fish
% w/climate indices

clear all
close all

%%
apath = '/Users/cpetrik/Dropbox/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'Climate_anomalies_annual_means.mat'],'canom');

tanom = canom;
clear canom

%% Map data
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
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

cmRB=flipud(cbrewer('div','RdBu',66,'pchip'));

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
ppath = [pp cfile '/corrs/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

load([fpath 'grid_linreg_climate_inputs_fish_FOSI_fished_',mod,'.mat'])

%% All types together with same lag
yr = 0:4; %lags

for j=10%1:11
    for k=1:length(yr) %Linear regression at diff lags
        t = yr(k);
        
        %Coeff of linear regression
        Scorr  = nan(ni,nj);
        Mcorr  = nan(ni,nj);
        Lcorr  = nan(ni,nj);
        Fcorr  = nan(ni,nj);
        Pcorr  = nan(ni,nj);
        Dcorr  = nan(ni,nj);
        Vcorr  = nan(ni,nj);
        Bcorr  = nan(ni,nj);
        TPcorr = nan(ni,nj);
        TBcorr = nan(ni,nj);
        DETcorr = nan(ni,nj);
        ZBcorr  = nan(ni,nj);
        ZLcorr  = nan(ni,nj);
        
        Scorr(ID)  = rS(:,j,k);
        Mcorr(ID)  = rM(:,j,k);
        Lcorr(ID)  = rL(:,j,k);
        Fcorr(ID)  = rF(:,j,k);
        Pcorr(ID)  = rP(:,j,k);
        Dcorr(ID)  = rD(:,j,k);
        Vcorr(ID)  = rV(:,j,k);
        Bcorr(ID)  = rB(:,j,k);
        TPcorr(ID) = rTp(:,j,k);
        TBcorr(ID) = rTb(:,j,k);
        DETcorr(ID) = rDet(:,j,k);
        ZBcorr(ID)  = rZb(:,j,k);
        ZLcorr(ID)  = rZl(:,j,k);
        
        %significance
        Ssig  = ones(ni,nj);
        Msig  = ones(ni,nj);
        Lsig  = ones(ni,nj);
        Fsig  = ones(ni,nj);
        Psig  = ones(ni,nj);
        Dsig  = ones(ni,nj);
        Vsig  = ones(ni,nj);
        Bsig  = ones(ni,nj);
        TPsig = ones(ni,nj);
        TBsig = ones(ni,nj);
        DETsig = ones(ni,nj);
        ZBsig  = ones(ni,nj);
        ZLsig  = ones(ni,nj);
        
        Ssig(ID)  = pS(:,j,k);
        Msig(ID)  = pM(:,j,k);
        Lsig(ID)  = pL(:,j,k);
        Fsig(ID)  = pF(:,j,k);
        Psig(ID)  = pP(:,j,k);
        Dsig(ID)  = pD(:,j,k);
        Vsig(ID)  = pV(:,j,k);
        Bsig(ID)  = pB(:,j,k);
        TPsig(ID) = pTp(:,j,k);
        TBsig(ID) = pTb(:,j,k);
        DETsig(ID) = pDet(:,j,k);
        ZBsig(ID)  = pZb(:,j,k);
        ZLsig(ID)  = pZl(:,j,k);
        
        Ssig   = (Ssig<0.05);
        Msig   = (Msig<0.05);
        Lsig   = (Lsig<0.05);
        Fsig   = (Fsig<0.05);
        Psig   = (Psig<0.05);
        Dsig   = (Dsig<0.05);
        Vsig   = (Vsig<0.05);
        Bsig   = (Bsig<0.05);
        TPsig  = (TPsig<0.05);
        TBsig  = (TBsig<0.05);
        DETsig = (DETsig<0.05);
        ZBsig  = (ZBsig<0.05);
        ZLsig  = (ZLsig<0.05);
        
        %% Fix seam
        [lat_s,lon_s,Scorr2] = cyclic_map_seam_cesm(TLAT,TLONG,Scorr);
        [~,~,Mcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Mcorr);
        [~,~,Lcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Lcorr);
        [~,~,Fcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Fcorr);
        [~,~,Pcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Pcorr);
        [~,~,Dcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Dcorr);
        [~,~,Bcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Bcorr);
        [~,~,Vcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Vcorr);
        [~,~,TPcorr2]  = cyclic_map_seam_cesm(TLAT,TLONG,TPcorr);
        [~,~,TBcorr2]  = cyclic_map_seam_cesm(TLAT,TLONG,TBcorr);
        [~,~,DETcorr2] = cyclic_map_seam_cesm(TLAT,TLONG,DETcorr);
        [~,~,ZBcorr2]  = cyclic_map_seam_cesm(TLAT,TLONG,ZBcorr);
        [~,~,ZLcorr2]  = cyclic_map_seam_cesm(TLAT,TLONG,ZLcorr);
        
        %% Plots
        close all
        f1 = figure('Units','inches','Position',[1 3 16 8]);
        
        %1,1 - forage
        subplot('Position',[0.01 0.675 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Fcorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,Fsig,'markersize',3,'density',150)
        text(0,1.675,'Forage fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %2,1 - all verts
        subplot('Position',[0.16 0.375 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Vcorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,Vsig,'markersize',3,'density',150)
        text(0,1.675,'All fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %3,1 - small
        subplot('Position',[0.01 0.075 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Scorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,Ssig,'markersize',3,'density',150)
        text(0,1.675,'Small fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %1,2 - lg pel
        subplot('Position',[0.32 0.675 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Pcorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,Psig,'markersize',3,'density',150)
        text(0,1.85,[tanom{j} ' lag ',num2str(t) ' yr'],'HorizontalAlignment','center','FontWeight','bold')
        text(0,1.675,'Large pelagics','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %2,2 -
        %         subplot('Position',[0.32 0.375 0.3 0.275])
        %         axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        %             'Grid','off','FLineWidth',1)
        %         surfm(lat_s,lon_s,pdiff_nc2)
        %         colormap(cmRB)
        %         caxis([-0.7 0.7])
        %         text(0,1.675,'CNRM','HorizontalAlignment','center')
        %         h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        %
        %3,2 - med
        subplot('Position',[0.32 0.075 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Mcorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,Msig,'markersize',3,'density',150)
        text(0,1.675,'Medium fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %1,3 - dems
        subplot('Position',[0.63 0.675 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Dcorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,Dsig,'markersize',3,'density',150)
        text(0,1.675,'Demersals','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %2,3 - benthos
        subplot('Position',[0.475 0.375 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Bcorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,Bsig,'markersize',3,'density',150)
        text(0,1.675,'Benthic inverts','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %3,3 - large
        subplot('Position',[0.63 0.075 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,Lcorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,Lsig,'markersize',3,'density',150)
        text(0,1.675,'Large fish','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        colorbar('Position',[0.325 0.035 0.3 0.02],'orientation','horizontal');
        print('-dpng',[ppath 'Map_climate_linereg_FOSI_',mod,'_',tanom{j},'_','lag',num2str(t),'_fishtypes.png'])
        
        
        %% INPUTS
        f2 = figure('Units','inches','Position',[1 3 16 8]);
        
        %1,1 - TP
        subplot('Position',[0.01 0.675 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,TPcorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,TPsig,'markersize',3,'density',150)
        text(0,1.675,'T pelagic','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %2,1 - TB
        subplot('Position',[0.16 0.375 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,TBcorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,TBsig,'markersize',3,'density',150)
        text(0,1.675,'T bottom','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %1,2 - ZB
        subplot('Position',[0.32 0.675 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,ZBcorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,ZBsig,'markersize',3,'density',150)
        text(0,1.85,[tanom{j} ' lag ',num2str(t) ' yr'],'HorizontalAlignment','center','FontWeight','bold')
        text(0,1.675,'Mesoz biom','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %1,3 - ZL
        subplot('Position',[0.63 0.675 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,ZLcorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,ZLsig,'markersize',3,'density',150)
        text(0,1.675,'Mesoz loss','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        
        %2,3 - det
        subplot('Position',[0.475 0.375 0.3 0.275])
        axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(lat_s,lon_s,DETcorr2)
        colormap(cmRB)
        caxis([-0.7 0.7])
        stipplem(TLAT,TLONG,DETsig,'markersize',3,'density',150)
        text(0,1.675,'Detritus','HorizontalAlignment','center')
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        colorbar('Position',[0.325 0.335 0.3 0.02],'orientation','horizontal');
        print('-dpng',[ppath 'Map_climate_linereg_FOSI_',mod,'_',tanom{j},'_','lag',num2str(t),'_inputs.png'])
        
        
    end
end


%% focus on 3 types cols, rows of lags -------------------------------
for j=1:11
    close all
    
    %% S,M,L
    ScorrK  = nan(ni,nj,length(yr));
    McorrK  = nan(ni,nj,length(yr));
    LcorrK  = nan(ni,nj,length(yr));
    
    SsigK  = ones(ni,nj,length(yr));
    MsigK  = ones(ni,nj,length(yr));
    LsigK  = ones(ni,nj,length(yr));
    
    for k=1:length(yr) %Linear regression at diff lags
        
        %Coeff of linear regression
        Scorr  = nan(ni,nj);
        Mcorr  = nan(ni,nj);
        Lcorr  = nan(ni,nj);
        
        Scorr(ID)  = rS(:,j,k);
        Mcorr(ID)  = rM(:,j,k);
        Lcorr(ID)  = rL(:,j,k);
        
        ScorrK(:,:,k) = Scorr;
        McorrK(:,:,k) = Mcorr;
        LcorrK(:,:,k) = Lcorr;
        
        %significance
        Ssig  = ones(ni,nj);
        Msig  = ones(ni,nj);
        Lsig  = ones(ni,nj);
        
        Ssig(ID)  = pS(:,j,k);
        Msig(ID)  = pM(:,j,k);
        Lsig(ID)  = pL(:,j,k);
        
        SsigK(:,:,k) = Ssig;
        MsigK(:,:,k) = Msig;
        LsigK(:,:,k) = Lsig;
    end
    
    SsigK   = (SsigK<0.05);
    MsigK   = (MsigK<0.05);
    LsigK   = (LsigK<0.05);
    
    clear Scorr Mcorr Lcorr Ssig Msig Lsig
        
    %% Fix seam
    [lat_s,lon_s,Scorr0] = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,1));
    [~,~,Scorr1]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,2));
    [~,~,Scorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,3));
    [~,~,Scorr3]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,4));
    [~,~,Scorr4]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,5));
    [~,~,Mcorr0]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,1));
    [~,~,Mcorr1]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,2));
    [~,~,Mcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,3));
    [~,~,Mcorr3]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,4));
    [~,~,Mcorr4]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,5));
    [~,~,Lcorr0]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,1));
    [~,~,Lcorr1]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,2));
    [~,~,Lcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,3));
    [~,~,Lcorr3]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,4));
    [~,~,Lcorr4]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,5));
    
    %%
    f1 = figure('Units','inches','Position',[1 3 7.5 10]);
    % Small
    %1 -
    subplot('Position',[0.01 0.78 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr0)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,1)),'markersize',3,'density',75)
    text(0,2.25,'Small fish','HorizontalAlignment','center','FontWeight','bold')
    text(-1.75,1.75,'0 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2 -
    subplot('Position',[0.01 0.605 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr1)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,2)),'markersize',3,'density',75)
    text(-1.75,1.75,'1 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3 -
    subplot('Position',[0.01 0.43 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr2)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,3)),'markersize',3,'density',75)
    text(-1.75,1.75,'2 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %4 -
    subplot('Position',[0.01 0.255 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr3)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,4)),'markersize',3,'density',75)
    text(-1.75,1.75,'3 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %5 -
    subplot('Position',[0.01 0.08 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr4)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,5)),'markersize',3,'density',75)
    text(-1.75,1.75,'4 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colorbar('Position',[0.25 0.05 0.5 0.02],'orientation','horizontal');
    
    % Medium
    %1 -
    subplot('Position',[0.33 0.78 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr0)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,1)),'markersize',3,'density',75)
    text(0,2.75,tanom{j},'HorizontalAlignment','center','FontWeight','bold')
    text(0,2.25,'Medium fish','HorizontalAlignment','center','FontWeight','bold')
    text(-1.95,1.75,'0 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2 -
    subplot('Position',[0.33 0.605 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr1)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,2)),'markersize',3,'density',75)
    text(-1.95,1.75,'1 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3 -
    subplot('Position',[0.33 0.43 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr2)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,3)),'markersize',3,'density',75)
    text(-1.95,1.75,'2 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %4 -
    subplot('Position',[0.33 0.255 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr3)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,4)),'markersize',3,'density',75)
    text(-1.95,1.75,'3 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %5 -
    subplot('Position',[0.33 0.08 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr4)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,5)),'markersize',3,'density',75)
    text(-1.95,1.75,'4 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    % Large
    % 1 -
    subplot('Position',[0.65 0.78 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr0)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,1)),'markersize',3,'density',75)
    text(0,2.25,'Large fish','HorizontalAlignment','center','FontWeight','bold')
    text(-1.95,1.75,'0 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2 -
    subplot('Position',[0.65 0.605 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr1)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,2)),'markersize',3,'density',75)
    text(-1.95,1.75,'1 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3 -
    subplot('Position',[0.65 0.43 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr2)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,3)),'markersize',3,'density',75)
    text(-1.95,1.75,'2 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %4 -
    subplot('Position',[0.65 0.255 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr3)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,4)),'markersize',3,'density',75)
    text(-1.95,1.75,'3 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %5 -
    subplot('Position',[0.65 0.08 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr4)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,5)),'markersize',3,'density',75)
    text(-1.95,1.75,'4 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    print('-dpng',[ppath 'Map_climate_linereg_FOSI_',mod,'_',tanom{j},'_fish_sizes.png'])
        
    
    %% F, P, D
    ScorrK  = nan(ni,nj,length(yr));
    McorrK  = nan(ni,nj,length(yr));
    LcorrK  = nan(ni,nj,length(yr));
    
    SsigK  = ones(ni,nj,length(yr));
    MsigK  = ones(ni,nj,length(yr));
    LsigK  = ones(ni,nj,length(yr));
    
    for k=1:length(yr) %Linear regression at diff lags
        
        %Coeff of linear regression
        Scorr  = nan(ni,nj);
        Mcorr  = nan(ni,nj);
        Lcorr  = nan(ni,nj);
        
        Scorr(ID)  = rF(:,j,k);
        Mcorr(ID)  = rP(:,j,k);
        Lcorr(ID)  = rD(:,j,k);
        
        ScorrK(:,:,k) = Scorr;
        McorrK(:,:,k) = Mcorr;
        LcorrK(:,:,k) = Lcorr;
        
        %significance
        Ssig  = ones(ni,nj);
        Msig  = ones(ni,nj);
        Lsig  = ones(ni,nj);
        
        Ssig(ID)  = pF(:,j,k);
        Msig(ID)  = pP(:,j,k);
        Lsig(ID)  = pD(:,j,k);
        
        SsigK(:,:,k) = Ssig;
        MsigK(:,:,k) = Msig;
        LsigK(:,:,k) = Lsig;
    end
    
    SsigK   = (SsigK<0.05);
    MsigK   = (MsigK<0.05);
    LsigK   = (LsigK<0.05);
    
    clear Scorr Mcorr Lcorr Ssig Msig Lsig
        
    %% Fix seam
    [lat_s,lon_s,Scorr0] = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,1));
    [~,~,Scorr1]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,2));
    [~,~,Scorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,3));
    [~,~,Scorr3]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,4));
    [~,~,Scorr4]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,5));
    [~,~,Mcorr0]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,1));
    [~,~,Mcorr1]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,2));
    [~,~,Mcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,3));
    [~,~,Mcorr3]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,4));
    [~,~,Mcorr4]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,5));
    [~,~,Lcorr0]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,1));
    [~,~,Lcorr1]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,2));
    [~,~,Lcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,3));
    [~,~,Lcorr3]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,4));
    [~,~,Lcorr4]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,5));
    
    %%
    f2 = figure('Units','inches','Position',[1 3 7.5 10]);
    % Small
    %1 -
    subplot('Position',[0.01 0.78 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr0)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,1)),'markersize',3,'density',75)
    text(0,2.25,'Forage','HorizontalAlignment','center','FontWeight','bold')
    text(-1.75,1.75,'0 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2 -
    subplot('Position',[0.01 0.605 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr1)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,2)),'markersize',3,'density',75)
    text(-1.75,1.75,'1 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3 -
    subplot('Position',[0.01 0.43 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr2)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,3)),'markersize',3,'density',75)
    text(-1.75,1.75,'2 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %4 -
    subplot('Position',[0.01 0.255 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr3)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,4)),'markersize',3,'density',75)
    text(-1.75,1.75,'3 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %5 -
    subplot('Position',[0.01 0.08 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr4)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,5)),'markersize',3,'density',75)
    text(-1.75,1.75,'4 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colorbar('Position',[0.25 0.05 0.5 0.02],'orientation','horizontal');
    
    % Medium
    %1 -
    subplot('Position',[0.33 0.78 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr0)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,1)),'markersize',3,'density',75)
    text(0,2.75,tanom{j},'HorizontalAlignment','center','FontWeight','bold')
    text(0,2.25,'Large pelagics','HorizontalAlignment','center','FontWeight','bold')
    text(-1.95,1.75,'0 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2 -
    subplot('Position',[0.33 0.605 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr1)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,2)),'markersize',3,'density',75)
    text(-1.95,1.75,'1 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3 -
    subplot('Position',[0.33 0.43 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr2)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,3)),'markersize',3,'density',75)
    text(-1.95,1.75,'2 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %4 -
    subplot('Position',[0.33 0.255 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr3)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,4)),'markersize',3,'density',75)
    text(-1.95,1.75,'3 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %5 -
    subplot('Position',[0.33 0.08 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr4)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,5)),'markersize',3,'density',75)
    text(-1.95,1.75,'4 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    % Large
    % 1 -
    subplot('Position',[0.65 0.78 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr0)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,1)),'markersize',3,'density',75)
    text(0,2.25,'Demersals','HorizontalAlignment','center','FontWeight','bold')
    text(-1.95,1.75,'0 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2 -
    subplot('Position',[0.65 0.605 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr1)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,2)),'markersize',3,'density',75)
    text(-1.95,1.75,'1 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3 -
    subplot('Position',[0.65 0.43 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr2)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,3)),'markersize',3,'density',75)
    text(-1.95,1.75,'2 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %4 -
    subplot('Position',[0.65 0.255 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr3)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,4)),'markersize',3,'density',75)
    text(-1.95,1.75,'3 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %5 -
    subplot('Position',[0.65 0.08 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr4)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,5)),'markersize',3,'density',75)
    text(-1.95,1.75,'4 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    print('-dpng',[ppath 'Map_climate_linereg_FOSI_',mod,'_',tanom{j},'_fish_types.png'])
       
    
    %% TP, Zb, Zl
    ScorrK  = nan(ni,nj,length(yr));
    McorrK  = nan(ni,nj,length(yr));
    LcorrK  = nan(ni,nj,length(yr));
    
    SsigK  = ones(ni,nj,length(yr));
    MsigK  = ones(ni,nj,length(yr));
    LsigK  = ones(ni,nj,length(yr));
    
    for k=1:length(yr) %Linear regression at diff lags
        
        %Coeff of linear regression
        Scorr  = nan(ni,nj);
        Mcorr  = nan(ni,nj);
        Lcorr  = nan(ni,nj);
        
        Scorr(ID)  = rTp(:,j,k);
        Mcorr(ID)  = rZb(:,j,k);
        Lcorr(ID)  = rZl(:,j,k);
        
        ScorrK(:,:,k) = Scorr;
        McorrK(:,:,k) = Mcorr;
        LcorrK(:,:,k) = Lcorr;
        
        %significance
        Ssig  = ones(ni,nj);
        Msig  = ones(ni,nj);
        Lsig  = ones(ni,nj);
        
        Ssig(ID)  = pTp(:,j,k);
        Msig(ID)  = pZb(:,j,k);
        Lsig(ID)  = pZl(:,j,k);
        
        SsigK(:,:,k) = Ssig;
        MsigK(:,:,k) = Msig;
        LsigK(:,:,k) = Lsig;
    end
    
    SsigK   = (SsigK<0.05);
    MsigK   = (MsigK<0.05);
    LsigK   = (LsigK<0.05);
    
    clear Scorr Mcorr Lcorr Ssig Msig Lsig
        
    %% Fix seam
    [lat_s,lon_s,Scorr0] = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,1));
    [~,~,Scorr1]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,2));
    [~,~,Scorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,3));
    [~,~,Scorr3]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,4));
    [~,~,Scorr4]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,5));
    [~,~,Mcorr0]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,1));
    [~,~,Mcorr1]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,2));
    [~,~,Mcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,3));
    [~,~,Mcorr3]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,4));
    [~,~,Mcorr4]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,5));
    [~,~,Lcorr0]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,1));
    [~,~,Lcorr1]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,2));
    [~,~,Lcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,3));
    [~,~,Lcorr3]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,4));
    [~,~,Lcorr4]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,5));
    
    %%
    f3 = figure('Units','inches','Position',[1 3 7.5 10]);
    % Small
    %1 -
    subplot('Position',[0.01 0.78 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr0)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,1)),'markersize',3,'density',75)
    text(0,2.25,'T pel','HorizontalAlignment','center','FontWeight','bold')
    text(-1.75,1.75,'0 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2 -
    subplot('Position',[0.01 0.605 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr1)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,2)),'markersize',3,'density',75)
    text(-1.75,1.75,'1 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3 -
    subplot('Position',[0.01 0.43 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr2)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,3)),'markersize',3,'density',75)
    text(-1.75,1.75,'2 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %4 -
    subplot('Position',[0.01 0.255 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr3)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,4)),'markersize',3,'density',75)
    text(-1.75,1.75,'3 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %5 -
    subplot('Position',[0.01 0.08 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr4)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,5)),'markersize',3,'density',75)
    text(-1.75,1.75,'4 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colorbar('Position',[0.25 0.05 0.5 0.02],'orientation','horizontal');
    
    % Medium
    %1 -
    subplot('Position',[0.33 0.78 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr0)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,1)),'markersize',3,'density',75)
    text(0,2.75,tanom{j},'HorizontalAlignment','center','FontWeight','bold')
    text(0,2.25,'Mesoz biom','HorizontalAlignment','center','FontWeight','bold')
    text(-1.95,1.75,'0 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2 -
    subplot('Position',[0.33 0.605 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr1)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,2)),'markersize',3,'density',75)
    text(-1.95,1.75,'1 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3 -
    subplot('Position',[0.33 0.43 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr2)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,3)),'markersize',3,'density',75)
    text(-1.95,1.75,'2 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %4 -
    subplot('Position',[0.33 0.255 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr3)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,4)),'markersize',3,'density',75)
    text(-1.95,1.75,'3 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %5 -
    subplot('Position',[0.33 0.08 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr4)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,5)),'markersize',3,'density',75)
    text(-1.95,1.75,'4 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    % Large
    % 1 -
    subplot('Position',[0.65 0.78 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr0)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,1)),'markersize',3,'density',75)
    text(0,2.25,'Mesoz loss','HorizontalAlignment','center','FontWeight','bold')
    text(-1.95,1.75,'0 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2 -
    subplot('Position',[0.65 0.605 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr1)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,2)),'markersize',3,'density',75)
    text(-1.95,1.75,'1 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3 -
    subplot('Position',[0.65 0.43 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr2)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,3)),'markersize',3,'density',75)
    text(-1.95,1.75,'2 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %4 -
    subplot('Position',[0.65 0.255 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr3)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,4)),'markersize',3,'density',75)
    text(-1.95,1.75,'3 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %5 -
    subplot('Position',[0.65 0.08 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr4)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,5)),'markersize',3,'density',75)
    text(-1.95,1.75,'4 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    print('-dpng',[ppath 'Map_climate_linereg_FOSI_',mod,'_',tanom{j},'_pel.png'])
       
    
    %% TB, Det, Bent
    ScorrK  = nan(ni,nj,length(yr));
    McorrK  = nan(ni,nj,length(yr));
    LcorrK  = nan(ni,nj,length(yr));
    
    SsigK  = ones(ni,nj,length(yr));
    MsigK  = ones(ni,nj,length(yr));
    LsigK  = ones(ni,nj,length(yr));
    
    for k=1:length(yr) %Linear regression at diff lags
        
        %Coeff of linear regression
        Scorr  = nan(ni,nj);
        Mcorr  = nan(ni,nj);
        Lcorr  = nan(ni,nj);
        
        Scorr(ID)  = rTb(:,j,k);
        Mcorr(ID)  = rDet(:,j,k);
        Lcorr(ID)  = rB(:,j,k);
        
        ScorrK(:,:,k) = Scorr;
        McorrK(:,:,k) = Mcorr;
        LcorrK(:,:,k) = Lcorr;
        
        %significance
        Ssig  = ones(ni,nj);
        Msig  = ones(ni,nj);
        Lsig  = ones(ni,nj);
        
        Ssig(ID)  = pTb(:,j,k);
        Msig(ID)  = pDet(:,j,k);
        Lsig(ID)  = pB(:,j,k);
        
        SsigK(:,:,k) = Ssig;
        MsigK(:,:,k) = Msig;
        LsigK(:,:,k) = Lsig;
    end
    
    SsigK   = (SsigK<0.05);
    MsigK   = (MsigK<0.05);
    LsigK   = (LsigK<0.05);
    
    clear Scorr Mcorr Lcorr Ssig Msig Lsig
        
    %% Fix seam
    [lat_s,lon_s,Scorr0] = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,1));
    [~,~,Scorr1]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,2));
    [~,~,Scorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,3));
    [~,~,Scorr3]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,4));
    [~,~,Scorr4]   = cyclic_map_seam_cesm(TLAT,TLONG,ScorrK(:,:,5));
    [~,~,Mcorr0]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,1));
    [~,~,Mcorr1]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,2));
    [~,~,Mcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,3));
    [~,~,Mcorr3]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,4));
    [~,~,Mcorr4]   = cyclic_map_seam_cesm(TLAT,TLONG,McorrK(:,:,5));
    [~,~,Lcorr0]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,1));
    [~,~,Lcorr1]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,2));
    [~,~,Lcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,3));
    [~,~,Lcorr3]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,4));
    [~,~,Lcorr4]   = cyclic_map_seam_cesm(TLAT,TLONG,LcorrK(:,:,5));
    
    %%
    f4 = figure('Units','inches','Position',[1 3 7.5 10]);
    % Small
    %1 -
    subplot('Position',[0.01 0.78 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr0)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,1)),'markersize',3,'density',75)
    text(0,2.25,'T bottom','HorizontalAlignment','center','FontWeight','bold')
    text(-1.75,1.75,'0 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2 -
    subplot('Position',[0.01 0.605 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr1)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,2)),'markersize',3,'density',75)
    text(-1.75,1.75,'1 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3 -
    subplot('Position',[0.01 0.43 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr2)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,3)),'markersize',3,'density',75)
    text(-1.75,1.75,'2 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %4 -
    subplot('Position',[0.01 0.255 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr3)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,4)),'markersize',3,'density',75)
    text(-1.75,1.75,'3 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %5 -
    subplot('Position',[0.01 0.08 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr4)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(SsigK(:,:,5)),'markersize',3,'density',75)
    text(-1.75,1.75,'4 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colorbar('Position',[0.25 0.05 0.5 0.02],'orientation','horizontal');
    
    % Medium
    %1 -
    subplot('Position',[0.33 0.78 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr0)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,1)),'markersize',3,'density',75)
    text(0,2.75,tanom{j},'HorizontalAlignment','center','FontWeight','bold')
    text(0,2.25,'Detritus','HorizontalAlignment','center','FontWeight','bold')
    text(-1.95,1.75,'0 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2 -
    subplot('Position',[0.33 0.605 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr1)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,2)),'markersize',3,'density',75)
    text(-1.95,1.75,'1 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3 -
    subplot('Position',[0.33 0.43 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr2)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,3)),'markersize',3,'density',75)
    text(-1.95,1.75,'2 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %4 -
    subplot('Position',[0.33 0.255 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr3)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,4)),'markersize',3,'density',75)
    text(-1.95,1.75,'3 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %5 -
    subplot('Position',[0.33 0.08 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr4)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(MsigK(:,:,5)),'markersize',3,'density',75)
    text(-1.95,1.75,'4 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    % Large
    % 1 -
    subplot('Position',[0.65 0.78 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr0)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,1)),'markersize',3,'density',75)
    text(0,2.25,'Benthos','HorizontalAlignment','center','FontWeight','bold')
    text(-1.95,1.75,'0 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %2 -
    subplot('Position',[0.65 0.605 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr1)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,2)),'markersize',3,'density',75)
    text(-1.95,1.75,'1 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %3 -
    subplot('Position',[0.65 0.43 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr2)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,3)),'markersize',3,'density',75)
    text(-1.95,1.75,'2 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %4 -
    subplot('Position',[0.65 0.255 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr3)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,4)),'markersize',3,'density',75)
    text(-1.95,1.75,'3 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    
    %5 -
    subplot('Position',[0.65 0.08 0.3 0.175])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr4)
    colormap(cmRB)
    caxis([-0.6 0.6])
    stipplem(TLAT,TLONG,squeeze(LsigK(:,:,5)),'markersize',3,'density',75)
    text(-1.95,1.75,'4 yr lag','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    print('-dpng',[ppath 'Map_climate_linereg_FOSI_',mod,'_',tanom{j},'_bent.png'])
       
    
end

%%

