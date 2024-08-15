% Map correlations of FEISTY FOSI
% w/ experiments

clear
close all


%% Map data
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

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
ppath = [pp cfile '/corrs/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

load([fpath 'grid_corrs_fish_FOSI_fished_expers.mat']);


%% Plots
close all

for j=2:4
    
    %%
    mod = sims{j};

    %Corr
    Scorr  = nan(ni,nj);
    Mcorr  = nan(ni,nj);
    Lcorr  = nan(ni,nj);
    Fcorr  = nan(ni,nj);
    Pcorr  = nan(ni,nj);
    Dcorr  = nan(ni,nj);
    Vcorr  = nan(ni,nj);
    Bcorr  = nan(ni,nj);

    Scorr(ID)  = rS(:,j);
    Mcorr(ID)  = rM(:,j);
    Lcorr(ID)  = rL(:,j);
    Fcorr(ID)  = rF(:,j);
    Pcorr(ID)  = rP(:,j);
    Dcorr(ID)  = rD(:,j);
    Vcorr(ID)  = rA(:,j);
    Bcorr(ID)  = rB(:,j);

    %significance
    Ssig  = ones(ni,nj);
    Msig  = ones(ni,nj);
    Lsig  = ones(ni,nj);
    Fsig  = ones(ni,nj);
    Psig  = ones(ni,nj);
    Dsig  = ones(ni,nj);
    Vsig  = ones(ni,nj);
    Bsig  = ones(ni,nj);

    Ssig(ID)  = pS(:,j);
    Msig(ID)  = pM(:,j);
    Lsig(ID)  = pL(:,j);
    Fsig(ID)  = pF(:,j);
    Psig(ID)  = pP(:,j);
    Dsig(ID)  = pD(:,j);
    Vsig(ID)  = pA(:,j);
    Bsig(ID)  = pB(:,j);

    Ssig   = (Ssig<0.05);
    Msig   = (Msig<0.05);
    Lsig   = (Lsig<0.05);
    Fsig   = (Fsig<0.05);
    Psig   = (Psig<0.05);
    Dsig   = (Dsig<0.05);
    Vsig   = (Vsig<0.05);
    Bsig   = (Bsig<0.05);

    %% Fix seam
    [lat_s,lon_s,Scorr2] = cyclic_map_seam_cesm(TLAT,TLONG,Scorr);
    [~,~,Mcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Mcorr);
    [~,~,Lcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Lcorr);
    [~,~,Fcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Fcorr);
    [~,~,Pcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Pcorr);
    [~,~,Dcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Dcorr);
    [~,~,Bcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Bcorr);
    [~,~,Vcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Vcorr);


    %% Figure by size & type
    f5 = figure('Units','inches','Position',[1 3 16 8]);

    %1,1 - forage
    subplot('Position',[0.01 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Fcorr2)
    cmocean('balance')
    caxis([-1 1])
    stipplem(TLAT,TLONG,Fsig,'markersize',3,'density',150)
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
    stipplem(TLAT,TLONG,Vsig,'markersize',3,'density',150)
        text(0,1.675,'All fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,1 - small
    subplot('Position',[0.01 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr2)
    cmocean('balance')
    caxis([-1 1])
    stipplem(TLAT,TLONG,Ssig,'markersize',3,'density',150)
        text(0,1.675,'Small fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %1,2 - lg pel
    subplot('Position',[0.32 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Pcorr2)
    cmocean('balance')
    caxis([-1 1])
    stipplem(TLAT,TLONG,Psig,'markersize',3,'density',150)
        %         text(-2.5,2.25,'b','FontWeight','bold','FontSize',14)
    text(0,1.85,mod,'HorizontalAlignment','center','FontWeight','bold')
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
    stipplem(TLAT,TLONG,Msig,'markersize',3,'density',150)
        text(0,1.675,'Medium fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %1,3 - dems
    subplot('Position',[0.63 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Dcorr2)
    cmocean('balance')
    caxis([-1 1])
    stipplem(TLAT,TLONG,Dsig,'markersize',3,'density',150)
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
    stipplem(TLAT,TLONG,Bsig,'markersize',3,'density',150)
        text(0,1.675,'Benthic inverts','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,3 - large
    subplot('Position',[0.63 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr2)
    cmocean('balance')
    caxis([-1 1])
    stipplem(TLAT,TLONG,Lsig,'markersize',3,'density',150)
    text(0,1.675,'Large fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    colorbar('Position',[0.325 0.035 0.3 0.02],'orientation','horizontal');

    print('-dpng',[ppath 'Map_corr_FOSI_full_',mod,'types.png'])

end

%% Calc and map W
cmap = crameri('batlow',49);
cmap(50,:) = [1 1 1];
cmap = flipud(cmap);

for j=2:4
    
    %%
    mod = sims{j};

    %Corr
    Scorr  = nan(ni,nj);
    Mcorr  = nan(ni,nj);
    Lcorr  = nan(ni,nj);
    Fcorr  = nan(ni,nj);
    Pcorr  = nan(ni,nj);
    Dcorr  = nan(ni,nj);
    Vcorr  = nan(ni,nj);
    Bcorr  = nan(ni,nj);

    Scorr(ID)  = 1-rS(:,j);
    Mcorr(ID)  = 1-rM(:,j);
    Lcorr(ID)  = 1-rL(:,j);
    Fcorr(ID)  = 1-rF(:,j);
    Pcorr(ID)  = 1-rP(:,j);
    Dcorr(ID)  = 1-rD(:,j);
    Vcorr(ID)  = 1-rA(:,j);
    Bcorr(ID)  = 1-rB(:,j);

    %significance
    Ssig  = ones(ni,nj);
    Msig  = ones(ni,nj);
    Lsig  = ones(ni,nj);
    Fsig  = ones(ni,nj);
    Psig  = ones(ni,nj);
    Dsig  = ones(ni,nj);
    Vsig  = ones(ni,nj);
    Bsig  = ones(ni,nj);

    Ssig(ID)  = pS(:,j);
    Msig(ID)  = pM(:,j);
    Lsig(ID)  = pL(:,j);
    Fsig(ID)  = pF(:,j);
    Psig(ID)  = pP(:,j);
    Dsig(ID)  = pD(:,j);
    Vsig(ID)  = pA(:,j);
    Bsig(ID)  = pB(:,j);

    Ssig   = (Ssig<0.05);
    Msig   = (Msig<0.05);
    Lsig   = (Lsig<0.05);
    Fsig   = (Fsig<0.05);
    Psig   = (Psig<0.05);
    Dsig   = (Dsig<0.05);
    Vsig   = (Vsig<0.05);
    Bsig   = (Bsig<0.05);

    %% Fix seam
    [lat_s,lon_s,Scorr2] = cyclic_map_seam_cesm(TLAT,TLONG,Scorr);
    [~,~,Mcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Mcorr);
    [~,~,Lcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Lcorr);
    [~,~,Fcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Fcorr);
    [~,~,Pcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Pcorr);
    [~,~,Dcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Dcorr);
    [~,~,Bcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Bcorr);
    [~,~,Vcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Vcorr);


    %% Figure by size & type
    f6 = figure('Units','inches','Position',[1 3 16 8]);

    %1,1 - forage
    subplot('Position',[0.01 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Fcorr2)
    colormap(cmap) %crameri('-batlow')
    caxis([0 2])
    stipplem(TLAT,TLONG,Fsig,'markersize',3,'density',150)
        %         text(-2.5,2.25,'a','FontWeight','bold','FontSize',14)
    %         text(0,2.2,'Historic corr2','HorizontalAlignment','center','FontWeight','bold')
    text(0,1.675,'Forage fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %2,1 - all verts
    subplot('Position',[0.16 0.375 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Vcorr2)
    colormap(cmap)
    caxis([0 2])
    stipplem(TLAT,TLONG,Vsig,'markersize',3,'density',150)
        text(0,1.675,'All fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,1 - small
    subplot('Position',[0.01 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr2)
    colormap(cmap)
    caxis([0 2])
    stipplem(TLAT,TLONG,Ssig,'markersize',3,'density',150)
        text(0,1.675,'Small fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %1,2 - lg pel
    subplot('Position',[0.32 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Pcorr2)
    colormap(cmap)
    caxis([0 2])
    stipplem(TLAT,TLONG,Psig,'markersize',3,'density',150)
        %         text(-2.5,2.25,'b','FontWeight','bold','FontSize',14)
    text(0,1.85,mod,'HorizontalAlignment','center','FontWeight','bold')
    text(0,1.675,'Large pelagics','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %2,2 -
    %         subplot('Position',[0.32 0.375 0.3 0.275])
    %         axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    %             'Grid','off','FLineWidth',1)
    %         surfm(lat_s,lon_s,pdiff_nc2)
    %         colormap(cmap)
    %         caxis([0 2])
    %         text(0,1.675,'CNRM','HorizontalAlignment','center')
    %         h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    %
    %3,2 - med
    subplot('Position',[0.32 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr2)
    colormap(cmap)
    caxis([0 2])
    stipplem(TLAT,TLONG,Msig,'markersize',3,'density',150)
        text(0,1.675,'Medium fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %1,3 - dems
    subplot('Position',[0.63 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Dcorr2)
    colormap(cmap)
    caxis([0 2])
    stipplem(TLAT,TLONG,Dsig,'markersize',3,'density',150)
        %         text(-2.5,2.25,'c','FontWeight','bold','FontSize',14)
    %         text(0,2.2,'% \Delta mesoz','HorizontalAlignment','center','FontWeight','bold')
    text(0,1.675,'Demersals','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %2,3 - benthos
    subplot('Position',[0.475 0.375 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Bcorr2)
    colormap(cmap)
    caxis([0 2])
    stipplem(TLAT,TLONG,Bsig,'markersize',3,'density',150)
        text(0,1.675,'Benthic inverts','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,3 - large
    subplot('Position',[0.63 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr2)
    colormap(cmap)
    caxis([0 2])
    stipplem(TLAT,TLONG,Lsig,'markersize',3,'density',150)
    text(0,1.675,'Large fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    colorbar('Position',[0.325 0.035 0.3 0.02],'orientation','horizontal');

    print('-dpng',[ppath 'Map_W_FOSI_full_',mod,'types.png'])

end

%% Calc and map W, NO STIPPLING
for j=2:4
    
    %%
    mod = sims{j};

    %Corr
    Scorr  = nan(ni,nj);
    Mcorr  = nan(ni,nj);
    Lcorr  = nan(ni,nj);
    Fcorr  = nan(ni,nj);
    Pcorr  = nan(ni,nj);
    Dcorr  = nan(ni,nj);
    Vcorr  = nan(ni,nj);
    Bcorr  = nan(ni,nj);

    Scorr(ID)  = 1-rS(:,j);
    Mcorr(ID)  = 1-rM(:,j);
    Lcorr(ID)  = 1-rL(:,j);
    Fcorr(ID)  = 1-rF(:,j);
    Pcorr(ID)  = 1-rP(:,j);
    Dcorr(ID)  = 1-rD(:,j);
    Vcorr(ID)  = 1-rA(:,j);
    Bcorr(ID)  = 1-rB(:,j);

    %significance
    Ssig  = ones(ni,nj);
    Msig  = ones(ni,nj);
    Lsig  = ones(ni,nj);
    Fsig  = ones(ni,nj);
    Psig  = ones(ni,nj);
    Dsig  = ones(ni,nj);
    Vsig  = ones(ni,nj);
    Bsig  = ones(ni,nj);

    Ssig(ID)  = pS(:,j);
    Msig(ID)  = pM(:,j);
    Lsig(ID)  = pL(:,j);
    Fsig(ID)  = pF(:,j);
    Psig(ID)  = pP(:,j);
    Dsig(ID)  = pD(:,j);
    Vsig(ID)  = pA(:,j);
    Bsig(ID)  = pB(:,j);

    Ssig   = (Ssig<0.05);
    Msig   = (Msig<0.05);
    Lsig   = (Lsig<0.05);
    Fsig   = (Fsig<0.05);
    Psig   = (Psig<0.05);
    Dsig   = (Dsig<0.05);
    Vsig   = (Vsig<0.05);
    Bsig   = (Bsig<0.05);

    %% Fix seam
    [lat_s,lon_s,Scorr2] = cyclic_map_seam_cesm(TLAT,TLONG,Scorr);
    [~,~,Mcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Mcorr);
    [~,~,Lcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Lcorr);
    [~,~,Fcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Fcorr);
    [~,~,Pcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Pcorr);
    [~,~,Dcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Dcorr);
    [~,~,Bcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Bcorr);
    [~,~,Vcorr2]   = cyclic_map_seam_cesm(TLAT,TLONG,Vcorr);


    %% Figure by size & type
    f6 = figure('Units','inches','Position',[1 3 16 8]);

    %1,1 - forage
    subplot('Position',[0.01 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Fcorr2)
    colormap(cmap) %crameri('-batlow')
    caxis([0 2])
    text(0,1.675,'Forage fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %2,1 - all verts
    subplot('Position',[0.16 0.375 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Vcorr2)
    colormap(cmap)
    caxis([0 2])
    text(0,1.675,'All fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,1 - small
    subplot('Position',[0.01 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Scorr2)
    colormap(cmap)
    caxis([0 2])
    text(0,1.675,'Small fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %1,2 - lg pel
    subplot('Position',[0.32 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Pcorr2)
    colormap(cmap)
    caxis([0 2])
    text(0,1.85,mod,'HorizontalAlignment','center','FontWeight','bold')
    text(0,1.675,'Large pelagics','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,2 - med
    subplot('Position',[0.32 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Mcorr2)
    colormap(cmap)
    caxis([0 2])
    text(0,1.675,'Medium fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %1,3 - dems
    subplot('Position',[0.63 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Dcorr2)
    colormap(cmap)
    caxis([0 2])
    text(0,1.675,'Demersals','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %2,3 - benthos
    subplot('Position',[0.475 0.375 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Bcorr2)
    colormap(cmap)
    caxis([0 2])
    text(0,1.675,'Benthic inverts','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,3 - large
    subplot('Position',[0.63 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(lat_s,lon_s,Lcorr2)
    colormap(cmap)
    caxis([0 2])
    text(0,1.675,'Large fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    colorbar('Position',[0.325 0.035 0.3 0.02],'orientation','horizontal');

    print('-dpng',[ppath 'Map_W_FOSI_full_',mod,'types_NOstipple.png'])

end

