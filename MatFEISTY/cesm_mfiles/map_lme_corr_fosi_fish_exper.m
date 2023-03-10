% Corr LME biomass of FEISTY experiments
% with full (and each other?)
% CESM FOSI

clear
close all

%% Map data
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

tlme = double(lme_mask);
% tlme = double(lme_mask(ID));
% tlme(tlme<-1e9) = nan;

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

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
ppath = [pp cfile '/corrs/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

load([fpath 'LME_corrs_FOSI_fished_expers.mat'])

simtex = {'climatol';'varTemp';'varFood'};

%% loop over expers
for j=1:3

    mod = simtex{j};

    %% on grid
    Scorr  = nan(ni,nj);
    Mcorr  = nan(ni,nj);
    Lcorr  = nan(ni,nj);
    Fcorr  = nan(ni,nj);
    Pcorr  = nan(ni,nj);
    Dcorr  = nan(ni,nj);
    Vcorr  = nan(ni,nj);
    Bcorr  = nan(ni,nj);

    %significance
    Ssig  = ones(ni,nj);
    Msig  = ones(ni,nj);
    Lsig  = ones(ni,nj);
    Fsig  = ones(ni,nj);
    Psig  = ones(ni,nj);
    Dsig  = ones(ni,nj);
    Vsig  = ones(ni,nj);
    Bsig  = ones(ni,nj);

    for L=1:66
        lid = find(tlme==L);

        Scorr(lid) = rS(L,j);
        Mcorr(lid) = rM(L,j);
        Lcorr(lid) = rL(L,j);
        Fcorr(lid) = rF(L,j);
        Pcorr(lid) = rP(L,j);
        Dcorr(lid) = rD(L,j);
        Vcorr(lid) = rA(L,j);
        Bcorr(lid) = rB(L,j);
        
        Ssig(lid) = pS(L,j);
        Msig(lid) = pM(L,j);
        Lsig(lid) = pL(L,j);
        Fsig(lid) = pF(L,j);
        Psig(lid) = pP(L,j);
        Dsig(lid) = pD(L,j);
        Vsig(lid) = pA(L,j);
        Bsig(lid) = pB(L,j);

    end

    %% Figure by size & type
    f5 = figure('Units','inches','Position',[1 3 16 8]);

    %1,1 - forage
    subplot('Position',[0.01 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Fcorr)
    cmocean('balance')
    caxis([-1 1])
    %stipplem(TLAT,TLONG,Fsig,'markersize',3,'density',150)
        %         text(-2.5,2.25,'a','FontWeight','bold','FontSize',14)
    %         text(0,2.2,'Historic corr2','HorizontalAlignment','center','FontWeight','bold')
    text(0,1.675,'Forage fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %2,1 - all verts
    subplot('Position',[0.16 0.375 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Vcorr)
    cmocean('balance')
    caxis([-1 1])
    %stipplem(TLAT,TLONG,Vsig,'markersize',3,'density',150)
        text(0,1.675,'All fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,1 - small
    subplot('Position',[0.01 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Scorr)
    cmocean('balance')
    caxis([-1 1])
    %stipplem(TLAT,TLONG,Ssig,'markersize',3,'density',150)
        text(0,1.675,'Small fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %1,2 - lg pel
    subplot('Position',[0.32 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Pcorr)
    cmocean('balance')
    caxis([-1 1])
    %stipplem(TLAT,TLONG,Psig,'markersize',3,'density',150)
        %         text(-2.5,2.25,'b','FontWeight','bold','FontSize',14)
    text(0,1.85,mod,'HorizontalAlignment','center','FontWeight','bold')
    text(0,1.675,'Large pelagics','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,2 - med
    subplot('Position',[0.32 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Mcorr)
    cmocean('balance')
    caxis([-1 1])
    %stipplem(TLAT,TLONG,Msig,'markersize',3,'density',150)
        text(0,1.675,'Medium fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %1,3 - dems
    subplot('Position',[0.63 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Dcorr)
    cmocean('balance')
    caxis([-1 1])
    %stipplem(TLAT,TLONG,Dsig,'markersize',3,'density',150)
        %         text(-2.5,2.25,'c','FontWeight','bold','FontSize',14)
    %         text(0,2.2,'% \Delta mesoz','HorizontalAlignment','center','FontWeight','bold')
    text(0,1.675,'Demersals','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %2,3 - benthos
    subplot('Position',[0.475 0.375 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Bcorr)
    cmocean('balance')
    caxis([-1 1])
    %stipplem(TLAT,TLONG,Bsig,'markersize',3,'density',150)
        text(0,1.675,'Benthic inverts','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,3 - large
    subplot('Position',[0.63 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Lcorr)
    cmocean('balance')
    caxis([-1 1])
    %stipplem(TLAT,TLONG,Lsig,'markersize',3,'density',150)
    text(0,1.675,'Large fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    colorbar('Position',[0.325 0.035 0.3 0.02],'orientation','horizontal');

    print('-dpng',[ppath 'Map_lme_corr_FOSI_full_',mod,'_types.png'])


    %% W by size & type
    f6 = figure('Units','inches','Position',[1 3 16 8]);

    %1,1 - forage
    subplot('Position',[0.01 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,1-Fcorr)
    crameri('-batlow')
    caxis([0 2])
    %stipplem(TLAT,TLONG,Fsig,'markersize',3,'density',150)
        %         text(-2.5,2.25,'a','FontWeight','bold','FontSize',14)
    %         text(0,2.2,'Historic corr2','HorizontalAlignment','center','FontWeight','bold')
    text(0,1.675,'Forage fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %2,1 - all verts
    subplot('Position',[0.16 0.375 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,1-Vcorr)
    crameri('-batlow')
    caxis([0 2])
    %stipplem(TLAT,TLONG,Vsig,'markersize',3,'density',150)
        text(0,1.675,'All fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,1 - small
    subplot('Position',[0.01 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,1-Scorr)
    crameri('-batlow')
    caxis([0 2])
    %stipplem(TLAT,TLONG,Ssig,'markersize',3,'density',150)
        text(0,1.675,'Small fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %1,2 - lg pel
    subplot('Position',[0.32 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,1-Pcorr)
    crameri('-batlow')
    caxis([0 2])
    %stipplem(TLAT,TLONG,Psig,'markersize',3,'density',150)
        %         text(-2.5,2.25,'b','FontWeight','bold','FontSize',14)
    text(0,1.85,mod,'HorizontalAlignment','center','FontWeight','bold')
    text(0,1.675,'Large pelagics','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,2 - med
    subplot('Position',[0.32 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,1-Mcorr)
    crameri('-batlow')
    caxis([0 2])
    %stipplem(TLAT,TLONG,Msig,'markersize',3,'density',150)
        text(0,1.675,'Medium fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %1,3 - dems
    subplot('Position',[0.63 0.675 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,1-Dcorr)
    crameri('-batlow')
    caxis([0 2])
    %stipplem(TLAT,TLONG,Dsig,'markersize',3,'density',150)
        %         text(-2.5,2.25,'c','FontWeight','bold','FontSize',14)
    %         text(0,2.2,'% \Delta mesoz','HorizontalAlignment','center','FontWeight','bold')
    text(0,1.675,'Demersals','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %2,3 - benthos
    subplot('Position',[0.475 0.375 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,1-Bcorr)
    crameri('-batlow')
    caxis([0 2])
    %stipplem(TLAT,TLONG,Bsig,'markersize',3,'density',150)
        text(0,1.675,'Benthic inverts','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    %3,3 - large
    subplot('Position',[0.63 0.075 0.3 0.275])
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,1-Lcorr)
    crameri('-batlow')
    caxis([0 2])
    %stipplem(TLAT,TLONG,Lsig,'markersize',3,'density',150)
    text(0,1.675,'Large fish','HorizontalAlignment','center')
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    colorbar('Position',[0.325 0.035 0.3 0.02],'orientation','horizontal');

    print('-dpng',[ppath 'Map_lme_W_FOSI_full_',mod,'_types.png'])

end

