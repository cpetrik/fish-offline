% Calc LME biomass of FEISTY
% CESM DPLE

clear all
close all

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/DPLE/';
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%pick year
StartYr = 2015;
%loop over members
submem = 1:40;

load([fpath 'LME_DPLE_fished_ensemble',cfile '.mat']);

%% Plot ts of Types together
    F = squeeze(lme_msfb(3,:,:)+lme_mmfb(3,:,:));
    P = squeeze(lme_mspb(3,:,:)+lme_mmpb(3,:,:)+lme_mlpb(3,:,:));
    D = squeeze(lme_msdb(3,:,:)+lme_mmdb(3,:,:)+lme_mldb(3,:,:));
    B = squeeze(lme_mbb(3,:,:));
    
    %%
    y=2015:2024;
    
    figure(2)
    subplot(2,2,1)
    plot(y,(F),'r','Linewidth',1); hold on;
    plot(y,nanmean(F,2),'color',[0.75 0 0.25],'Linewidth',2); hold on;
    xlim([y(1) y(end)])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    title('Forage fish')
    
    subplot(2,2,2)
    plot(y,(P),'b','Linewidth',1); hold on;
    plot(y,nanmean(P,2),'color',[0 0.25 0.5],'Linewidth',2); hold on;
    xlim([y(1) y(end)])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    title('Large Pelagic fish')
    
    subplot(2,2,3)
    plot(y,(D),'color',[0.25 0.25 0.25],'Linewidth',1); hold on;
    plot(y,nanmean(D,2),'k','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    title('Demersal fish')
    
    subplot(2,2,4)
    plot(y,(B),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
    plot(y,nanmean(B,2),'color',[0.35 0.35 0.35],'Linewidth',2); hold on;
    xlim([y(1) y(end)])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    title('Benthic Invertebrates')
    stamp(exper)
    print('-dpng',[ppath 'DPLE_v14_Y2015_All_fish03_CCLME_all_types.png'])
