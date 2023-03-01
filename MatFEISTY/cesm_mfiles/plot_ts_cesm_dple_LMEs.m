% Visualize output of FEISTY
% CESM DPLE
% Time series plots and maps

clear
close all

%% Map data
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/DPLE/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/DPLE/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%pick year
StartYr = 1954;

exper2 = ['v15_Y' num2str(StartYr) '_All_fish03' ];

load([fpath 'Time_Means_DPLE_LME_' exper2 '.mat'])

%%
Ftmean = lme_mSF + lme_mMF;
Ptmean = lme_mSP + lme_mMP + lme_mLP;
Dtmean = lme_mSD + lme_mMD + lme_mLD;

%% Plots in time
t = 1:122; %time;
y = (t/12) + (10/12) + StartYr;

% All types
figure(1)
subplot(2,2,1)
plot(y,(Ftmean),'r','Linewidth',0.5); hold on;
plot(y,(mean(Ftmean)),'color',[0.75 0 0],'Linewidth',3); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title('F')

subplot(2,2,2)
plot(y,(Ptmean),'b','Linewidth',0.5); hold on;
plot(y,(mean(Ptmean)),'color',[0 0 0.75],'Linewidth',3); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title('P')

subplot(2,2,3)
plot(y,(Dtmean),'color',[0 0.75 0.5],'Linewidth',0.5); hold on;
plot(y,(mean(Dtmean)),'color',[0 0.5 0.25],'Linewidth',3); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title('D')

subplot(2,2,4)
plot(y,squeeze(lme_mB(3,:,:)),'color',[0.5 0.5 0.5],'Linewidth',0.5); hold on;
plot(y,mean(squeeze(lme_mB(3,:,:))),'k','Linewidth',3); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title('B')
stamp(exper2)
print('-dpng',[ppath 'DPLE_',exper2,'all_types.png'])


