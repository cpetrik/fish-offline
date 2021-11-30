% Group output of FEISTY
% CESM DPLE members by start year
% Time series plots

clear all
close all

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

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
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%pick year
StartYr = 1954;
%loop over members
submem = [1,10:20];

tF = nan*ones(length(submem),120);
tP = tF;
tD = tF;
tB = tF;
sF = nan*ones(ni,nj,length(submem));
sP = sF;
sD = sF;
sS = sF;
sM = sF;
sL = sF;

for mem=1:length(submem) %will loop over
    Member = submem(mem);
    exper = ['v14_Y' num2str(StartYr) '_M' num2str(Member) '_All_fish03_' ];
    
    load([fpath 'Plot_Means_DPLE_' exper cfile '.mat'],'F','P','D','B',...
        'AllF','AllP','AllD','AllS','AllM','AllL');
    
    tF(mem,:) = F;
    tP(mem,:) = P;
    tD(mem,:) = D;
    tB(mem,:) = B;
    sF(:,:,mem) = AllF;
    sP(:,:,mem) = AllP;
    sD(:,:,mem) = AllD;
    sS(:,:,mem) = AllS;
    sM(:,:,mem) = AllM;
    sL(:,:,mem) = AllL;
    
end
exper2 = ['v14_Y' num2str(StartYr) '_allM_All_fish03_' ];
save([fpath 'Plot_Means_DPLE_' exper2 cfile '.mat'],'tF','tP','tD','tB',...
        'sF','sP','sD','sS','sM','sL')
    
%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black
    
set(groot,'defaultAxesColorOrder',cm10);

cmBP50=cbrewer('seq','BuPu',50,'PCHIP');

%% Plots in time
t = 1:120; %time;
y = t/12;

% All types
figure(1)
subplot(2,2,1)
plot(y,log10(tF),'r','Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('F')

subplot(2,2,2)
plot(y,log10(tP),'b','Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('P')

subplot(2,2,3)
plot(y,log10(tD),'k','Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('D')

subplot(2,2,4)
plot(y,log10(tB),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('B')
stamp(exper)
%     print('-dpng',[ppath 'DPLE_',exper,'all_types.png'])