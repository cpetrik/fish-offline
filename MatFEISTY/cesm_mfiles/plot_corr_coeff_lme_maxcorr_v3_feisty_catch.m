% Plot corr coeffs of driver-catch corrs
% Lag with max R2
% For all 63 LMEs
% Only positive corrs with prey, fish
% Use calc corr of catch with forcing
% v3 = consistent corrs only calc once

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',cfile,'/corrs/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

%% All corrs
CFtab = nan*ones(63,6,5);
PFtab = nan*ones(63,6,5);
CPtab = CFtab;
PPtab = CFtab;
CDtab = CFtab;
PDtab = CFtab;
CAtab = CFtab;
PAtab = CFtab;

%% sat & driver
load([spath 'LMEs_corr_catch_satyrs_driver_lags.mat'])
stex = tanom(1:4);
%dim: LME x driver x lags

%drivers & sat
CAtab(:,1:4,:) = AtabC(:,1:4,:);
CFtab(:,1:4,:) = FtabC(:,1:4,:);
CPtab(:,1:4,:) = PtabC(:,1:4,:);
CDtab(:,1:4,:) = DtabC(:,1:4,:);

PAtab(:,1:4,:) = AtabP(:,1:4,:);
PFtab(:,1:4,:) = FtabP(:,1:4,:);
PPtab(:,1:4,:) = PtabP(:,1:4,:);
PDtab(:,1:4,:) = DtabP(:,1:4,:);

clear AtabC AtabP FtabC FtabP PtabC PtabP DtabC DtabP tanom

%%
load([spath 'LMEs_corr_catch_satyrs_feisty_lags.mat'])
ftex = tanom;
%dim: LME x driver x lags

%fish - const effort
CAtab(:,5:6,1:4) = AtabC(:,1:2,:);
CFtab(:,5:6,1:4) = FtabC(:,1:2,:);
CPtab(:,5:6,1:4) = PtabC(:,1:2,:);
CDtab(:,5:6,1:4) = DtabC(:,1:2,:);

PAtab(:,5:6,1:4) = AtabP(:,1:2,:);
PFtab(:,5:6,1:4) = FtabP(:,1:2,:);
PPtab(:,5:6,1:4) = PtabP(:,1:2,:);
PDtab(:,5:6,1:4) = DtabP(:,1:2,:);

clear AtabC AtabP FtabC FtabP PtabC PtabP DtabC DtabP tanom

%%
ctex = [stex,ftex];

%All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

%%  -------- Find lag with maxcorr ---------------------

yr = 0:4;

Fcor = nan*ones(63,6);
Fsig = nan*ones(63,6);
Pcor = Fcor;
Psig = Fcor;
Dcor = Fcor;
Dsig = Fcor;
Acor = Fcor;
Asig = Fcor;

for L = 1:length(lid)
    for d=1:length(ctex)
        AtabC = squeeze(CAtab(L,d,:));
        AtabP = squeeze(PAtab(L,d,:));
        maxC = max(abs(AtabC));
        pid = find(abs(AtabC)==maxC);
        Acor(L,d) = AtabC(pid);
        Asig(L,d) = AtabP(pid);
        clear maxC pid

        if(L~=61)
            if(L~=55)
                FtabC = squeeze(CFtab(L,d,:));
                FtabP = squeeze(PFtab(L,d,:));
                maxC = max(abs(FtabC));
                pid = find(abs(FtabC)==maxC);
                Fcor(L,d) = FtabC(pid);
                Fsig(L,d) = FtabP(pid);
                clear maxC pid
            end
        end

        if(L~=61)
            if(L~=63)
                PtabC = squeeze(CPtab(L,d,:));
                PtabP = squeeze(PPtab(L,d,:));
                maxC = max(abs(PtabC));
                pid = find(abs(PtabC)==maxC);
                Pcor(L,d) = PtabC(pid);
                Psig(L,d) = PtabP(pid);
                clear maxC pid
            end
        end

        DtabC = squeeze(CDtab(L,d,:));
        DtabP = squeeze(PDtab(L,d,:));
        maxC = max(abs(DtabC));
        pid = find(abs(DtabC)==maxC);
        Dcor(L,d) = DtabC(pid);
        Dsig(L,d) = DtabP(pid);
        clear maxC pid
    end
end

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
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%% Correlations w/ drivers map
% Correlation value map
cmR = cbrewer('seq','Reds',10,'PCHIP');
cmRB = cbrewer('div','RdBu',21,'PCHIP');
cmRB = flipud(cmRB);

tanom = ctex;

%Now these are fn type biom & prod, correlated with Catch fn types,

figure(1)
for j=1:2 %pos & neg w/temp
    clf
    Fcorr  = nan(ni,nj);
    Pcorr  = nan(ni,nj);
    Dcorr  = nan(ni,nj);
    Acorr  = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (Fsig(i,j) <= 0.05)
            Fcorr(id) = Fcor(i,j);
        end
        if (Psig(i,j) <= 0.05)
            Pcorr(id) = Pcor(i,j);
        end
        if (Dsig(i,j) <= 0.05)
            Dcorr(id) = Dcor(i,j);
        end
        if (Asig(i,j) <= 0.05)
            Acorr(id) = Acor(i,j);
        end
    end

    subplot('Position',[0 0.53 0.5 0.5])
    %F
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Fcorr)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 1]);
    set(gcf,'renderer','painters')
    title(['corr Catch Forage with ' tanom{j}])

    %P
    subplot('Position',[0.5 0.53 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Pcorr)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 1]);
    set(gcf,'renderer','painters')
    title(['corr Catch LgPel with ' tanom{j}])

    %D
    subplot('Position',[0.0 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Dcorr)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 1]);
    set(gcf,'renderer','painters')
    title(['corr Catch Dem with ' tanom{j}])

    %All
    subplot('Position',[0.5 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Acorr)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 1]);
    colorbar('Position',[0.2 0.5 0.6 0.05],'orientation','horizontal','AxisLocation','in')
    set(gcf,'renderer','painters')
    title(['corr Catch All with ' tanom{j}])
    %stamp('')

    print('-dpng',[ppath 'Map_LMEs_corr_coeff_v3_',tanom{j},'_catch_fntypes.png'])
end

%% Just pos prey, fish corrs
for j=3:6
    Fcorr  = nan(ni,nj);
    Pcorr  = nan(ni,nj);
    Dcorr  = nan(ni,nj);
    Acorr  = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (Fsig(i,j) <= 0.05)
            if (Fcor(i,j) >= 0)
                Fcorr(id) = Fcor(i,j);
            end
        end
        if (Psig(i,j) <= 0.05)
            if (Pcor(i,j) >= 0)
                Pcorr(id) = Pcor(i,j);
            end
        end
        if (Dsig(i,j) <= 0.05)
            if (Dcor(i,j) >= 0)
                Dcorr(id) = Dcor(i,j);
            end
        end
        if (Asig(i,j) <= 0.05)
            if (Acor(i,j) >= 0)
                Acorr(id) = Acor(i,j);
            end
        end
    end

    figure(2)
    clf
    subplot('Position',[0 0.53 0.5 0.5])
    %F
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Fcorr)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colormap(cmR)
    clim([0 1]);
    set(gcf,'renderer','painters')
    title(['corr Catch Forage with ' tanom{j}])

    %P
    subplot('Position',[0.5 0.53 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Pcorr)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colormap(cmR)
    clim([0 1]);
    set(gcf,'renderer','painters')
    title(['corr Catch LgPel with ' tanom{j}])

    %D
    subplot('Position',[0.0 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Dcorr)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colormap(cmR)
    clim([0 1]);
    set(gcf,'renderer','painters')
    title(['corr Catch Dem with ' tanom{j}])

    %All
    subplot('Position',[0.5 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Acorr)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colormap(cmR)
    clim([0 1]);
    colorbar('Position',[0.2 0.5 0.6 0.05],'orientation','horizontal','AxisLocation','in')
    set(gcf,'renderer','painters')
    title(['corr Catch All with ' tanom{j}])
    %stamp('')
    print('-dpng',[ppath 'Map_LMEs_pos_corr_coeff_v3_',tanom{j},'_catch_fntypes.png'])

    % %All
    % figure(5)
    % clf
    % axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    %     'Grid','off','FLineWidth',1)
    % surfm(TLAT,TLONG,Acorr)
    % h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    % colormap(cmR)
    % clim([0 1]);
    % colorbar
    % set(gcf,'renderer','painters')
    % title(['corr Catch All with ' tanom{j}])
    % %stamp('')
    % print('-dpng',[ppath 'Map_LMEs_pos_corr_coeff_v3_',tanom{j},'_catch_all.png'])

end


%% Just pos fish corrs
for j=5:6
    Fcorr  = nan(ni,nj);
    Pcorr  = nan(ni,nj);
    Dcorr  = nan(ni,nj);
    Acorr  = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (Fsig(i,j) <= 0.05)
            Fcorr(id) = Fcor(i,j);
        end
        if (Psig(i,j) <= 0.05)
            Pcorr(id) = Pcor(i,j);
        end
        if (Dsig(i,j) <= 0.05)
            Dcorr(id) = Dcor(i,j);
        end
        if (Asig(i,j) <= 0.05)
            Acorr(id) = Acor(i,j);
        end
    end

    figure(3)
    clf
    subplot('Position',[0 0.53 0.5 0.5])
    %F
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Fcorr.^2)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colormap(cmR)
    clim([0 1]);
    set(gcf,'renderer','painters')
    title(['R^2 of Catch Forage corr with ' tanom{j}])

    %P
    subplot('Position',[0.5 0.53 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Pcorr.^2)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colormap(cmR)
    clim([0 1]);
    set(gcf,'renderer','painters')
    title(['R^2 of Catch LgPel corr with ' tanom{j}])

    %D
    subplot('Position',[0.0 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Dcorr.^2)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colormap(cmR)
    clim([0 1]);
    set(gcf,'renderer','painters')
    title(['R^2 of Catch Dem corr with ' tanom{j}])

    %All
    subplot('Position',[0.5 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Acorr.^2)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colormap(cmR)
    clim([0 1]);
    colorbar('Position',[0.2 0.5 0.6 0.05],'orientation','horizontal','AxisLocation','in')
    set(gcf,'renderer','painters')
    title(['R^2 of Catch All corr with ' tanom{j}])
    %stamp('')
    print('-dpng',[ppath 'Map_LMEs_pos_corr_coeff_v3_',tanom{j},'_catch_fntypes_R2.png'])

end

%% Just all fish Catch with ind drivers
TPcorr  = nan(ni,nj);
TBcorr  = nan(ni,nj);
Dcorr  = nan(ni,nj);
Zcorr  = nan(ni,nj);

for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);

    if (Asig(i,1) <= 0.05)
        TPcorr(id) = Acor(i,1);
    end
    if (Asig(i,2) <= 0.05)
        TBcorr(id) = Acor(i,2);
    end
    if (Asig(i,3) <= 0.05)
        Dcorr(id) = Acor(i,3);
    end
    if (Asig(i,4) <= 0.05)
        Zcorr(id) = Acor(i,4);
    end
end

figure(7)
subplot('Position',[0 0.53 0.5 0.5])
%TP
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,TPcorr)
colormap(cmRB)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
set(gcf,'renderer','painters')
title('TP')

%ZmL
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,Zcorr)
colormap(cmRB)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
set(gcf,'renderer','painters')
title('ZmLoss')

%TB
subplot('Position',[0.0 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,TBcorr)
colormap(cmRB)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
set(gcf,'renderer','painters')
title('TB')

%Det
subplot('Position',[0.5 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,Dcorr)
colormap(cmRB)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Det')
print('-dpng',[ppath 'Map_LMEs_corr_coeff_v3_drivers_posfood_catch_allfish.png'])

%% Just all fish Catch with ind drivers & biom, prod
TPcorr  = nan(ni,nj);
TBcorr  = nan(ni,nj);
Dcorr  = nan(ni,nj);
Zcorr  = nan(ni,nj);
Bcorr  = nan(ni,nj);
Pcorr  = nan(ni,nj);

for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);

    if (Asig(i,1) <= 0.05)
        TPcorr(id) = Acor(i,1);
    end
    if (Asig(i,2) <= 0.05)
        TBcorr(id) = Acor(i,2);
    end
    if (Asig(i,3) <= 0.05)
        Dcorr(id) = Acor(i,3);
    end
    if (Asig(i,4) <= 0.05)
        Zcorr(id) = Acor(i,4);
    end
    if (Asig(i,5) <= 0.05)
        Bcorr(id) = Acor(i,5);
    end
    if (Asig(i,6) <= 0.05)
        Pcorr(id) = Acor(i,6);
    end
end

%% Only pos corrs for fish
Zcorr2  = Zcorr;
Dcorr2  = Dcorr;
Bcorr2  = Bcorr;
Pcorr2  = Pcorr;
Zcorr2(Zcorr2<0) = nan;
Dcorr2(Dcorr2<0) = nan;
Bcorr2(Bcorr2<0) = nan;
Pcorr2(Pcorr2<0) = nan;

f9 = figure('Units','inches','Position',[1 3 7.5 4.5]);

subplot('Position',[0.01 0.55 0.32 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,TPcorr)
colormap(cmRB)
clim([-1 1]);
title('TP')

subplot('Position',[0.335 0.55 0.32 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Zcorr2)
colormap(cmRB)
clim([-1 1]);
title('ZmLoss')

subplot('Position',[0.66 0.55 0.32 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Bcorr2)
colormap(cmRB)
clim([-1 1]);
title('Biomass')

subplot('Position',[0.01 0.01 0.32 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,TBcorr)
colormap(cmRB)
clim([-1 1]);
title('TB')

subplot('Position',[0.335 0.01 0.32 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Dcorr2)
colormap(cmRB)
clim([-1 1]);
title('Det')

subplot('Position',[0.66 0.01 0.32 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Pcorr2)
colormap(cmRB)
clim([-1 1]);
title('Production')
colorbar('Position',[0.25 0.5 0.5 0.03],'Orientation','horizontal','AxisLocation','in')
print('-dpng',[ppath 'Map_LMEs_pos_corr_coeff_v3_drivers_catch_allfish.png'])

