% Plot coefs of driver-cpue SLR
% Lag with max R2
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',cfile,'/corrs/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

load([spath,'LME_linreg_maxR2lag_driver_cpue.mat'])

%% colorblind friendly
load('paul_tol_cmaps.mat')

% mcol = [238/255 102/255 119/255;... %red - TP
%     0/255 68/255 136/255;...    %blue - TB
%     34/255 136/255 51/255;...   %green - Det
%     51/255 187/255 238/255;...  %cyan - ZL
%     ];

%colorblind friendly - subselect & re-order drainbow
%ctex = {'TP','TB','Det','ZmLoss','Biom','Prod','SST','chl'};
% orange, dk blue, grey, lt blue, dk purp, lt purp, red, green
mcol(1,:) = drainbow(12,:)/255; % orange
mcol(2,:) = drainbow(4,:)/255; %dk blue
mcol(3,:) = drainbow(15,:)/255; %grey
mcol(4,:) = drainbow(6,:)/255; %lt blue
mcol(5,:) = drainbow(3,:)/255; %dk purp
mcol(6,:) = drainbow(1,:)/255; %lt purp
mcol(7,:) = drainbow(14,:)/255; %red
mcol(8,:) = drainbow(7,:)/255; %green

%% Map
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

lid = LAtab(:,7);

%These are LRs with CPUE fn types by fn type

figure(1)
for j=1:6
    clf
    Fcorr  = nan(ni,nj);
    Pcorr  = nan(ni,nj);
    Dcorr  = nan(ni,nj);
    Acorr  = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (LFsig(i,j) <= 0.05)
            Fcorr(id) = LFtab(i,j);
        end
        if (LPsig(i,j) <= 0.05)
            Pcorr(id) = LPtab(i,j);
        end
        if (LDsig(i,j) <= 0.05)
            Dcorr(id) = LDtab(i,j);
        end
        if (LAsig(i,j) <= 0.05)
            Acorr(id) = LAtab(i,j);
        end
    end

    subplot('Position',[0 0.53 0.5 0.5])
    %F
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Fcorr)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 1]);
    set(gcf,'renderer','painters')
    title(['LR CPUE Forage with ' tanom{j}])

    %P
    subplot('Position',[0.5 0.53 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Pcorr)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 1]);
    set(gcf,'renderer','painters')
    title(['LR CPUE LgPel with ' tanom{j}])

    %D
    subplot('Position',[0.0 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Dcorr)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 1]);
    set(gcf,'renderer','painters')
    title(['LR CPUE Dem with ' tanom{j}])

    %All
    subplot('Position',[0.5 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Acorr)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 1]);
    colorbar('Position',[0.2 0.5 0.6 0.05],'orientation','horizontal','AxisLocation','in')
    set(gcf,'renderer','painters')
    title(['LR CPUE All with ' tanom{j}])
    %stamp('')

    print('-dpng',[ppath 'Map_LMEs_SLR_coeff_',tanom{j},'_cpue_fntypes.png'])
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

        if (LFsig(i,j) <= 0.05)
            if (LFtab(i,j) >= 0)
                Fcorr(id) = LFtab(i,j);
            end
        end
        if (LPsig(i,j) <= 0.05)
            if (LPtab(i,j) >= 0)
                Pcorr(id) = LPtab(i,j);
            end
        end
        if (LDsig(i,j) <= 0.05)
            if (LDtab(i,j) >= 0)
                Dcorr(id) = LDtab(i,j);
            end
        end
        if (LAsig(i,j) <= 0.05)
            if (LAtab(i,j) >= 0)
                Acorr(id) = LAtab(i,j);
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
        caxis([0 1]);
        set(gcf,'renderer','painters')
        title(['LR CPUE Forage with ' tanom{j}])
    
        %P
        subplot('Position',[0.5 0.53 0.5 0.5])
        axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(TLAT,TLONG,Pcorr)
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        colormap(cmR)
        caxis([0 1]);
        set(gcf,'renderer','painters')
        title(['LR CPUE LgPel with ' tanom{j}])
    
        %D
        subplot('Position',[0.0 0.0 0.5 0.5])
        axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(TLAT,TLONG,Dcorr)
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        colormap(cmR)
        caxis([0 1]);
        colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
        set(gcf,'renderer','painters')
        title(['LR CPUE Dem with ' tanom{j}])
    
        %All
        subplot('Position',[0.5 0.0 0.5 0.5])
        axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
            'Grid','off','FLineWidth',1)
        surfm(TLAT,TLONG,Acorr)
        h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
        colormap(cmR)
        caxis([0 1]);
        colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
        set(gcf,'renderer','painters')
        title(['LR CPUE All with ' tanom{j}])
        %stamp('')
        print('-dpng',[ppath 'Map_LMEs_pos_LR_coeff_',tanom{j},'_cpue_fntypes.png'])

    %All
    figure(5)
    clf
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Acorr)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colormap(cmR)
    caxis([0 1]);
    colorbar
    set(gcf,'renderer','painters')
    title(['LR CPUE All with ' tanom{j}])
    %stamp('')
    print('-dpng',[ppath 'Map_LMEs_pos_LR_coeff_',tanom{j},'_cpue_all.png'])

end

%% pos & neg
for j=5:6
    Acorr  = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (LAsig(i,j) <= 0.05)
            Acorr(id) = LAtab(i,j);
        end
    end


    %All
    figure(6)
    clf
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Acorr)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colormap(cmRB)
    caxis([-1 1]);
    colorbar
    set(gcf,'renderer','painters')
    title(['LR CPUE All with ' tanom{j}])
    %stamp('')
    print('-dpng',[ppath 'Map_LMEs_LR_coeff_',tanom{j},'_cpue_all.png'])

end

%% Just all fish CPUE with ind drivers
TPcorr  = nan(ni,nj);
TBcorr  = nan(ni,nj);
Dcorr  = nan(ni,nj);
Zcorr  = nan(ni,nj);

for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);

    if (LAsig(i,1) <= 0.05)
        TPcorr(id) = LAtab(i,1);
    end
    if (LAsig(i,2) <= 0.05)
        TBcorr(id) = LAtab(i,2);
    end
    if (LAsig(i,3) <= 0.05)
        Dcorr(id) = LAtab(i,3);
    end
    if (LAsig(i,4) <= 0.05)
        Zcorr(id) = LAtab(i,4);
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
caxis([-1 1]);
set(gcf,'renderer','painters')
title('TP')

%ZmL
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,Zcorr)
colormap(cmRB)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('ZmLoss')

%TB
subplot('Position',[0.0 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,TBcorr)
colormap(cmRB)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('TB')

%Det
subplot('Position',[0.5 0.0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,Dcorr)
colormap(cmRB)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('Det')
print('-dpng',[ppath 'Map_LMEs_LR_coeff_drivers_cpue_allfish.png'])

%% Just all fish CPUE with ind drivers & biom, prod
TPcorr  = nan(ni,nj);
TBcorr  = nan(ni,nj);
Dcorr  = nan(ni,nj);
Zcorr  = nan(ni,nj);
Bcorr  = nan(ni,nj);
Pcorr  = nan(ni,nj);

for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);

    if (LAsig(i,1) <= 0.05)
        TPcorr(id) = LAtab(i,1);
    end
    if (LAsig(i,2) <= 0.05)
        TBcorr(id) = LAtab(i,2);
    end
    if (LAsig(i,3) <= 0.05)
        Dcorr(id) = LAtab(i,3);
    end
    if (LAsig(i,4) <= 0.05)
        Zcorr(id) = LAtab(i,4);
    end
    if (LAsig(i,5) <= 0.05)
        Bcorr(id) = LAtab(i,5);
    end
    if (LAsig(i,6) <= 0.05)
        Pcorr(id) = LAtab(i,6);
    end
end

%% All corrs
f8 = figure('Units','inches','Position',[1 3 7.5 4.5]);

subplot('Position',[0.01 0.55 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,TPcorr)
colormap(cmRB)
caxis([-1 1]);
title('TP')

subplot('Position',[0.335 0.55 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Zcorr)
colormap(cmRB)
caxis([-1 1]);
title('ZmLoss')

subplot('Position',[0.66 0.55 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Bcorr)
colormap(cmRB)
caxis([-1 1]);
title('Biomass')

subplot('Position',[0.01 0.01 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,TBcorr)
colormap(cmRB)
caxis([-1 1]);
title('TB')

subplot('Position',[0.335 0.01 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Dcorr)
colormap(cmRB)
caxis([-1 1]);
title('Det')

subplot('Position',[0.66 0.01 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Pcorr)
colormap(cmRB)
caxis([-1 1]);
title('Production')
colorbar('Position',[0.25 0.5 0.5 0.03],'Orientation','horizontal','AxisLocation','in')
print('-dpng',[ppath 'Map_LMEs_LR_coeff_drivers_fish_cpue_allfish.png'])

%% Only pos corrs for fish
Bcorr2  = Bcorr;
Pcorr2  = Pcorr;
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
caxis([-1 1]);
title('TP')

subplot('Position',[0.335 0.55 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Zcorr)
colormap(cmRB)
caxis([-1 1]);
title('ZmLoss')

subplot('Position',[0.66 0.55 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Bcorr2)
colormap(cmRB)
caxis([-1 1]);
title('Biomass')

subplot('Position',[0.01 0.01 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,TBcorr)
colormap(cmRB)
caxis([-1 1]);
title('TB')

subplot('Position',[0.335 0.01 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Dcorr)
colormap(cmRB)
caxis([-1 1]);
title('Det')

subplot('Position',[0.66 0.01 0.32 0.4]) 
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Pcorr2)
colormap(cmRB)
caxis([-1 1]);
title('Production')
colorbar('Position',[0.25 0.5 0.5 0.03],'Orientation','horizontal','AxisLocation','in')
print('-dpng',[ppath 'Map_LMEs_LR_coeff_drivers_posfish_cpue_allfish.png'])

