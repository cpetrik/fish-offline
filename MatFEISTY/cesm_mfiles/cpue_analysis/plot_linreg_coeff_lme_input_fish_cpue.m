% Plot max SLR coeff of driver-cpue corrs
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

load([spath,'LMEs_SLR_cpue_driver_feisty_norec_maxR2s.mat']) %table of coeff of each driver at best lag

%%  ---------------------------------------------------------
cnam = {'coef','p','lag','idriver','driver'};
ctex = {'TP','TB','Det','ZmLoss','Biom','Prod'};

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

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
cmR = cbrewer('seq','Reds',9,'PCHIP');
cmRB = cbrewer('div','RdBu',21,'PCHIP');
cmRB = flipud(cmRB);

lid = LAtab(:,8);

figure(1)
for j=1:7
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
    title(['corr CPUE Forage with ' tanom{j}])

    %P
    subplot('Position',[0.5 0.53 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Pcorr)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 1]);
    set(gcf,'renderer','painters')
    title(['corr CPUE LgPel with ' tanom{j}])

    %D
    subplot('Position',[0.0 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Dcorr)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 1]);
    colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title(['corr CPUE Dem with ' tanom{j}])

    %All
    subplot('Position',[0.5 0.0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Acorr)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 1]);
    colorbar('Position',[0.2 0.485 0.6 0.05],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title(['corr CPUE All with ' tanom{j}])
    %stamp('')

    print('-dpng',[ppath 'Map_LMEs_corr_coeff_',tanom{j},'_cpue_fntypes.png'])
end

%% Just pos fish corrs
for j=5:7
    figure(2)
    clf
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

    subplot('Position',[0 0.53 0.5 0.5])
    %F
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Fcorr)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colormap(cmR)
    caxis([0 1]);
    set(gcf,'renderer','painters')
    title(['corr CPUE Forage with ' tanom{j}])

    %P
    subplot('Position',[0.5 0.53 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(TLAT,TLONG,Pcorr)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colormap(cmR)
    caxis([0 1]);
    set(gcf,'renderer','painters')
    title(['corr CPUE LgPel with ' tanom{j}])

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
    title(['corr CPUE Dem with ' tanom{j}])

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
    title(['corr CPUE All with ' tanom{j}])
    %stamp('')
    print('-dpng',[ppath 'Map_LMEs_pos_corr_coeff_',tanom{j},'_cpue_fntypes.png'])

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
    title(['corr CPUE All with ' tanom{j}])
    %stamp('')
    print('-dpng',[ppath 'Map_LMEs_pos_corr_coeff_',tanom{j},'_cpue_all.png'])

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

figure(3)
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
print('-dpng',[ppath 'Map_LMEs_corr_coeff_drivers_cpue_allfish.png'])

%% Dominant driver map - FEISTY only
figure(6)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Acorr  = nan(ni,nj);
    Acorr(id) = LAtab(i,4);

    if (LAtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Acorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',LAtab(i,1))] , ['lag=' num2str(LAtab(i,3))]};
    loc = round(length(id)/2);
    %     textm(TLAT(id(loc)),TLONG(id(loc)),str,...
    %         'HorizontalAlignment','center')
end
colormap(fcol)
colorbar('Ticks',1:8,'TickLabels',ctex)
title('CPUE corr All fishes')
stamp('CPUE corr')
%print('-dpng',[ppath 'Map_LMEs_cpue_feisty_maxcorr_allfish_v2.png'])

%% Bar plot of CPUE and each driver
figure(1)
for i=1:length(lid)
    L=lid(i);
    if (LAtab(i,2) <= 0.05)
        b=bar(L,LAtab(i,1),'EdgeColor',mcol(LAtab(i,4),:),'FaceColor',mcol(LAtab(i,4),:));
        hold on
    else
        b=bar(L,LAtab(i,1),'EdgeColor',mcol(LAtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(LAtab(i,3));
    if (LAtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1 1])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
xlabel('LME')
ylabel('Corr Coeff')
title('CPUE All fishes')
%print('-dpng',[ppath 'Bar_LMEs_cpue_sat_driver_feisty_maxcorr_allfish_v2.png'])

%%
f2 = figure('Units','inches','Position',[1 3 7.5 10]);
subplot('Position',[0.1 0.24 0.856 0.17])
% Get fake colors first for legend
for i=1:length(ctex)
    L=lid(i);
    b=bar(L,LAtab(i,1),'EdgeColor','none','FaceColor',mcol(i,:));
    hold on
end
%legend of colors and shapes
lgd = legend(ctex,'Location','eastoutside');
lgd.AutoUpdate = 'off';
for i=1:length(lid)
    L=lid(i);
    if (LAtab(i,2) <= 0.05)
        b=bar(L,LAtab(i,1),'EdgeColor',mcol(LAtab(i,4),:),'FaceColor',mcol(LAtab(i,4),:));
        hold on
    else
        b=bar(L,LAtab(i,1),'EdgeColor',mcol(LAtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(LAtab(i,3));
    if (LAtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('All fishes')

subplot('Position',[0.1 0.81 0.7 0.17])
for i=1:length(lid)
    L=lid(i);
    if (LFtab(i,2) <= 0.05)
        b=bar(L,LFtab(i,1),'EdgeColor',mcol(LFtab(i,4),:),'FaceColor',mcol(LFtab(i,4),:));
        hold on
    else
        b=bar(L,LFtab(i,1),'EdgeColor',mcol(LFtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(LFtab(i,3));
    if (LFtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('Forage fishes')
title('CPUE corr')

subplot('Position',[0.1 0.62 0.7 0.17])
for i=1:length(lid)
    L=lid(i);
    if(i~=61)
        if (LPtab(i,2) <= 0.05)
            b=bar(L,LPtab(i,1),'EdgeColor',mcol(LPtab(i,4),:),'FaceColor',mcol(LPtab(i,4),:));
            hold on
        else
            b=bar(L,LPtab(i,1),'EdgeColor',mcol(LPtab(i,4),:),'FaceColor',[1 1 1]);
            hold on
        end
        xtips1 = b(1).XEndPoints;
        ytips1 = b(1).YEndPoints;
        labels = num2str(LPtab(i,3));
        if (LPtab(i,1) < 0)
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','top')
        else
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
        end
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('Large pelagics')

subplot('Position',[0.1 0.43 0.7 0.17])
for i=1:length(lid)
    L=lid(i);
    if (LDtab(i,2) <= 0.05)
        b=bar(L,LDtab(i,1),'EdgeColor',mcol(LDtab(i,4),:),'FaceColor',mcol(LDtab(i,4),:));
        hold on
    else
        b=bar(L,LDtab(i,1),'EdgeColor',mcol(LDtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(LDtab(i,3));
    if (LDtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('Demersals')
stamp('CPUE corr')

%print('-dpng',[ppath 'Bar_LMEs_cpue_sat_driver_feisty_maxcorr_fntypes.png'])

