% Plot max corr of driver-fish corrs
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
%spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',cfile,'/corrs/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';
load([spath,'LMEs_corr_driver_maxcorrs.mat'])

%%  ---------------------------------------------------------
cnam = {'coef','p','lag','idriver','driver'};
ctex = {'TP','TB','Det','ZmLoss'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

% mcol = [238/255 102/255 119/255;... %red - TP
%     0/255 68/255 136/255;...    %blue - TB
%     34/255 136/255 51/255;...   %green - Det
%     51/255 187/255 238/255;...  %cyan - ZL
%     ];

load('paul_tol_cmaps.mat')

mcol(1,:) = drainbow(12,:)/255; % orange
mcol(2,:) = drainbow(4,:)/255; %dk blue
mcol(3,:) = drainbow(15,:)/255; %grey 
mcol(4,:) = drainbow(6,:)/255; %lt blue

%colorblind friendly
% cb = [238/255 102/255 119/255;... %red
%     170/255 51/255 119/255;...  %purple
%     34/255 136/255 51/255;...   %green
%     0/255 68/255 136/255;...    %blue
%     51/255 187/255 238/255;...  %cyan
%     ];

%%
figure(1)
% Get fake colors first for legend
for i=1:4
    L=lid(i);
    b=bar(L,LAtab(i,1),'EdgeColor','none','FaceColor',mcol(i,:));
    hold on
end
%legend of colors and shapes
lgd = legend(ctex,'Location','southeast');
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
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
xlabel('LME')
ylabel('Biomass Corr Coeff')
title('All fishes')
print('-dpng',[ppath 'Bar_LMEs_driver_biom_maxcorr_allfish.png'])

%%
f2 = figure('Units','inches','Position',[1 3 7.5 10]);
subplot('Position',[0.1 0.24 0.856 0.17])
% Get fake colors first for legend
for i=1:4
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
title('Biomass corr')

subplot('Position',[0.1 0.62 0.7 0.17])
for i=1:length(lid)
    L=lid(i);
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

subplot('Position',[0.1 0.05 0.7 0.17])
for i=1:length(lid)
    L=lid(i);
    if (LBtab(i,2) <= 0.05)
        b=bar(L,LBtab(i,1),'EdgeColor',mcol(LBtab(i,4),:),'FaceColor',mcol(LBtab(i,4),:));
        hold on
    else
        b=bar(L,LBtab(i,1),'EdgeColor',mcol(LBtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(LBtab(i,3));
    if (LBtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
xlabel('LME')
ylabel('Benthos')

print('-dpng',[ppath 'Bar_LMEs_driver_biom_maxcorr_fntypes.png'])

%% Map
%cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
cpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';

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
clatlim=[plotminlat plotmaxlat];
clonlim=[plotminlon plotmaxlon];

load coastlines;

%%
% on grid
Fcorr  = nan(ni,nj);
Pcorr  = nan(ni,nj);
Dcorr  = nan(ni,nj);
Acorr  = nan(ni,nj);
Bcorr  = nan(ni,nj);

%significance
Fsig  = ones(ni,nj);
Psig  = ones(ni,nj);
Dsig  = ones(ni,nj);
Vsig  = ones(ni,nj);
Bsig  = ones(ni,nj);

%%
figure(3)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
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
colormap(mcol)
colorbar('Ticks',1:5,'TickLabels',ctex)
title('Biomass corr All fishes')
stamp('Biomass corr')
print('-dpng',[ppath 'Map_LMEs_driver_biom_maxcorr_allfish.png'])

%%
f4 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Fcorr(id) = LFtab(i,4);

    if (LFtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Fcorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',LFtab(i,1))] , ['lag=' num2str(LFtab(i,3))]};
    %     textm(TLAT(id(1)),TLONG(id(1)),str,...
    %         'HorizontalAlignment','center')
    %     textm(TLAT(id(2)),TLONG(id(2)),['lag=' num2str(LFtab(i,3))],...
    %         'HorizontalAlignment','center')
end
colormap(mcol)
title('Forage fishes')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Pcorr(id) = LPtab(i,4);

    if (LPtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Pcorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',LPtab(i,1))] , ['lag=' num2str(LPtab(i,3))]};
end
colormap(mcol)
title('Large pelagics')

subplot('Position',[0.65 0.575 0.32 0.4]) %A
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
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
colormap(mcol)
title('All fishes')

subplot('Position',[0.01 0.10 0.32 0.4]) %D
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Dcorr(id) = LDtab(i,4);

    if (LDtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Dcorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',LDtab(i,1))] , ['lag=' num2str(LDtab(i,3))]};

end
colormap(mcol)
title('Demersals')

subplot('Position',[0.33 0.10 0.32 0.4]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Bcorr(id) = LBtab(i,4);

    if (LBtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Bcorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',LBtab(i,1))] , ['lag=' num2str(LBtab(i,3))]};

end
colormap(mcol)
colorbar('Position',[0.7 0.2 0.03 0.25],'Ticks',1:5,'TickLabels',ctex)
title('Benthos')
stamp('Biomass corr')
print('-dpng',[ppath 'Map_LMEs_driver_biom_maxcorr_fntypes.png'])


%% Correlation value map
cmR = cbrewer('seq','Reds',9,'PCHIP');
cmRB = cbrewer('div','RdBu',21,'PCHIP');
cmRB = flipud(cmRB);

% put on grid only significant ones
Rcorr  = nan(ni,nj);
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    if (LAtab(i,2) <= 0.05)
        Rcorr(id) = LAtab(i,1);
    end
end

%%
figure(7)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Rcorr)
hold on
colormap(cmRB)
colorbar
clim([-1 1])
title('Biomass corr coeff All fishes')
stamp('Biomass corr')
print('-dpng',[ppath 'Map_LMEs_driver_biom_maxcorr_coeff_allfish_v2.png'])

figure(8)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Rcorr.^2)
hold on
colormap(cmR)
colorbar
clim([0 0.9])
title('Biomass R^2 All fishes')
stamp('Biomass R^2')
print('-dpng',[ppath 'Map_LMEs_driver_biom_maxcorrR2_allfish_v2.png'])



