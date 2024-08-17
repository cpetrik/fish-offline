% Plot corrs of driver-catch corrs
% Barplots and maps
% Lag with max R2
% For all 63 LMEs
% Obs fishing effort

clear
close all

%% ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regress_cpue/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',...
    cfile,'/corrs_cpue/'];

%mod = 'v15_All_fish03';
mod = 'v15_obsfish';

%%
load([spath,'LMEs_corr_catch_chlyrs_inputs_obsfish_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab','lid')

cnam = {'coef','p','lag','idriver','driver'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

Atab = LAtab;
Ftab = LFtab;
Ptab = LPtab;
Dtab = LDtab;

% Catch of Lg Pel in LME 61 is 0, fix nans
LPtab(61,1) = 0;
LPtab(61,2) = 1;
%LPtab(61,3) = 0;
LPtab(61,4) = 0;

%% colorblind friendly
load('paul_tol_cmaps.mat')

% mcol = [238/255 102/255 119/255;... %red - TP
%     0/255 68/255 136/255;...    %blue - TB
%     34/255 136/255 51/255;...   %green - Det
%     51/255 187/255 238/255;...  %cyan - ZL
%     ];

%colorblind friendly - subselect & re-order drainbow
ctex = {'TP','TB','Det','ZmLoss','SST','Chl','Biom','Prod','Yield'};
% orange, dk blue, grey, lt blue, dk purp, lt purp, red, green
mcol(1,:) = drainbow(12,:)/255; % orange
mcol(2,:) = drainbow(4,:)/255; %dk blue
mcol(3,:) = drainbow(15,:)/255; %grey
mcol(4,:) = drainbow(6,:)/255; %lt blue
mcol(5,:) = drainbow(14,:)/255; %red
mcol(6,:) = drainbow(7,:)/255; %green
mcol(7,:) = drainbow(3,:)/255; %dk purp
mcol(8,:) = drainbow(1,:)/255; %lt purp
mcol(9,:) = drainbow(9,:)/255; %lt green

%% 
f1 = figure('Units','inches','Position',[1 3 7.5 10]);

subplot('Position',[0.1 0.52 0.856 0.22])
% Get fake colors first for legend
for i=1:length(ctex)
    L=lid(i);
    b=bar(L,LPtab(i,1),'EdgeColor','none','FaceColor',mcol(i,:));
    hold on
end
%legend of colors and shapes
lgd = legend(ctex,'Location','eastoutside');
lgd.AutoUpdate = 'off';
for i=1:length(lid)
    L=lid(i);
    if(i~=61)
        if(i~=63)
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
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('Large pelagics')

subplot('Position',[0.1 0.75 0.7 0.22])
for i=1:length(lid)
    L=lid(i);
    if(i~=55)
        if(i~=61)
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
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('Forage fishes')
title('Catch corr')

subplot('Position',[0.1 0.29 0.7 0.22])
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

subplot('Position',[0.1 0.05 0.7 0.22])
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
ylabel('All fishes')

print('-dpng',[ppath 'Bar_LMEs_chlyr_catch_driver_obsfish_maxcorr_fntypes.png'])

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
clatlim=[plotminlat plotmaxlat];
clonlim=[plotminlon plotmaxlon];

load coastlines;


%% Dominant driver map
f6 = figure('Units','inches','Position',[1 3 7.5 5]);

subplot('Position',[0.05 0.555 0.45 0.4]) %Top L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Fcorr  = nan(ni,nj);
    Fcorr(id) = LFtab(i,4);

    if (LFtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Fcorr)
        hold on
    end
    
end
colormap(mcol)
title('Catch corr Forage')

subplot('Position',[0.5 0.555 0.45 0.4]) %Top R
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Pcorr  = nan(ni,nj);
    Pcorr(id) = LPtab(i,4);

    if (LPtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Pcorr)
        hold on
    end
    
end
colormap(mcol)
title('Catch corr Lg Pelagics')

subplot('Position',[0.05 0.105 0.45 0.4]) %Bottom L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Dcorr  = nan(ni,nj);
    Dcorr(id) = LDtab(i,4);

    if (LDtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Dcorr)
        hold on
    end
    
end
colormap(mcol)
title('Catch corr Demersals')

subplot('Position',[0.5 0.105 0.45 0.4]) %Bottom R
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
    
end
colormap(mcol)
colorbar('TickLabels',ctex,'Position',[0.25 0.055 0.5 0.03],'orientation','horizontal')
title('Catch corr All fishes')

print('-dpng',[ppath 'Map_LMEs_chlyr_catch_driver_obsfish_maxcorr_fntypes.png'])



