% Plot max SLR coeff of driver-fish biom & nu corrs together
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',cfile,'/corrs/'];

mod = 'v15_All_fish03';

load([spath,'LMEs_SLR_biom_driver_maxR2s.mat'],'LAtab','LFtab','LPtab','LDtab',... %table of coeff of best driver at best lag
    'lid')

BAtab = LAtab;
BFtab = LFtab;
BPtab = LPtab;
BDtab = LDtab;

clear LAtab LFtab LPtab LDtab

%%
load([spath,'LMEs_SLR_nu_driver_maxR2s.mat'],'LAtab','LFtab','LPtab','LDtab') %table of coeff of best driver at best lag

PAtab = LAtab;
PFtab = LFtab;
PPtab = LPtab;
PDtab = LDtab;

clear LAtab LFtab LPtab LDtab

%%  ---------------------------------------------------------
cnam = {'coef','p','R2','lag','idriver','driver'};
ctex = {'TP','TB','Det','Mzoo'};
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

%%
f1 = figure('Units','inches','Position',[1 3 7.5 10]);

% Biomass
subplot('Position',[0.1 0.74 0.4 0.22])
for i=1:length(lid)
    L=lid(i);
    if (BAtab(i,2) <= 0.05)
        b=bar(L,BAtab(i,1),'EdgeColor',mcol(BAtab(i,5),:),'FaceColor',mcol(BAtab(i,5),:));
        hold on
    else
        b=bar(L,BAtab(i,1),'EdgeColor',mcol(BAtab(i,5),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(BAtab(i,4));
    if (BAtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:5:66,'XTickLabel','')
ylabel('All fishes')
title('Biomass')

subplot('Position',[0.1 0.51 0.4 0.22])
for i=1:length(lid)
    L=lid(i);
    if (BFtab(i,2) <= 0.05)
        b=bar(L,BFtab(i,1),'EdgeColor',mcol(BFtab(i,5),:),'FaceColor',mcol(BFtab(i,5),:));
        hold on
    else
        b=bar(L,BFtab(i,1),'EdgeColor',mcol(BFtab(i,5),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(BFtab(i,4));
    if (BFtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:5:66,'XTickLabel','')
ylabel('Forage fishes')

subplot('Position',[0.1 0.28 0.4 0.22])
for i=1:length(lid)
    L=lid(i);
    if (BPtab(i,2) <= 0.05)
        b=bar(L,BPtab(i,1),'EdgeColor',mcol(BPtab(i,5),:),'FaceColor',mcol(BPtab(i,5),:));
        hold on
    else
        b=bar(L,BPtab(i,1),'EdgeColor',mcol(BPtab(i,5),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(BPtab(i,4));
    if (BPtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:5:66,'XTickLabel','')
ylabel('Large pelagics')

subplot('Position',[0.1 0.05 0.4 0.22])
for i=1:length(lid)
    L=lid(i);
    if (BDtab(i,2) <= 0.05)
        b=bar(L,BDtab(i,1),'EdgeColor',mcol(BDtab(i,5),:),'FaceColor',mcol(BDtab(i,5),:));
        hold on
    else
        b=bar(L,BDtab(i,1),'EdgeColor',mcol(BDtab(i,5),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(BDtab(i,4));
    if (BDtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:5:66,'XTickLabel',1:5:66)
xlabel('LME')
ylabel('Demersals')

% Prod
subplot('Position',[0.55 0.74 0.4 0.22])
for i=1:length(lid)
    L=lid(i);
    if (PAtab(i,2) <= 0.05)
        b=bar(L,PAtab(i,1),'EdgeColor',mcol(PAtab(i,5),:),'FaceColor',mcol(PAtab(i,5),:));
        hold on
    else
        b=bar(L,PAtab(i,1),'EdgeColor',mcol(PAtab(i,5),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(PAtab(i,4));
    if (PAtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:5:66,'XTickLabel','')
title('Production')

subplot('Position',[0.55 0.51 0.4 0.22])
% Get fake colors first for legend
for i=1:4
    L=lid(i);
    b=bar(L,PAtab(i,1),'EdgeColor','none','FaceColor',mcol(i,:));
    hold on
end
%legend of colors and shapes
lgd = legend(ctex,'Location','southwest');
lgd.AutoUpdate = 'off';
for i=1:length(lid)
    L=lid(i);
    if (PFtab(i,2) <= 0.05)
        b=bar(L,PFtab(i,1),'EdgeColor',mcol(PFtab(i,5),:),'FaceColor',mcol(PFtab(i,5),:));
        hold on
    else
        b=bar(L,PFtab(i,1),'EdgeColor',mcol(PFtab(i,5),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(PFtab(i,4));
    if (PFtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:5:66,'XTickLabel','')

subplot('Position',[0.55 0.28 0.4 0.22])
for i=1:length(lid)
    L=lid(i);
    if (PPtab(i,2) <= 0.05)
        b=bar(L,PPtab(i,1),'EdgeColor',mcol(PPtab(i,5),:),'FaceColor',mcol(PPtab(i,5),:));
        hold on
    else
        b=bar(L,PPtab(i,1),'EdgeColor',mcol(PPtab(i,5),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(PPtab(i,4));
    if (PPtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:5:66,'XTickLabel','')

subplot('Position',[0.55 0.05 0.4 0.22])
for i=1:length(lid)
    L=lid(i);
    if (PDtab(i,2) <= 0.05)
        b=bar(L,PDtab(i,1),'EdgeColor',mcol(PDtab(i,5),:),'FaceColor',mcol(PDtab(i,5),:));
        hold on
    else
        b=bar(L,PDtab(i,1),'EdgeColor',mcol(PDtab(i,5),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(PDtab(i,4));
    if (PDtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:5:66,'XTickLabel',1:5:66)
xlabel('LME')
stamp('LR')

print('-dpng',[ppath 'Bar_LMEs_biom_nu_driver_SLR_maxR2s_fntypes.png'])

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
clatlim=[plotminlat plotmaxlat];
clonlim=[plotminlon plotmaxlon];

load coastlines;

%% Correlation value map
cmR = cbrewer('seq','Reds',9,'PCHIP');
cmRB = cbrewer('div','RdBu',21,'PCHIP');
cmRB = flipud(cmRB);

% put on grid only significant ones
BRcorr  = nan(ni,nj);
PRcorr  = nan(ni,nj);
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    if (BAtab(i,2) <= 0.05)
        BRcorr(id) = BAtab(i,3);
    end
    if (PAtab(i,2) <= 0.05)
        PRcorr(id) = PAtab(i,3);
    end
end

%%
f2 = figure('Units','inches','Position',[1 3 7.5 2.75]);
subplot('Position',[0.01 0.05 0.45 0.8])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,BRcorr)
hold on
colormap(cmR)
clim([0 0.9])
title('Biomass')

subplot('Position',[0.47 0.05 0.45 0.8])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,PRcorr)
hold on
colormap(cmR)
clim([0 0.9])
colorbar('Position',[0.925 0.1 0.03 0.7])
title('Production')
print('-dpng',[ppath 'Map_LMEs_driver_biom_nu_SLR_maxR2s_allfish_R2.png'])


