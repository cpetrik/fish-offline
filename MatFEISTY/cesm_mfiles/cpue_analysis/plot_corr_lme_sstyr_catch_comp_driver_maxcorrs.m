% Plot corrs of ind driver combos with catch
% Subplots together for comparison
% Restricted analysis to sst yrs

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regress_cpue/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',...
    cfile,'/corrs_cpue/'];

mod = 'v15_All_fish03';
mod2 = 'v15_obsfish';

%%
load([spath,'LMEs_corr_catch_sstyrs_inputs_feisty_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab','lid')

AtabC = LAtab;
FtabC = LFtab;
PtabC = LPtab;
DtabC = LDtab;

clear LAtab LFtab LPtab LDtab

%%
load([spath,'LMEs_corr_catch_sstyrs_inputs_obsfish_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

AtabO = LAtab;
FtabO = LFtab;
PtabO = LPtab;
DtabO = LDtab;

clear LAtab LFtab LPtab LDtab

%%
load([spath,'LMEs_corr_catch_sstyrs_inputs_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

AtabD = LAtab(lid,:);
FtabD = LFtab(lid,:);
PtabD = LPtab(lid,:);
DtabD = LDtab(lid,:);

clear LAtab LFtab LPtab LDtab

%%
load([spath,'LMEs_corr_catch_sstyrs_maxcorr.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

AtabS = LAtab(lid,:);
FtabS = LFtab(lid,:);
PtabS = LPtab(lid,:);
DtabS = LDtab(lid,:);

sst = (AtabS(:,4)==1);
AtabS(sst,4) = 5;

clear LAtab LFtab LPtab LDtab

%%  ---------------------------------------------------------
cnam = {'coef','p','lag','idriver','driver'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

%% colorblind friendly
load('paul_tol_cmaps.mat')

%colorblind friendly - subselect & re-order drainbow
ctex = {'TP','TB','Det','ZmLoss','SST','Biom','Prod','Yield'};
% orange, dk blue, grey, lt blue, dk purp, lt purp, red, green
mcol(1,:) = drainbow(12,:)/255; % orange
mcol(2,:) = drainbow(4,:)/255; %dk blue
mcol(3,:) = drainbow(15,:)/255; %grey
mcol(4,:) = drainbow(6,:)/255; %lt blue
mcol(5,:) = drainbow(14,:)/255; %red
mcol(6,:) = drainbow(3,:)/255; %dk purp
mcol(7,:) = drainbow(1,:)/255; %lt purp
mcol(8,:) = drainbow(9,:)/255; %lt green

ccol = mcol(1:7,:);

scol = mcol(5,:);

dcol = mcol(1:5,:);

%% Barplot 4x1
f1 = figure('Units','inches','Position',[1 3 7.5 10]);

%Bottom
subplot('Position',[0.1 0.05 0.856 0.22])
% Get fake colors first for legend
for i=1:length(ctex)
    L=lid(i);
    b=bar(L,AtabO(i,1),'EdgeColor','none','FaceColor',mcol(i,:));
    hold on
end
%legend of colors and shapes
lgd = legend(ctex,'Location','eastoutside');
lgd.AutoUpdate = 'off';
for i=1:length(lid)
    L=lid(i);
    if (AtabO(i,2) <= 0.05)
        b=bar(L,AtabO(i,1),'EdgeColor',mcol(AtabO(i,4),:),'FaceColor',mcol(AtabO(i,4),:));
        hold on
    else
        b=bar(L,AtabO(i,1),'EdgeColor',mcol(AtabO(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(AtabO(i,3));
    if (AtabO(i,1) < 0)
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
ylabel('Obs Effort')

%Top
subplot('Position',[0.1 0.75 0.7 0.22])
for i=1:length(lid)
    L=lid(i);
    if (AtabS(i,2) <= 0.05)
        b=bar(L,AtabS(i,1),'EdgeColor',mcol(AtabS(i,4),:),'FaceColor',mcol(AtabS(i,4),:));
        hold on
    else
        b=bar(L,AtabS(i,1),'EdgeColor',mcol(AtabS(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(AtabS(i,3));
    if (AtabS(i,1) < 0)
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
ylabel('Sat')
title('Catch corr')

%2nd
subplot('Position',[0.1 0.52 0.7 0.22])
for i=1:length(lid)
    L=lid(i);
    if (AtabD(i,2) <= 0.05)
        b=bar(L,AtabD(i,1),'EdgeColor',mcol(AtabD(i,4),:),'FaceColor',mcol(AtabD(i,4),:));
        hold on
    else
        b=bar(L,AtabD(i,1),'EdgeColor',mcol(AtabD(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(AtabD(i,3));
    if (AtabD(i,1) < 0)
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
ylabel('Sat+OBGC')

%3rd
subplot('Position',[0.1 0.29 0.7 0.22])
for i=1:length(lid)
    L=lid(i);
    if (AtabC(i,2) <= 0.05)
        b=bar(L,AtabC(i,1),'EdgeColor',mcol(AtabC(i,4),:),'FaceColor',mcol(AtabC(i,4),:));
        hold on
    else
        b=bar(L,AtabC(i,1),'EdgeColor',mcol(AtabC(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(AtabC(i,3));
    if (AtabC(i,1) < 0)
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
ylabel('Const Effort')
stamp('')

print('-dpng',[ppath 'Bar_LMEs_sstyr_catch_driver_comp_maxcorr_allfish.png'])

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
f2 = figure('Units','inches','Position',[1 3 7.5 5]);

axs = subplot('Position',[0.1 0.525 0.35 0.45]); %Top L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Acorr  = nan(ni,nj);
    Acorr(id) = AtabS(i,4);

    if (AtabS(i,2) <= 0.05)
        surfm(TLAT,TLONG,Acorr)
        hold on
    end
end
colormap(axs,scol);
title('Catch Sat')

axd = subplot('Position',[0.5 0.525 0.35 0.45]); %Top R
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Acorr  = nan(ni,nj);
    Acorr(id) = AtabD(i,4);

    if (AtabD(i,2) <= 0.05)
        surfm(TLAT,TLONG,Acorr)
        hold on
    end
end
colormap(axd, dcol);
title('Catch Sat+OBGC')

axc = subplot('Position',[0.1 0.05 0.35 0.45]); %Bottom L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Acorr  = nan(ni,nj);
    Acorr(id) = AtabC(i,4);

    if (AtabC(i,2) <= 0.05)
        surfm(TLAT,TLONG,Acorr)
        hold on
    end
end
colormap(axc,ccol)
title('Catch Const Effort')

axo = subplot('Position',[0.5 0.05 0.43 0.45]); %Bottom R
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Acorr  = nan(ni,nj);
    Acorr(id) = AtabO(i,4);

    if (AtabO(i,2) <= 0.05)
        surfm(TLAT,TLONG,Acorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',AtabO(i,1))] , ['lag=' num2str(AtabO(i,3))]};
    loc = round(length(id)/2);
end
colormap(axo,mcol)
colorbar('eastoutside','Ticks',1:9,'TickLabels',ctex)
%colorbar('Ticks',1:9,'TickLabels',ctex)
title('Catch Obs Effort')

print('-dpng',[ppath 'Map_LMEs_sstyr_catch_driver_comp_maxcorr_allfish.png'])

