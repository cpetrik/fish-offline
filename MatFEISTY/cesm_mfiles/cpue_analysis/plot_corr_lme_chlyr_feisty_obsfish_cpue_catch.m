% Plot corrs of driver with cpue & catch, const & obs effort
% Subplots together for comparison
% Restricted analysis to chl yrs

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

%PC:cpue const, PO:cpue obs, CC:catch const, CO:catch obs

%%
load([spath,'LMEs_corr_cpue_chlyrs_inputs_feisty_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab','lid')

AtabPC = LAtab;
FtabPC = LFtab;
PtabPC = LPtab;
DtabPC = LDtab;

clear LAtab LFtab LPtab LDtab

%%
load([spath,'LMEs_corr_cpue_chlyrs_driver_obsfish_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

AtabPO = LAtab;
FtabPO = LFtab;
PtabPO = LPtab;
DtabPO = LDtab;

clear LAtab LFtab LPtab LDtab

%%
load([spath,'LMEs_corr_catch_chlyrs_inputs_feisty_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

AtabCC = LAtab;
FtabCC = LFtab;
PtabCC = LPtab;
DtabCC = LDtab;

clear LAtab LFtab LPtab LDtab

%%
load([spath,'LMEs_corr_catch_chlyrs_inputs_obsfish_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

AtabCO = LAtab;
FtabCO = LFtab;
PtabCO = LPtab;
DtabCO = LDtab;

clear LAtab LFtab LPtab LDtab

%%  ---------------------------------------------------------
cnam = {'coef','p','lag','idriver','driver'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

%% colorblind friendly
load('paul_tol_cmaps.mat')

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

%% Barplot 2x2
f1 = figure('Units','inches','Position',[1 3 12 6]);

subplot('Position',[0.51 0.53 0.48 0.43])
%subplot(2,2,2)
% Get fake colors first for legend
for i=1:length(ctex)
    L=lid(i);
    b=bar(L,AtabCC(i,1),'EdgeColor','none','FaceColor',mcol(i,:));
    hold on
end
%legend of colors and shapes
lgd = legend(ctex,'Location','eastoutside');
lgd.AutoUpdate = 'off';
for i=1:length(lid)
    L=lid(i);
    if (AtabCC(i,2) <= 0.05)
        b=bar(L,AtabCC(i,1),'EdgeColor',mcol(AtabCC(i,4),:),'FaceColor',mcol(AtabCC(i,4),:));
        hold on
    else
        b=bar(L,AtabCC(i,1),'EdgeColor',mcol(AtabCC(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(AtabCC(i,3));
    if (AtabCC(i,1) < 0)
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
title('Catch Const Effort')

subplot('Position',[0.035 0.53 0.39 0.43])
%subplot(2,2,1)
for i=1:length(lid)
    L=lid(i);
    if (AtabPC(i,2) <= 0.05)
        b=bar(L,AtabPC(i,1),'EdgeColor',mcol(AtabPC(i,4),:),'FaceColor',mcol(AtabPC(i,4),:));
        hold on
    else
        b=bar(L,AtabPC(i,1),'EdgeColor',mcol(AtabPC(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(AtabPC(i,3));
    if (AtabPC(i,1) < 0)
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
title('CPUE Const Effort')

subplot('Position',[0.035 0.05 0.39 0.43])
%subplot(2,2,3)
for i=1:length(lid)
    L=lid(i);
    if (AtabPO(i,2) <= 0.05)
        b=bar(L,AtabPO(i,1),'EdgeColor',mcol(AtabPO(i,4),:),'FaceColor',mcol(AtabPO(i,4),:));
        hold on
    else
        b=bar(L,AtabPO(i,1),'EdgeColor',mcol(AtabPO(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(AtabPO(i,3));
    if (AtabPO(i,1) < 0)
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
title('CPUE Obs Effort')

subplot('Position',[0.51 0.05 0.39 0.43])
%subplot(2,2,4)
for i=1:length(lid)
    L=lid(i);
    if (AtabCO(i,2) <= 0.05)
        b=bar(L,AtabCO(i,1),'EdgeColor',mcol(AtabCO(i,4),:),'FaceColor',mcol(AtabCO(i,4),:));
        hold on
    else
        b=bar(L,AtabCO(i,1),'EdgeColor',mcol(AtabCO(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(AtabCO(i,3));
    if (AtabCO(i,1) < 0)
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
title('Catch Obs Effort')

print('-dpng',[ppath 'Bar_LMEs_chlyr_cpue_catch_feisty_obsfish_maxcorr_allfish.png'])

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
mcol8 = mcol(1:8,:);

f2 = figure('Units','inches','Position',[1 3 7.5 5]);

ax1=subplot('Position',[0.1 0.525 0.35 0.45]); %Top L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Acorr  = nan(ni,nj);
    Acorr(id) = AtabPC(i,4);

    if (AtabPC(i,2) <= 0.05)
        surfm(TLAT,TLONG,Acorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',AtabPC(i,1))] , ['lag=' num2str(AtabPC(i,3))]};
    loc = round(length(id)/2);
end
colormap(ax1,mcol8)
title('CPUE Const Effort')

ax2=subplot('Position',[0.5 0.525 0.35 0.45]); %Top R
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Acorr  = nan(ni,nj);
    Acorr(id) = AtabCC(i,4);

    if (AtabCC(i,2) <= 0.05)
        surfm(TLAT,TLONG,Acorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',AtabCC(i,1))] , ['lag=' num2str(AtabCC(i,3))]};
    loc = round(length(id)/2);
end
colormap(ax2,mcol8)
title('Catch Const Effort')

ax3=subplot('Position',[0.1 0.05 0.35 0.45]); %Bottom L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Acorr  = nan(ni,nj);
    Acorr(id) = AtabPO(i,4);

    if (AtabPO(i,2) <= 0.05)
        surfm(TLAT,TLONG,Acorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',AtabPO(i,1))] , ['lag=' num2str(AtabPO(i,3))]};
    loc = round(length(id)/2);
end
colormap(ax3,mcol8)
title('CPUE Obs Effort')

ax4=subplot('Position',[0.5 0.05 0.43 0.45]); %Bottom R
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Acorr  = nan(ni,nj);
    Acorr(id) = AtabCO(i,4);

    if (AtabCO(i,2) <= 0.05)
        surfm(TLAT,TLONG,Acorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',AtabCC(i,1))] , ['lag=' num2str(AtabCC(i,3))]};
    loc = round(length(id)/2);
end
colormap(ax4,mcol)
colorbar('eastoutside','TickLabels',ctex,'Direction','reverse')
title('Catch Obs Effort')

print('-dpng',[ppath 'Map_LMEs_chlyr_cpue_catch_feisty_obsfish_pos_maxcorr_allfish.png'])

