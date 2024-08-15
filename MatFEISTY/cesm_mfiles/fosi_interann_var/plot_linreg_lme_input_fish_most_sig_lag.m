% Plot max coeffs of driver-fish regressions
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/' cfile '/corrs/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

load([spath,'LMEs_regress_driver_ALLdiv2SD_maxcoefs.mat'])

%%  ---------------------------------------------------------
cnam = {'coef','p','lag','idriver','driver'};
ctex = {'TP','TB','Det','Zmeso','ZmLoss'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

% 1= AMO, 2=AO, 3=NAO, 4=Nino34, 5=PDO
mtype = {'o','s','^','v','d','p'}; %for lag
mcol = [238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    34/255 136/255 51/255;...   %green
    0/255 68/255 136/255;...    %blue
    51/255 187/255 238/255;...  %cyan
    ];


%colorblind friendly
% cb=[34/255 136/255 51/255;...   %green
%     153/255 153/255 51/255;...  %olive
%     51/255 187/255 238/255;...  %cyan
%     0/255 68/255 136/255;...    %blue
%     238/255 102/255 119/255;... %red
%     170/255 51/255 119/255;...  %purple
%     0 0 0;...                   %black
%     0.25 0.25 0.25;...             %dk grey
%     0.50 0.50 0.50;...             % grey
%     0.75 0.75 0.75];               %lt grey

%%
figure(1)
% Get fake colors first for legend
for i=1:5
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
ylabel('Coeff')
title('All fishes')
print('-dpng',[ppath 'Bar_LMEs_driver_lrcoef_allfish.png'])

%%
f2 = figure('Units','inches','Position',[1 3 7.5 10]);
subplot('Position',[0.1 0.24 0.856 0.17])
% Get fake colors first for legend
for i=1:5
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

print('-dpng',[ppath 'Bar_LMEs_driver_lrcoef_fntypes.png'])

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
title('All fishes')
print('-dpng',[ppath 'Map_LMEs_driver_lrcoef_allfish.png'])

%%
f1 = figure('Units','inches','Position',[1 3 7.5 5]);
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
print('-dpng',[ppath 'Map_LMEs_driver_lrcoef_fntypes.png'])

%%
% on grid
Fcorr  = nan(ni,nj);
Pcorr  = nan(ni,nj);
Dcorr  = nan(ni,nj);
Acorr  = nan(ni,nj);
Bcorr  = nan(ni,nj);

%% coef = color, text of driver and lag
figure(4)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Acorr(id) = LAtab(i,1);

    if (LAtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Acorr)
        hold on
    end
    str = {ctex{LAtab(i,4)} , ['lag=' num2str(LAtab(i,3))]};
    loc = round(length(id)/2);
    if (abs(LAtab(i,1)) >= 0.2)
        textm(TLAT(id(loc)),TLONG(id(loc)),str,...
            'HorizontalAlignment','center')
    end
end
cmocean('balance')
caxis([-1 1])
colorbar
title('All fishes')
print('-dpng',[ppath 'Map_LMEs_driver_lrcoef_allfish_v2.png'])

%% coef = color, colormap varies by driver
cmRB=cbrewer('div','RdBu',11,'PCHIP'); %tp
cmPY=cbrewer('div','PiYG',11,'PCHIP'); %zoo
cmPG=cbrewer('div','PRGn',11,'PCHIP'); %zlos
cmBB=cbrewer('div','BrBG',11,'PCHIP'); %det
cmPO=cbrewer('div','PuOr',11,'PCHIP'); %tb


tpi = find(LAtab(:,4)==1);
tbi = find(LAtab(:,4)==2);
dti = find(LAtab(:,4)==3);
zmi = find(LAtab(:,4)==4);
zli = find(LAtab(:,4)==5);

Acorr  = nan(ni,nj);
Dcorr  = nan(ni,nj);
Zcorr  = nan(ni,nj);

%%
for i=1:length(lid)
    L=lid((i));
    id = find(tlme==L);
    ALLcorr(id) = LAtab(i,1);
end

figure(5)
%TP
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
for i=1:length(tpi)
    L=lid(tpi(i));
    id = find(tlme==L);
    Acorr(id) = LAtab(i,1);
    if (LAtab(i,2) <= 0.05)
        h(1) = surfm(TLAT,TLONG,Acorr); hold on;
    end
end
%colormap(ax3,cmRB)
%Det
for i=1:length(dti)
    L=lid(dti(i));
    id = find(tlme==L);
    Dcorr(id) = LAtab(i,1);
    if (LAtab(i,2) <= 0.05)
        h(2) = surfm(TLAT,TLONG,Dcorr); hold on;
    end
end
%colormap(ax3,cmBB)
%Zm
for i=1:length(zmi)
    L=lid(zmi(i));
    id = find(tlme==L);
    Zcorr(id) = LAtab(i,1);
    if (LAtab(i,2) <= 0.05)
        h(3) = surfm(TLAT,TLONG,Zcorr); hold on;
    end
end
%colormap(ax3,cmPY)
cmap = [cmRB;cmBB;cmPY];
colormap(cmap)
tl = 1:-0.2:-1;
tl3=[tl,tl,tl];
cb = colorbar;
cb.Ticks = -1:0.062:1;
cb.TickLabels = tl3;
%,1:2:33,'TickLabels',tl3(1):2:tl3(end))

%% 3 colorpmarp
% Project Data to the different planes.
h(1) = surf(X,Y,ax(5)+0*Z,Z); % X-Y Plane Projection
h(2) = surf(X,ax(4)+0*Y,Z);   % X-Z Plane Projection
h(3) = surf(ax(2)+0*X,Y,Z);   % Y-Z Plane Projection
set(h,'FaceColor','interp','EdgeColor','interp')
% Build a colormap that consists of three separate
% colormaps.
cmapX = bone(32);
cmapY = cool(32);
cmapZ = jet(32);
cmap = [cmapX;cmapY;cmapZ];
colormap(cmap)
% Map the CData of each surface plot to a contiguous, 
% nonoverlapping set of data.  Each CData must have
% the same range.
zmin = min(Z(:));
zmax = max(Z(:));
% CDX ranges from 1 to 32.
cdx = min(32,round(31*(Z-zmin)/(zmax-zmin))+1);
% CDY ranges from 33 to 64.
cdy = cdx+32;
% CDZ ranges from 65 to 96.
cdz = cdy+32;
% Update the CDatas.
set(h(1),'CData',cdx)
set(h(3),'CData',cdy)
set(h(2),'CData',cdz)
% Change CLim (Color Limits) so that it spans all the CDatas
caxis([min(cdx(:)) max(cdz(:))])


title('All fishes')
%print('-dpng',[ppath 'Map_LMEs_driver_lrcoef_allfish_v3.png'])
