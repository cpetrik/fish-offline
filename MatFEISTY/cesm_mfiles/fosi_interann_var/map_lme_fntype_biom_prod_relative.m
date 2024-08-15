% Map fish biomass and production
% Relative to other fn types
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/' cfile '/corrs/'];

mod = 'v15_All_fish03_';

%% inputs
load([cpath 'lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'],...
    'lme_tp_fosi','lme_tb_fosi','lme_det_fosi','lme_mz_fosi','lme_mzloss_fosi')

% fish
load([fpath 'LME_fosi_fished_',mod,cfile '.mat'],'lme_area','lme_mtype',...
    'lme_mprod');

adult = [4;7;8];
lme_prod = lme_mprod(:,adult);

%%  ---------------------------------------------------------
ftex = {'F','P','D','A','B'};
ctex = {'TP','TB','Det','Zmeso','ZmLoss'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

% US LMEs
lid = [54,1:2,10,3,5:7,65]; %ADD 65 = Aleutian Islands
lname = {'CHK','EBS','GAK','HI','CCE','GMX','SE','NE','AI'};

% mcol = [238/255 102/255 119/255;... %red
%     170/255 51/255 119/255;...  %purple
%     34/255 136/255 51/255;...   %green
%     0/255 68/255 136/255;...    %blue
%     51/255 187/255 238/255;...  %cyan
%     ];

mcol = [238/255 102/255 119/255;... %red
    0/255 68/255 136/255;...    %blue
    34/255 136/255 51/255;...   %green
    51/255 187/255 238/255;...  %cyan
    170/255 51/255 119/255;...  %purple
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

set(groot,'defaultAxesColorOrder',mcol);

cmB=cbrewer('seq','Blues',50,'PCHIP');
cmR=cbrewer('seq','Reds',50,'PCHIP');
cmG=cbrewer('seq','Greens',50,'PCHIP');

%% Relative amounts
lme_type = lme_mtype(:,1:3);
lme_btot = sum(lme_type,2);
lme_ptot = sum(lme_prod,2);

lbio = lme_type ./ repmat(lme_btot,1,3);
lprod = lme_prod ./ repmat(lme_ptot,1,3);

%% Biomass & prod bar plots
figure(1)
subplot(3,2,1)
bar(lme_mtype(lid,1:3))
ylabel('Biomass (gWW/m2)')
set(gca,'XTickLabel',lname)
legend(ftex(1:3))
legend('location','northeast')

subplot(3,2,2)
bar(lme_prod(lid,:))
ylabel('Productivity (gWW/m2/y)')
set(gca,'XTickLabel',lname)

subplot(3,2,3)
bar(lme_mtype(lid,1:3),'stacked')
ylabel('Biomass')
set(gca,'XTickLabel',lname)

subplot(3,2,4)
bar(lme_prod(lid,:),'stacked')
ylabel('Productivity')
set(gca,'XTickLabel',lname)

subplot(3,2,5)
bar(lbio(lid,:),'stacked')
ylabel('Biomass')
set(gca,'XTickLabel',lname)
ylim([0 1])

subplot(3,2,6)
bar(lprod(lid,:),'stacked')
ylabel('Productivity')
set(gca,'XTickLabel',lname)
ylim([0 1])
print('-dpng',[ppath 'Bar_USLMEs_biom_prod_fntypes.png'])

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
Fbiom  = nan(ni,nj);
Pbiom  = nan(ni,nj);
Dbiom  = nan(ni,nj);

Fprod  = nan(ni,nj);
Pprod  = nan(ni,nj);
Dprod  = nan(ni,nj);

for i=1:66
    L=i;
    id = find(tlme==L);

    Fbiom(id) = lbio(i,1);
    Pbiom(id) = lbio(i,2);
    Dbiom(id) = lbio(i,3);
    Fprod(id) = lprod(i,1);
    Pprod(id) = lprod(i,2);
    Dprod(id) = lprod(i,3);
end

%%
f1 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
ax1=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Fbiom)
colormap(ax1,cmR)
caxis([0 1]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('Biomass F')
%title('Forage fishes')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
ax2=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Pbiom)
colormap(ax2,cmB)
caxis([0 1]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('P')

subplot('Position',[0.65 0.575 0.32 0.4]) %D
ax3=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Dbiom)
colormap(ax3,cmG)
caxis([0 1]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('D')

subplot('Position',[0.01 0.10 0.32 0.4]) %F Prod
ax4=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Fprod)
colormap(ax4,cmR)
caxis([0 1]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('Productivity F')
colorbar('Position',[0.01 0.075 0.3 0.025],'orientation','horizontal')

subplot('Position',[0.33 0.10 0.32 0.4]) %
ax5=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Pprod)
colormap(ax5,cmB)
caxis([0 1]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('P')
colorbar('Position',[0.33 0.075 0.3 0.025],'orientation','horizontal')

subplot('Position',[0.65 0.10 0.32 0.4])
ax6=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Dprod)
colormap(ax6,cmG)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('D')
caxis([0 1]);
colorbar('Position',[0.65 0.075 0.3 0.025],'orientation','horizontal')

print('-dpng',[ppath 'Map_LMEs_rel_biom_prod_fntypes.png'])

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
print('-dpng',[ppath 'Map_LMEs_driver_corrcoef_allfish_v2.png'])

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
%print('-dpng',[ppath 'Map_LMEs_driver_corrcoef_allfish_v3.png'])
