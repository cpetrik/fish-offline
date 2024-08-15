% Plot max corr of driver-catch corrs
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

load([spath,'LMEs_corr_catch_sat_driver_feisty_maxcorrs.mat'])

%%  ---------------------------------------------------------
cnam = {'coef','p','lag','idriver','driver'};
ctex = {'TP','TB','Det','ZmLoss','Biom','Prod','Rec','SST','Chl'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

% colorblind friendly
load('paul_tol_cmaps.mat')

%mcol = drainbow(1:9,:)/255;

for i=1:15
    plot(1:5,i*(1:5),'color',drainbow(i,:)/255,'LineWidth',2); hold on;
end

% mcol = [238/255 102/255 119/255;... %red - TP
%     0/255 68/255 136/255;...    %blue - TB
%     34/255 136/255 51/255;...   %green - Det
%     51/255 187/255 238/255;...  %cyan - ZL
%     ];


%colorblind friendly - subselect & re-order drainbow
%ctex = {'TP','TB','Det','ZmLoss','Biom','Prod','Rec','SST','chl'};
% orange, dk blue, dk purp, lt green, grey, lt blue, lt purp, red, green
mcol(1,:) = drainbow(12,:)/255;
mcol(2,:) = drainbow(4,:)/255;
mcol(3,:) = drainbow(3,:)/255;
mcol(4,:) = drainbow(9,:)/255;
mcol(5,:) = drainbow(15,:)/255;
mcol(6,:) = drainbow(6,:)/255;
mcol(7,:) = drainbow(1,:)/255;
mcol(8,:) = drainbow(14,:)/255;
mcol(9,:) = drainbow(7,:)/255;

%%
close all
figure(1)
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
ylim([-1 1])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
xlabel('LME')
ylabel('Corr Coeff')
title('Catch All fishes')
print('-dpng',[ppath 'Bar_LMEs_catch_sat_driver_feisty_maxcorr_allfish.png'])


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

%significance
Vsig  = ones(ni,nj);
 
%%
figure(2)
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
colormap(mcol)
colorbar('Ticks',1:9,'TickLabels',ctex)
title('Catch corr All fishes')
stamp('Catch ')
print('-dpng',[ppath 'Map_LMEs_catch_sat_driver_feisty_maxcorr_allfish.png'])

%% Correlation value map
cmR = cbrewer('seq','Reds',10,'PCHIP');
cmRB = cbrewer('div','RdBu',21,'PCHIP');

% put on grid only significant ones
Rcorr  = nan(ni,nj);
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    if (LAtab(i,2) <= 0.05)
        Rcorr(id) = LAtab(i,1);
    end
end

figure(3)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Rcorr)
hold on
colormap(cmRB)
colorbar
caxis([-1 1])
title('Catch corr coeff All fishes')
stamp('Catch corr')
print('-dpng',[ppath 'Map_LMEs_catch_sat_driver_feisty_maxcorr_coeff_allfish.png'])

figure(4)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Rcorr.^2)
hold on
colormap(cmR)
colorbar
caxis([0 1])
title('Catch R^2 All fishes')
stamp('Catch ')
print('-dpng',[ppath 'Map_LMEs_catch_sat_driver_feisty_maxcorr_r2_allfish.png'])




