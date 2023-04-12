% Plot max coeffs of climate-driver or climate-fish regressions
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

load([spath,'LMEs_regress_seasonal_climate_div2SD_maxcoefs.mat'])

%%  ---------------------------------------------------------
cnam = {'coef','p','lag','iclimate','climate'};
ctex = {'AMO','AO','NAO','Nino34','PDO'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

% 1= AMO, 2=AO, 3=NAO, 4=Nino34, 5=PDO
mtype = {'o','s','^','v','d','p'}; %for lag
mcol = [238/255 102/255 119/255;... %red
    34/255 136/255 51/255;...   %green
    170/255 51/255 119/255;...  %purple
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255];    %blue


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
plot(0:67,zeros(68,1),'--k'); hold on;
for i=1:length(lid)
    L=lid(i);
    if (LAtab(i,2) <= 0.05)
        plot(L,LAtab(i,1),mtype{(LAtab(i,3))+1},'Color',mcol(LAtab(i,4),:),'MarkerFaceColor',mcol(LAtab(i,4),:));
        hold on
    else
        plot(L,LAtab(i,1),mtype{(LAtab(i,3))+1},'Color',mcol(LAtab(i,4),:));
        hold on
    end
end
xlim([0 67])
ylim([-0.8 0.8])
xlabel('LME')
ylabel('Coeff')
%legend(tanom2(:,1))

%%
figure(2)
% Get fake colors first for legend
for i=1:5
    L=lid(i);
    b=bar(L,LAtab(i,1),'EdgeColor','none','FaceColor',mcol(i,:));
    hold on
end
%legend of colors and shapes
lgd = legend({'AMO','AO','NAO','Nino34','PDO'},'Location','northwest');
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
ylim([-0.8 0.8])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
xlabel('LME')
ylabel('Coeff')
title('All fishes')
print('-dpng',[ppath 'Bar_LMEs_climate_seasons_corrcoef_allfish.png'])

%%
figure(3)
subplot(3,2,1)
% Get fake colors first for legend
for i=1:5
    L=lid(i);
    b=bar(L,LFtab(i,1),'EdgeColor','none','FaceColor',mcol(i,:));
    hold on
end
%legend of colors and shapes
lgd = legend({'AMO','AO','NAO','Nino34','PDO'},'Location','northwest');
lgd.AutoUpdate = 'off';
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
ylim([-0.8 0.8])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
ylabel('Coeff')
title('Forage fishes')

subplot(3,2,3)
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
ylim([-0.8 0.8])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
ylabel('Coeff')
title('Large pelagics')

subplot(3,2,4)
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
ylim([-0.8 0.8])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
ylabel('Coeff')
title('Demersals')

subplot(3,2,5)
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
ylim([-0.8 0.8])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
ylabel('Coeff')
xlabel('LME')
title('All fishes')

subplot(3,2,6)
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
ylim([-0.8 0.8])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
xlabel('LME')
ylabel('Coeff')
title('Benthos')

print('-dpng',[ppath 'Bar_LMEs_climate_seasons_corrcoef_fntypes.png'])

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
figure(4)
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
    textm(TLAT(id(1)),TLONG(id(1)),str,...
        'HorizontalAlignment','center') 
%     textm(TLAT(id(2)),TLONG(id(2)),['lag=' num2str(LAtab(i,3))],...
%         'HorizontalAlignment','center') 
end
colormap(mcol)
title('All fishes')
print('-dpng',[ppath 'Map_LMEs_climate_seasons_corrcoef_allfish.png'])

%%
figure(5)
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
print('-dpng',[ppath 'Map_LMEs_climate_seasons_corrcoef_F.png'])

%%
figure(6)
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
%     textm(TLAT(id(1)),TLONG(id(1)),str,...
%         'HorizontalAlignment','center') 
%     textm(TLAT(id(2)),TLONG(id(2)),['lag=' num2str(LPtab(i,3))],...
%         'HorizontalAlignment','center') 
end
colormap(mcol)
title('Large pelagics')
print('-dpng',[ppath 'Map_LMEs_climate_seasons_corrcoef_P.png'])

%%
figure(7)
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
%     textm(TLAT(id(1)),TLONG(id(1)),str,...
%         'HorizontalAlignment','center') 
%     textm(TLAT(id(2)),TLONG(id(2)),['lag=' num2str(LDtab(i,3))],...
%         'HorizontalAlignment','center') 
end
colormap(mcol)
title('Demersals')
print('-dpng',[ppath 'Map_LMEs_climate_seasons_corrcoef_D.png'])

%%
figure(8)
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
%     textm(TLAT(id(1)),TLONG(id(1)),str,...
%         'HorizontalAlignment','center') 
%     textm(TLAT(id(2)),TLONG(id(2)),['lag=' num2str(LBtab(i,3))],...
%         'HorizontalAlignment','center') 
end
colormap(mcol)
title('Benthos')
print('-dpng',[ppath 'Map_LMEs_climate_seasons_corrcoef_B.png'])

%%
% on grid
Fcorr  = nan(ni,nj);
Pcorr  = nan(ni,nj);
Dcorr  = nan(ni,nj);
Acorr  = nan(ni,nj);
Bcorr  = nan(ni,nj);

%%
figure(9)
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
caxis([-0.7 0.7])
colorbar
title('All fishes')
print('-dpng',[ppath 'Map_LMEs_climate_seasons_corrcoef_allfish_v2.png'])
