% Plot max coeffs of driver-fish mult linear regressions
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

load([spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_maxcoefs.mat'])

%%  ---------------------------------------------------------
cnam = {'coef','p','idriver','driver'};
ctex = {'Det','TB','TP','Zmeso','ZmLoss','DetTB','ZmesoTP','ZmLossTP'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

% 1= AMO, 2=AO, 3=NAO, 4=Nino34, 5=PDO
mtype = {'o','s','^','v','d','p'}; %for lag
mcol = [
    34/255 136/255 51/255;...   %green
    170/255 51/255 119/255;...  %purple
    238/255 102/255 119/255;... %red
    0/255 68/255 136/255;...    %blue
    51/255 187/255 238/255;...  %cyan
    153/255 153/255 51/255;...  %olive
    0 0 0;...                   %black
    0.50 0.50 0.50;...          % grey
    ];


%colorblind friendly
% cb=[34/255 136/255 51/255;...   %green
%     153/255 153/255 51/255;...  %olive
%     51/255 187/255 238/255;...  %cyan
%     0/255 68/255 136/255;...    %blue
%     238/255 102/255 119/255;... %red
%     170/255 51/255 119/255;...  %purple
%     0 0 0;...                   %black
%     0.25 0.25 0.25;...          %dk grey
%     0.50 0.50 0.50;...          % grey
%     0.75 0.75 0.75];            %lt grey

%%
cid = [15,61,1,5,3];

figure(1)
% Get fake colors first for legend
for i=1:length(cid)
    L=cid(i);
    b=bar(L,LAtab(i,1),'EdgeColor','none','FaceColor',mcol(i,:));
    hold on
end
%legend of colors and shapes
lgd = legend(ctex,'Location','southwest');
lgd.AutoUpdate = 'off';
for i=1:length(lid)
    L=lid(i);
    if (LAtab(i,2) <= 0.05)
        b=bar(L,LAtab(i,1),'EdgeColor',mcol(LAtab(i,3),:),'FaceColor',mcol(LAtab(i,3),:));
        hold on
    else
        b=bar(L,LAtab(i,1),'EdgeColor',[1 1 1],'FaceColor',[1 1 1]);
        hold on
    end
    
end
xlim([0 67])
%ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
xlabel('LME')
ylabel('Coeff')
title('All fishes')
%lgd = legend(ctex,'Location','southeast');
print('-dpng',[ppath 'Bar_LMEs_driver_mlr_nu_maxcoef_allfish.png'])

%%
cid = [6,1,11,2,3,57,21,12];

f2 = figure('Units','inches','Position',[1 3 7.5 10]);
subplot('Position',[0.1 0.43 0.856 0.17])
% Get fake colors first for legend
for i=1:length(cid)
    L=cid(i);
    b=bar(L,LDtab(i,1),'EdgeColor','none','FaceColor',mcol(i,:));
    hold on
end
%legend of colors and shapes
lgd = legend(ctex,'Location','eastoutside');
lgd.AutoUpdate = 'off';

for i=1:length(lid)
    L=lid(i);
    if (LDtab(i,2) <= 0.05)
        b=bar(L,LDtab(i,1),'EdgeColor',mcol(LDtab(i,3),:),'FaceColor',mcol(LDtab(i,3),:));
        hold on
    else
        b=bar(L,LDtab(i,1),'EdgeColor',[1 1 1],'FaceColor',[1 1 1]);
        hold on
    end
    
end
xlim([0 67])
%ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('Demersals')

subplot('Position',[0.1 0.24 0.7 0.17])
for i=1:length(lid)
    L=lid(i);
    if (LAtab(i,2) <= 0.05)
        b=bar(L,LAtab(i,1),'EdgeColor',mcol(LAtab(i,3),:),'FaceColor',mcol(LAtab(i,3),:));
        hold on
    else
        b=bar(L,LAtab(i,1),'EdgeColor',[1 1 1],'FaceColor',[1 1 1]);
        hold on
    end
    
end
xlim([0 67])
%ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('All fishes')

subplot('Position',[0.1 0.81 0.7 0.17])
for i=1:length(lid)
    L=lid(i);
    if (LFtab(i,2) <= 0.05)
        b=bar(L,LFtab(i,1),'EdgeColor',mcol(LFtab(i,3),:),'FaceColor',mcol(LFtab(i,3),:));
        hold on
    else
        b=bar(L,LFtab(i,1),'EdgeColor',[1 1 1],'FaceColor',[1 1 1]);
        hold on
    end
    
end
xlim([0 67])
%ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('Forage fishes')

subplot('Position',[0.1 0.62 0.7 0.17])
for i=1:length(lid)
    L=lid(i);
    if (LPtab(i,2) <= 0.05)
        b=bar(L,LPtab(i,1),'EdgeColor',mcol(LPtab(i,3),:),'FaceColor',mcol(LPtab(i,3),:));
        hold on
    else
        b=bar(L,LPtab(i,1),'EdgeColor',[1 1 1],'FaceColor',[1 1 1]);
        hold on
    end
    
end
xlim([0 67])
%ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('Large pelagics')

print('-dpng',[ppath 'Bar_LMEs_driver_mlr_nu_maxcoef_fntypes.png'])

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
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.85 0.85 0.85]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Acorr(id) = LAtab(i,3);

    if (LAtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Acorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',LAtab(i,1))] , ['lag=' num2str(LAtab(i,3))]};
    loc = round(length(id)/2);
    %     textm(TLAT(id(loc)),TLONG(id(loc)),str,...
    %         'HorizontalAlignment','center')
end
colormap(mcol(1:5,:))
colorbar('Ticks',1:5,'TickLabels',ctex(1:5))
title('All fishes')
print('-dpng',[ppath 'Map_LMEs_driver_mlr_nu_maxcoef_allfish.png'])

%%
f1 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.10 0.32 0.4]) %D
axd=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.85 0.85 0.85]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Dcorr(id) = LDtab(i,3);

    if (LDtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Dcorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',LDtab(i,1))] , ['lag=' num2str(LDtab(i,3))]};

end
colormap(axd,mcol)
title('Demersals')
clb=colorbar('Position',[0.7 0.2 0.03 0.25],'Ticks',1:8,'TickLabels',ctex);
%clb.AutoUpdate = 'off';

subplot('Position',[0.01 0.575 0.32 0.4]) %F
axf=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.85 0.85 0.85]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Fcorr(id) = LFtab(i,3);

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
colormap(axf,mcol(1:7,:))
title('Forage fishes')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axp=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.85 0.85 0.85]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Pcorr(id) = LPtab(i,3);

    if (LPtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Pcorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',LPtab(i,1))] , ['lag=' num2str(LPtab(i,3))]};
end
colormap(axp,mcol)
title('Large pelagics')

subplot('Position',[0.33 0.10 0.32 0.4]) %A
axa=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.85 0.85 0.85]);
hold on
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    Acorr(id) = LAtab(i,3);

    if (LAtab(i,2) <= 0.05)
        surfm(TLAT,TLONG,Acorr)
        hold on
    end
    str = {['coef=' sprintf('%0.2f',LAtab(i,1))] , ['lag=' num2str(LAtab(i,3))]};
    loc = round(length(id)/2);
    %     textm(TLAT(id(loc)),TLONG(id(loc)),str,...
    %         'HorizontalAlignment','center')
end
colormap(axa,mcol(1:5,:))
title('All fishes')

print('-dpng',[ppath 'Map_LMEs_driver_mlr_nu_maxcoef_fntypes.png'])



