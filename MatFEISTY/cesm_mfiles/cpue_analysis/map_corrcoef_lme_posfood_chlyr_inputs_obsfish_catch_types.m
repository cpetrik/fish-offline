% Mapp corr coefs of ind driver combos with catch
% Subplots together for comparison
% Restricted analysis to chl yrs
% Obs effort

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
load([spath 'LMEs_corr_catch_chlyrs_inputs_obsfish_mostsiglag_posfood.mat'])
%dim: LME x driver x (corr, p-val, lag)

%% Map data
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
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

cmR = cbrewer('seq','Reds',10,'PCHIP');
cmRB = cbrewer('div','RdBu',21,'PCHIP');
cmRB = flipud(cmRB);

subpos = [0.015 0.75 0.43 0.25;...
    0.015 0.5 0.43 0.25;...
    0.015 0.25 0.43 0.25;...
    0.015 0.0 0.43 0.25;...
    0.48 0.75 0.43 0.25;...
    0.48 0.5 0.43 0.25;...
    0.48 0.25 0.43 0.25;...
    0.48 0.0 0.43 0.25];

%% All fish 
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
for j=1:8 
    Cmat = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (LAtab(i,j,2) <= 0.05)
            Cmat(id) = LAtab(i,j,1);
        end
    end

    subplot('Position',subpos(j,:))
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(TLAT,TLONG,Cmat)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 1])
    set(gcf,'renderer','painters')
    text(0,1.75,tanom{j},'HorizontalAlignment','center')

end
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
print('-dpng',[ppath 'Map_LMEs_chlyr_catch_driver_obsfish_corrcoef_All.png'])

%% Forage
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
for j=1:8 
    Cmat = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (LFtab(i,j,2) <= 0.05)
            Cmat(id) = LFtab(i,j,1);
        end
    end

    subplot('Position',subpos(j,:))
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(TLAT,TLONG,Cmat)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 1])
    set(gcf,'renderer','painters')
    text(0,1.75,tanom{j},'HorizontalAlignment','center')

end
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
print('-dpng',[ppath 'Map_LMEs_chlyr_catch_driver_obsfish_corrcoef_F.png'])

%% Lg Pel
f3 = figure('Units','inches','Position',[1 3 6.5 8]);
for j=1:8 
    Cmat = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (LPtab(i,j,2) <= 0.05)
            Cmat(id) = LPtab(i,j,1);
        end
    end

    subplot('Position',subpos(j,:))
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(TLAT,TLONG,Cmat)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 1])
    set(gcf,'renderer','painters')
    text(0,1.75,tanom{j},'HorizontalAlignment','center')

end
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
print('-dpng',[ppath 'Map_LMEs_chlyr_catch_driver_obsfish_corrcoef_P.png'])

%% Demersal
f4 = figure('Units','inches','Position',[1 3 6.5 8]);
for j=1:8 
    Cmat = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (LDtab(i,j,2) <= 0.05)
            Cmat(id) = LDtab(i,j,1);
        end
    end

    subplot('Position',subpos(j,:))
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(TLAT,TLONG,Cmat)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 1])
    set(gcf,'renderer','painters')
    text(0,1.75,tanom{j},'HorizontalAlignment','center')

end
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
print('-dpng',[ppath 'Map_LMEs_chlyr_catch_driver_obsfish_corrcoef_D.png'])

%% 4-plot of catch corrs
Amat = nan(ni,nj);
Fmat = nan(ni,nj);
Pmat = nan(ni,nj);
Dmat = nan(ni,nj);

for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);

    if (LAtab(i,9,2) <= 0.05)
        Amat(id) = LAtab(i,9,1);
    end
    if (LFtab(i,9,2) <= 0.05)
        Fmat(id) = LFtab(i,9,1);
    end
    if (LPtab(i,9,2) <= 0.05)
        Pmat(id) = LPtab(i,9,1);
    end
    if (LDtab(i,9,2) <= 0.05)
        Dmat(id) = LDtab(i,9,1);
    end
end

f5 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.05 0.555 0.45 0.4]) %Top L
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Fmat)
colormap(cmR);
clim([0 1])
title('Forage')

subplot('Position',[0.5 0.555 0.45 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Pmat)
colormap(cmR);
clim([0 1])
title('Lg Pelagics')

subplot('Position',[0.05 0.105 0.45 0.4]) %Bottom L
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Dmat)
colormap(cmR);
clim([0 1])
title('Demersals')

subplot('Position',[0.5 0.105 0.45 0.4]) %Bottom R
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,Amat)
colormap(cmR);
clim([0 1])
colorbar('Position',[0.25 0.055 0.5 0.03],'orientation','horizontal')
title('All fishes')
print('-dpng',[ppath 'Map_LMEs_chlyr_catch_driver_obsfish_corrcoef_yield.png'])

%% 8plot by driver
% f1 = figure('Units','inches','Position',[1 3 6.5 8]);
% %f1.Units = 'inches';
% 
% %A - F biom
% subplot('Position',[0.015 0.75 0.43 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(TLAT2,TLON2,log10(AllF2))
% colormap(cmBP50)
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-1.5 2.0])
% set(gcf,'renderer','painters')
% text(0,1.75,'F biom','HorizontalAlignment','center')
% 
% %B - P biom
% subplot('Position',[0.015 0.5 0.43 0.25])
% 
% %C - D biom
% subplot('Position',[0.015 0.25 0.43 0.25])
% 
% %D - A biom
% subplot('Position',[0.015 0.0 0.43 0.25])
% 
% %E - F prod
% subplot('Position',[0.51 0.75 0.43 0.25])
% 
% %F - P prod
% subplot('Position',[0.51 0.5 0.43 0.25])
% 
% %G - D prod
% subplot('Position',[0.51 0.25 0.43 0.25])
% 
% %H - All prod
% subplot('Position',[0.51 0.0 0.43 0.25])