% CESM FOSI output
% map mean & interann variability by grid cell 
% no Zmeso
% v15 of Lfrac calc & ZmLoss calc

clear 
close all

%% Paths
%fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
fpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6_noSeas.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');
load([fpath 'LME-mask-POP_gx1v6.mat']);

%% FEISTY Inputs
fpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([fpath 'CESM_FOSI_v15_interann_var_forcings.mat'],...
    'tp_mean','tb_mean','dety_mean','mzly_mean',...
    'cvtp','cvtb','cvdety','cvzlosy');

%% fix seam
[TLAT2,TLON2,mTp] = cyclic_map_seam_cesm(TLAT,TLONG,tp_mean);
[~,~,mTb]  = cyclic_map_seam_cesm(TLAT,TLONG,tb_mean);
[~,~,mDet] = cyclic_map_seam_cesm(TLAT,TLONG,dety_mean);
[~,~,mZmL] = cyclic_map_seam_cesm(TLAT,TLONG,mzly_mean);

[~,~,vTp]  = cyclic_map_seam_cesm(TLAT,TLONG,cvtp);
[~,~,vTb]  = cyclic_map_seam_cesm(TLAT,TLONG,cvtb);
[~,~,vDet] = cyclic_map_seam_cesm(TLAT,TLONG,cvdety);
[~,~,vZmL] = cyclic_map_seam_cesm(TLAT,TLONG,cvzlosy);

%% map info
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

cmYOR=cbrewer('seq','YlOrRd',66,'pchip');
cmYOB=cbrewer('seq','YlOrBr',66,'pchip');
cmOR=cbrewer('seq','OrRd',66,'pchip');

%% map of CV
f1 = figure('Units','inches','Position',[1 5 7.5 5.5]);

subplot('Position',[0.01 0.575 0.45 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,vTp)
colormap(cmYOR)
caxis([0 0.5])
%colorbar('Position',[0.0451 0.565 0.2455 0.03],'orientation','horizontal')
text(0.2,1.65,'TP','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.575 0.45 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,(vZmL))
colormap(cmYOR)
caxis([0 0.5])
%colorbar('Position',[0.6855 0.565 0.2455 0.03],'orientation','horizontal')
text(0.2,1.65,'ZmLoss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.125 0.45 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,vTb)
colormap(cmYOR)
caxis([0 0.5])
%colorbar('Position',[0.0451 0.09 0.2455 0.03],'orientation','horizontal')
text(0.2,1.65,'TB','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.125 0.45 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,(vDet))
colormap(cmYOR)
caxis([0 0.5])
colorbar('Position',[0.25 0.075 0.5 0.03],'orientation','horizontal')
text(0.2,1.65,'Det','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%subplot('Position',[0.65 0.575 0.32 0.4])

print('-dpng',[pp 'Map_CESM_FOSI_v15_interann_coeffvar_forcings_ms.png'])


%% map of biomass/flux
f2 = figure('Units','inches','Position',[1 5 7.5 5.5]);

subplot('Position',[0.01 0.575 0.45 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,mTp)
cmocean('thermal')
caxis([0 35])
%colorbar('Position',[0.0451 0.565 0.2455 0.03],'orientation','horizontal')
text(0.2,1.65,'TP','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.575 0.45 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,log10(mZmL))
cmocean('tempo')
caxis([0 4])
%colorbar('Position',[0.55 0.565 0.35 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 ZmLoss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.125 0.45 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,mTb)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.075 0.35 0.03],'orientation','horizontal')
text(0.2,1.65,'TB','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.125 0.45 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,log10(mDet))
cmocean('tempo')
caxis([0 4])
colorbar('Position',[0.55 0.075 0.35 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Det','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CESM_FOSI_v15_mean_forcings_ms.png'])

