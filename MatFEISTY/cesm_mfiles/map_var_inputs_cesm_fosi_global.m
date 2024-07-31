% CESM FOSI output
% calc interann variability by grid cell 

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
    'cvtp','cvtb','cvdet','cvzoo','cvzlos');

%% fix seam
[TLAT2,TLON2,Tp] = cyclic_map_seam_cesm(TLAT,TLONG,cvtp);
[~,~,Tb]  = cyclic_map_seam_cesm(TLAT,TLONG,cvtb);
[~,~,Det] = cyclic_map_seam_cesm(TLAT,TLONG,cvdet);
[~,~,Zm]  = cyclic_map_seam_cesm(TLAT,TLONG,cvzoo);
[~,~,ZmL] = cyclic_map_seam_cesm(TLAT,TLONG,cvzlos);

%% map info
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

cmYOR=cbrewer('seq','YlOrRd',66,'pchip');
cmYOB=cbrewer('seq','YlOrBr',66,'pchip');
cmOR=cbrewer('seq','OrRd',66,'pchip');

%% map
f1 = figure('Units','inches','Position',[1 3 6.5 3]);

subplot('Position',[0.01 0.575 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,Tp)
colormap(cmYOR)
caxis([0 0.5])
%colorbar('Position',[0.0451 0.565 0.2455 0.03],'orientation','horizontal')
text(0.2,1.65,'TP','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.33 0.575 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,(Zm))
colormap(cmYOR)
caxis([0 0.5])
%colorbar('Position',[0.368 0.565 0.245 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zmeso','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.65 0.575 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,(ZmL))
colormap(cmYOR)
caxis([0 0.5])
%colorbar('Position',[0.6855 0.565 0.2455 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 ZmLoss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.125 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,Tb)
colormap(cmYOR)
caxis([0 0.5])
%colorbar('Position',[0.0451 0.09 0.2455 0.03],'orientation','horizontal')
text(0.2,1.65,'TB','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.33 0.125 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,(Det))
colormap(cmYOR)
caxis([0 0.5])
colorbar('Position',[0.3 0.095 0.37 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Det','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CESM_FOSI_v15_interann_coeffvar_forcings.png'])

