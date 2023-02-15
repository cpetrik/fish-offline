% CESM FOSI output

clear 
close all

%% Paths

fpath='/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6_noSeas.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');

%% FEISTY Inputs
% load means
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FEISTY_forcing_means_correct_units.mat'],...
    'TLAT','TLONG','cTp','cTb','cD','cZ','cZl');

%% fix seam
[TLAT2,TLON2,Tp] = cyclic_map_seam_cesm(TLAT,TLONG,cTp);
[~,~,Tb]  = cyclic_map_seam_cesm(TLAT,TLONG,cTb);
[~,~,Det] = cyclic_map_seam_cesm(TLAT,TLONG,cD);
[~,~,Zm]  = cyclic_map_seam_cesm(TLAT,TLONG,cZ);
[~,~,ZmL] = cyclic_map_seam_cesm(TLAT,TLONG,cZl);

%%
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

%%
f1 = figure('Units','inches','Position',[1 3 6.5 3]);

subplot('Position',[0.01 0.575 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,Tp)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.0451 0.565 0.2455 0.03],'orientation','horizontal')
text(0.2,1.65,'TP','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.33 0.575 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,log10(Zm))
cmocean('tempo')
caxis([0 2])
colorbar('Position',[0.368 0.565 0.245 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zmeso','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.65 0.575 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,log10(ZmL))
cmocean('tempo')
caxis([-1 1])
colorbar('Position',[0.6855 0.565 0.2455 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 ZmLoss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.10 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,Tb)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.0451 0.09 0.2455 0.03],'orientation','horizontal')
text(0.2,1.65,'TB','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.33 0.10 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT2,TLON2,log10(Det))
cmocean('tempo')
caxis([-2.5 1.5])
colorbar('Position',[0.368 0.09 0.245 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Det','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CESM_FOSI_mean_forcings.png'])

