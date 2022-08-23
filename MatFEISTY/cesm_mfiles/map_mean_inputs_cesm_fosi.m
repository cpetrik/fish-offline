% CESM FOSI output

clear all
close all

%% Paths

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6_noSeas.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');

%% FEISTY Inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'TLAT','TLONG','TAREA','time','yr');
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo_totloss_allphytoC.mat']);

%% nans & zeros
TEMP_150m = double(TEMP_150m);
TEMP_bottom = double(TEMP_bottom);
POC_FLUX_IN_bottom = double(POC_FLUX_IN_bottom);

TEMP_bottom(TEMP_bottom >= 9.9e+36) = nan;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;

LzooC_150m(LzooC_150m<0) = 0.0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0.0;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom<0) = 0.0;

%% Units
%poc flux: mmol/m^3 cm/s
%zoo loss: mmol/m^3/s cm
%zoo: mmolC/m^3 cm
%tp: degC
%tb: degC

% meso zoo: nmolC cm-2 to g(WW) m-2
% 1e9 nmol in 1 mol C
% 1e4 cm2 in 1 m2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
LzooC_150m = LzooC_150m * 1e-9 * 1e4 * 12.01 * 9.0;

% meso zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
% detrital flux to benthos: nmolC cm-2 s-1 to g(WW) m-2 d-1
% 1e9 nmol in 1 mol C
% 1e4 cm2 in 1 m2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
Lzoo_loss_150m = Lzoo_loss_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
%Lzoo_quad_150m = Lzoo_quad_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
POC_FLUX_IN_bottom = POC_FLUX_IN_bottom * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

%quad already in correct units

%%
cTp = mean(TEMP_150m,3);
cTb = mean(TEMP_bottom,3);
cD = mean(POC_FLUX_IN_bottom,3);
cZ = mean(LzooC_150m,3);
cZl = mean(Lzoo_loss_150m,3);
%cZq = mean(Lzoo_quad_150m,3);

%% save means
save([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FEISTY_forcing_means_correct_units.mat'],...
    'TLAT','TLONG','cTp','cTb','cD','cZ','cZl');

%%
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

%%
figure(1)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cTp)
cmocean('thermal')
caxis([0 35])
colorbar%('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Tp','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cTb)
cmocean('thermal')
caxis([0 35])
colorbar%('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Tb','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZ))
cmocean('tempo')
caxis([0 2])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 MesoZ','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cD))
cmocean('tempo')
caxis([-2.5 1.5])
colorbar%('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Det','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZl))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 MesoZ loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CESM_FOSI_mean_forcings.png'])

%%
figure(2)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cTp)
cmocean('thermal')
caxis([0 35])
colorbar%('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Tp','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cTb)
cmocean('thermal')
caxis([0 35])
colorbar%('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Tb','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZ))
cmocean('tempo')
caxis([0 2])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cD))
cmocean('tempo')
caxis([-2.5 1.5])
colorbar%('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Det','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% subplot('Position',[0.01 0.06 0.4 0.3])
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(TLAT,TLONG,log10(cZq))
% cmocean('tempo')
% caxis([-1 1])
% colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
% text(0.2,1.65,'log_1_0 Zoo quad loss','HorizontalAlignment','center','FontWeight','bold')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% 
% print('-dpng',[pp 'Map_CESM_FOSI_mean_forcings_quad.png'])

%%
figure(3)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cTp)
cmocean('thermal')
caxis([0 35])
colorbar%('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Tp','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cTb)
cmocean('thermal')
caxis([0 35])
colorbar%('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Tb','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZ))
cmocean('tempo')
caxis([0 2])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cD))
cmocean('tempo')
caxis([-2.5 1.5])
colorbar%('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Det','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZl))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% subplot('Position',[0.41 0.06 0.4 0.3])
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(TLAT,TLONG,log10(cZq))
% cmocean('tempo')
% caxis([-1 1])
% colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
% text(0.2,1.65,'log_1_0 Zoo quad loss','HorizontalAlignment','center','FontWeight','bold')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% print('-dpng',[pp 'Map_CESM_FOSI_mean_forcings_all.png'])

%% Map where zloss = 0
nzLoss = cZl;
nzLoss(cZl>0) = 1;

% nzQuad = cZq;
% nzQuad(cZq>0) = 1;

%%
figure(4)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,nzLoss)
%cmocean('tempo')
caxis([0 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% subplot(2,1,2)
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(TLAT,TLONG,nzQuad)
% %cmocean('tempo')
% caxis([0 1])
% colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
% text(0.2,1.65,'log_1_0 Zoo quad loss','HorizontalAlignment','center','FontWeight','bold')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% print('-dpng',[pp 'Map_CESM_FOSI_mean_losses_zeros.png'])

%%
figure(5)
subplot(2,2,1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cZ)
cmocean('balance')
caxis([-0.1 0.1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Zoo biom','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot(2,2,2)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cZl)
cmocean('balance')
caxis([-0.1 0.1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Zoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% subplot(2,2,4)
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(TLAT,TLONG,cZq)
% cmocean('balance')
% caxis([-0.1 0.1])
% colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
% text(0.2,1.65,'Zoo quad loss','HorizontalAlignment','center','FontWeight','bold')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% %print('-dpng',[pp 'Map_CESM_FOSI_mean_zoop_negs.png'])


