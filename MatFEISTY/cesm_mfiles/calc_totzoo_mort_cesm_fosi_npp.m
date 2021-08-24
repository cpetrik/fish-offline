% CESM FOSI output
% Split zoo into Large zoo and zoo loss
% Use NPP instead of biomass to calc frac large

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6.mat'],'GRD');

%% Plankton inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'zooC_150m','zoo_loss_150m',...
    'FillValue','missing_value','TLAT','TLONG','TAREA','time');

load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.phytoNPP_FIESTY-forcing.mat'],...
    'yr','photoC_diaz_150m','photoC_diat_150m','photoC_sp_150m');

%%
% zooC_150m
zooC_150m_units = 'mmol/m^3 cm';
zooC_150m_name  = 'zoo biomass';

% zoo_loss_150m
zoo_loss_150m_units = 'mmol/m^3/s cm';
zoo_loss_150m_name  = 'zoo loss';

% photoC_diaz_150m
photoC_diaz_150m_units = 'mmol/m^3/s cm';
photoC_diaz_150m_name  = 'thickness of layer k';

% photoC_diat_150m
photoC_diat_150m_units = 'mmol/m^3/s cm';
photoC_diat_150m_name  = 'thickness of layer k';

% photoC_sp_150m
photoC_sp_150m_units = 'mmol/m^3/s cm';
photoC_sp_150m_name  = 'thickness of layer k';

%% doubles
zooC_150m = double(zooC_150m); 
zoo_loss_150m = double(zoo_loss_150m); 
photoC_diaz_150m = double(photoC_diaz_150m); 
photoC_diat_150m = double(photoC_diat_150m); 
photoC_sp_150m = double(photoC_sp_150m);

%% fraction large phyto -> frac large zoo
fracL = photoC_diat_150m ./ (photoC_diat_150m + photoC_diaz_150m + photoC_sp_150m + eps);

LzooC_150m = fracL .* zooC_150m;
Lzoo_loss_150m = fracL .* zoo_loss_150m;

%%
save([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo_v5.mat'],...
    'fracL','LzooC_150m','Lzoo_loss_150m','zooC_150m_units',...
    'zoo_loss_150m_units');

%% Units
% meso zoo: nmolC cm-2 to g(WW) m-2
zooC_150m = zooC_150m * 1e-9 * 1e4 * 12.01 * 9.0;
LzooC_150m = LzooC_150m * 1e-9 * 1e4 * 12.01 * 9.0;

% meso zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
zoo_loss_150m = zoo_loss_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
Lzoo_loss_150m = Lzoo_loss_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

%%
cFL = mean(fracL,3);
cZ  = mean(LzooC_150m,3);
cZl = mean(Lzoo_loss_150m,3);
cT  = mean(zooC_150m,3);
cTl = mean(zoo_loss_150m,3);

cZ(cZ<=0) = eps;
cZl(cZl<=0) = eps;
cT(cT<=0) = eps;
cTl(cTl<=0) = eps;

%%
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

figure(2)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cT))
cmocean('tempo')
caxis([0 2])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 TotZoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cTl))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 TotZoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%subplot('Position',[0.01 0.37 0.4 0.3])

%subplot('Position',[0.41 0.37 0.4 0.3])
subplot('Position',[0.21 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(cFL))
cmocean('balance')
caxis([0 1])
colorbar%('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Frac Lg','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZ))
cmocean('tempo')
caxis([0 2])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 LgZoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZl))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 LgZoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CESM_FOSI_mean_zoo_loss_v5_Psplit.png'])

