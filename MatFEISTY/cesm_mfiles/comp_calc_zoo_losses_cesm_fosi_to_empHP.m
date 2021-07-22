% CESM FOSI output
% Back calc zoo quad mort = loss - linear mort

clear all
close all

%% Paths

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6.mat'],'GRD');

%% FEISTY Inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units',...
    'TLAT','TLONG','TAREA','time','yr');
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo.mat'],...
    'LzooC_150m','Lzoo_loss_150m');

%% nans & zeros
TEMP_150m = double(TEMP_150m);
LzooC_150m(LzooC_150m<0) = 0.0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0.0;

%% Units
%zoo loss: mmol/m^3/s cm
%zoo: mmolC/m^3 cm
%tp: degC

% meso zoo: nmolC cm-2 to g(WW) m-2
LzooC_150m = LzooC_150m * 1e-9 * 1e4 * 12.01 * 9.0;

% meso zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
Lzoo_loss_150m = Lzoo_loss_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

%% Temp fn
T = min(TEMP_150m(:)):max(TEMP_150m(:));
tfn = exp(-4000 .* ( (1./(T+273.15)) - (1./303.15) ));

figure(1)
plot(T,tfn);

Tfn = exp(-4000 .* ( (1./(TEMP_150m+273.15)) - (1./303.15) ));

%% Zprime
Zprime = max((LzooC_150m - 0.01),0);

%% Z quad
Lzoo_quad_150m = Lzoo_loss_150m - (Tfn .* 0.1 .* Zprime);
Lzoo_quad_150m = max(Lzoo_quad_150m,0);

%%
% HP loss is an empirical fitted fn of biomass and temp
dZm = 10 .^ (-2.925 + 1.964.*log10(LzooC_150m+eps) + 1.958e-2.*TEMP_150m);
dZm = max(dZm,0);

%%
cTp = mean(TEMP_150m,3);
cZ = mean(LzooC_150m,3);
cZl = mean(Lzoo_loss_150m,3);
cZq = mean(Lzoo_quad_150m,3);
cHP = mean(dZm,3);

%% theoretical
testHP = 10 .^ (-2.925 + 1.964.*log10(nanmean(LzooC_150m(:))) + 1.958e-2.*T);
testQ = nanmean(Lzoo_loss_150m(:)) - (tfn .* 0.1 .* nanmean(Zprime(:)));

Z = 1:0.5:9;
zooHP = 10 .^ (-2.925 + 1.964.*log10(Z) + 1.958e-2.*nanmean(TEMP_150m(:)));
zooQ = nanmean(Lzoo_loss_150m(:)) - (nanmean(TEMP_150m(:)) .* 0.1 .* Z);

figure(2)
subplot(1,2,1)
plot(T,testHP,'b','LineWidth',2); hold on;
plot(T,testQ,'k','LineWidth',2)
xlabel('Temp')
ylabel('Loss')

subplot(1,2,2)
plot(Z,zooHP,'b','LineWidth',2); hold on;
plot(Z,zooQ,'k','LineWidth',2)
xlabel('ZooC')
ylabel('Loss')
print('-dpng',[pp 'Theoretical_CESM_FOSI_COBALTempHP_zoo_loss_quad.png'])

%%
figure(3)
plot(cTp(:),cZq(:),'k.'); hold on;
plot(cTp(:),cHP(:),'b.'); hold on;
ylim([0 1])

%%
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

figure(3)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cHP))
cmocean('tempo')
caxis([-2 1])
colorbar%('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 empHP Zoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZ))
cmocean('tempo')
caxis([0 2])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZl))
cmocean('tempo')
caxis([-2 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo tot loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZq))
cmocean('tempo')
caxis([-2 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo quad loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%subplot('Position',[0.01 0.06 0.4 0.3])
print('-dpng',[pp 'Map_CESM_FOSI_mean_zoo_loss_quad_empHP.png'])

