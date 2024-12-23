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

% Units - NEED TO KEEP IN NATIVE UNITS FOR CONSTANTS TO APPLY

%% Temp fn
T = min(TEMP_150m(:)):max(TEMP_150m(:));
Q10 = 2.0;
T0Kelv = 273.15;
Tref = 30.0;
tfn = Q10 .^ ( ( (T + T0Kelv) - (Tref + T0Kelv) ) / 10.0 );

% figure(1)
% plot(T,tfn);

Tfn = Q10 .^ ( ( (TEMP_150m + T0Kelv) - (Tref + T0Kelv) ) ./ 10.0 );

%% Zprime
avg_f_loss_thres = 0.9167;
loss_thres_zoo = 0.2; 
C_loss_thres = avg_f_loss_thres * loss_thres_zoo;

Zprime = max((LzooC_150m - C_loss_thres),0.0);

zid = find(Zprime(:)==0);
length(zid)/(320*384*816) %0.2998 = 30% of time/space

zid2 = find(Lzoo_loss_150m(:)==0);
length(zid2)/(320*384*816) %0.0112 = 1.1% of time/space

%% Z quad
spd = 86400; %seconds per day
dps = 1 ./ spd; %days per second
parm_z_mort = 0.08 * dps;
parm_z_mort2 = 0.42 * dps;

zmort = parm_z_mort .* Tfn;
zmort2 = parm_z_mort2 .* Tfn;
%zoo_loss = (zmort2 .* Zprime.^1.4) + (zmort .* Zprime);
Lzoo_quad_150m = Lzoo_loss_150m - (zmort .* Zprime);

nid = find(Lzoo_quad_150m(:)<0);
length(nid)/(320*384*816) %0.0717 = 7.2% of time/space

nid2 = find(Lzoo_quad_150m(:)==0);
length(nid2)/(320*384*816) %0.0014 = 1.4% of time/space

Lzoo_quad_150m = max(Lzoo_quad_150m,0);

%%
save([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo.mat'],...
    'Lzoo_quad_150m','-append');

%%
cTp = mean(TEMP_150m,3);
cZ = mean(LzooC_150m,3);
cZl = mean(Lzoo_loss_150m,3);
cZq = mean(Lzoo_quad_150m,3);

%% Convert Units to plot
%zoo loss: mmol/m^3/s cm
%zoo: mmolC/m^3 cm
%tp: degC

% meso zoo: nmolC cm-2 to g(WW) m-2
cZ = cZ * 1e-9 * 1e4 * 12.01 * 9.0;

% meso zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
cZl = cZl * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
cZq = cZq * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

%%
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

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
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo tot loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZq))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo quad loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%subplot('Position',[0.01 0.06 0.4 0.3])
print('-dpng',[pp 'Map_CESM_FOSI_mean_zoo_loss_quad_v2.png'])

