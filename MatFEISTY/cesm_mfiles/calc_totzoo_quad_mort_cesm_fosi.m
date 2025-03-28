% CESM FOSI output
% Calc zoo quad mort
% Then split into Large zoo and zoo loss

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6.mat'],'GRD');

%% Plankton inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'zooC_150m','zoo_loss_150m','diatC_150m','spC_150m',...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units',...
    'TLAT','TLONG','TAREA','time','yr');

%%
% Size:       320x384x816
% Dimensions: nlon,nlat,time
FillValue   = NaN;
missing_value = 9.969209968386869e+36;

% zooC_150m
zooC_150m_units = 'mmol/m^3 cm';
zooC_150m_name  = 'zoo biomass';

% zoo_loss_150m
zoo_loss_150m_units = 'mmol/m^3/s cm';
zoo_loss_150m_name  = 'zoo loss';

% diatC_150m
diatC_150m_units = 'mmol/m^3 cm';
diatC_150m_name  = 'diatom biomass';

% spC_150m
spC_150m_units = 'mmol/m^3 cm';
spC_150m_name  = 'SP biomass';

%% doubles, nans & zeros
zooC_150m = double(zooC_150m); 
zoo_loss_150m = double(zoo_loss_150m); 
diatC_150m = double(diatC_150m); 
spC_150m = double(spC_150m);  
TEMP_150m = double(TEMP_150m);

zooC_150m(isnan(TEMP_150m)) = nan; 
zoo_loss_150m(isnan(TEMP_150m)) = nan; 
diatC_150m(isnan(TEMP_150m)) = nan; 
spC_150m(isnan(TEMP_150m)) = nan; 

%% Back-calc quad mort
% Units - NEED TO KEEP IN NATIVE UNITS FOR CONSTANTS TO APPLY

% Temp fn
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

Zprime = max((zooC_150m - C_loss_thres),0.0);

zid = find(Zprime(:)==0);
length(zid)/(320*384*816) %0.2998 = 30% of time/space

zid2 = find(zoo_loss_150m(:)==0);
length(zid2)/(320*384*816) %0.0112 = 1% of time/space

%% Z quad
spd = 86400; %seconds per day
dps = 1 ./ spd; %days per second
parm_z_mort = 0.08 * dps;
parm_z_mort2 = 0.42 * dps;

zmort = parm_z_mort .* Tfn;
zmort2 = parm_z_mort2 .* Tfn;
%zoo_loss = (zmort2 .* Zprime.^1.4) + (zmort .* Zprime);
zoo_lin_150m = zmort .* Zprime;
zoo1_quad_150m = zoo_loss_150m - (zmort .* Zprime);
zoo2_quad_150m = (zmort2 .* Zprime.^1.4);

nid = find(zoo2_quad_150m(:)<0);
length(nid)/(320*384*816) %0 = 0% of time/space

nid2 = find(zoo2_quad_150m(:)==0);
length(nid2)/(320*384*816) %0.0014 = 1.4% of time/space

zoo_lin_150m(zoo_lin_150m<0) = 0;
zoo1_quad_150m(zoo1_quad_150m<0) = 0;
zoo2_quad_150m(zoo2_quad_150m<0) = 0;

%% fraction large phyto -> frac large zoo
fracL = diatC_150m ./ (diatC_150m + spC_150m + eps);

LzooC_150m = fracL .* zooC_150m;
Lzoo_loss_150m = fracL .* zoo_loss_150m;
Lzoo_lin_150m = fracL .* zoo_lin_150m;
Lzoo1_quad_150m = fracL .* zoo1_quad_150m;
Lzoo2_quad_150m = fracL .* zoo2_quad_150m;

%%
save([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo_v4.mat'],...
    'fracL','LzooC_150m','Lzoo_loss_150m','zooC_150m_units',...
    'zoo_loss_150m_units','Lzoo_quad_150m');

%% Units
% meso zoo: nmolC cm-2 to g(WW) m-2
zooC_150m = zooC_150m * 1e-9 * 1e4 * 12.01 * 9.0;
LzooC_150m = LzooC_150m * 1e-9 * 1e4 * 12.01 * 9.0;

% meso zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
zoo_loss_150m = zoo_loss_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
Lzoo_loss_150m = Lzoo_loss_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
Lzoo_lin_150m = Lzoo_lin_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
Lzoo1_quad_150m = Lzoo1_quad_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
Lzoo2_quad_150m = Lzoo2_quad_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

%%
cTp = mean(TEMP_150m,3);
cTf = mean(Tfn,3);
cFL = mean(fracL,3);
cZ  = mean(LzooC_150m,3);
cZt = mean(Lzoo_loss_150m,3);
cZl = mean(Lzoo_lin_150m,3);
cZq1 = mean(Lzoo1_quad_150m,3);
cZq2 = mean(Lzoo2_quad_150m,3);
cT  = mean(zooC_150m,3);
cTt = mean(zoo_loss_150m,3);

%% 
cZ(cZ<=0) = eps;
cT(cT<=0) = eps;
cTt(cTt<=0) = eps;
cZt(cZt<=0) = eps;
cZl(cZl<=0) = eps;
cZq1(cZq1<=0) = eps;
cZq2(cZq2<=0) = eps;

%%
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

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
surfm(TLAT,TLONG,cTf)
cmocean('balance')
caxis([0 1])
colorbar%('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Tfn','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZ))
cmocean('tempo')
caxis([0 2])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 LgZoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
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
surfm(TLAT,TLONG,log10(cZt))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 LgZoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZq2))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 LgZoo quad loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CESM_FOSI_mean_zoo_loss_quad_v4_calc_before_Bsplit.png'])

%%
figure(2)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZ))
cmocean('tempo')
caxis([0 2])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 LgZoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(cFL))
cmocean('balance')
caxis([0 1])
colorbar%('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Frac Lg','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZt))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 LgZoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZl))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 LgZoo lin loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZq1))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 LgZoo quad backcalc','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZq2))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 LgZoo quad calc','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CESM_FOSI_mean_zoo_loss_all_v4_calc_before_Bsplit.png'])

%%
figure(3)
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
surfm(TLAT,TLONG,log10(cTt))
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
surfm(TLAT,TLONG,log10(cZt))
cmocean('tempo')
caxis([-1 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 LgZoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CESM_FOSI_mean_zoo_loss_v2_Bsplit.png'])
