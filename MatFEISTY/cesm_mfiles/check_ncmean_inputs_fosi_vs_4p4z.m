% check inputs, because CORE zoo biomass too low

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/4P4Z/';
% FOSI
fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([fpath 'gridspec_POP_gx1v6.mat']);
load([fpath 'Data_grid_POP_gx1v6.mat'],'GRD');
% 4P4Z
zpath='/Volumes/MIP/GCM_DATA/CESM/4P4Z/';

%% FOSI
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'time','yr');
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo.mat'],...
    'LzooC_150m','Lzoo_loss_150m');

%% nans
TEMP_bottom(TEMP_bottom >= 9.9e+36) = nan;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;

%zeros
LzooC_150m(LzooC_150m<0) = 0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0;

%%
tp = double(TEMP_150m);
tb = double(TEMP_bottom);
det_btm = double(POC_FLUX_IN_bottom);
mz = double(LzooC_150m);
hploss_mz = double(Lzoo_loss_150m);

% means in space
tp_mean_fosi = mean(tp,3);
tb_mean_fosi = mean(tb,3);
det_mean_fosi = mean(det_btm,3);
mz_mean_fosi = mean(mz,3);
mzloss_mean_fosi = mean(hploss_mz,3);
mzloss_min_fosi = min(hploss_mz,[],3);

Ftime = time;
Fyr = yr;

clear TEMP_150m TEMP_bottom POC_FLUX_IN_bottom LzooC_150m Lzoo_loss_150m
clear tp tb det_btm mz hploss_mz time yr

%% Compare to 4P4Z 
load([zpath 'g.e22a06.G1850ECOIAF_JRA_PHYS_DEV.TL319_g17.4p4z.004.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'time','yr','zoo3C_150m','zoo3_loss_150m',...
    'zoo4C_150m','zoo4_loss_150m');

% doubles
TEMP_150m = double(TEMP_150m);
TEMP_bottom = double(TEMP_bottom);
POC_FLUX_IN_bottom = double(POC_FLUX_IN_bottom);
zoo3C_150m = double(zoo3C_150m);
zoo4C_150m = double(zoo4C_150m);
zoo3_loss_150m = double(zoo3_loss_150m);
zoo4_loss_150m = double(zoo4_loss_150m);

%zeros
zoo3C_150m(zoo3C_150m<0) = 0;
zoo4C_150m(zoo4C_150m<0) = 0;
zoo3_loss_150m(zoo3_loss_150m<0) = 0;
zoo4_loss_150m(zoo4_loss_150m<0) = 0;

zooC_150m = zoo3C_150m + zoo4C_150m;
zoo_loss_150m = zoo3_loss_150m + zoo4_loss_150m;

%%
tp_mean_4p4z = mean(TEMP_150m,3);
tb_mean_4p4z = mean(TEMP_bottom,3);
det_mean_4p4z = mean(POC_FLUX_IN_bottom,3);
mz_mean_4p4z = mean(zoo3C_150m,3);
lz_mean_4p4z = mean(zoo4C_150m,3);
z_mean_4p4z = mean(zooC_150m,3);
mzloss_mean_4p4z = mean(zoo3_loss_150m,3);
lzloss_mean_4p4z = mean(zoo4_loss_150m,3);
zloss_mean_4p4z = mean(zoo_loss_150m,3);
mzloss_min_4p4z = min(zoo3_loss_150m,[],3);
lzloss_min_4p4z = min(zoo4_loss_150m,[],3);
zloss_min_4p4z = min(zoo_loss_150m,[],3);

Ztime = time;
Zyr = yr;

%%
clear TEMP_150m TEMP_bottom POC_FLUX_IN_bottom zoo3C_150m zoo3_loss_150m
clear zoo4C_150m zoo4_loss_150m
clear tp tb det_btm mz hploss_mz time yr

%% how many grid cells with zero min?
fid = (mzloss_min_fosi(:)==0);
id3 = (mzloss_min_4p4z(:)==0);
id4 = (lzloss_min_4p4z(:)==0);
zid = (zloss_min_4p4z(:)==0);

sum(fid)
sum(id3)
sum(id4)
sum(zid)

fid2 = double(mzloss_min_fosi==0);
id32 = double(mzloss_min_4p4z==0);
id42 = double(lzloss_min_4p4z==0);
zid2 = double(zloss_min_4p4z==0);

%%
clatlim=[-90 90];
clonlim=[-280 80];

LAT = double(TLAT);
LON = double(TLONG);

%% Zoo biomass
figure(1)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(mz_mean_fosi))
cmocean('tempo')
caxis([3 5])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('FOSI LZ')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(mz_mean_4p4z))
cmocean('tempo')
caxis([3 5])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('4P4Z Z3')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(z_mean_4p4z))
cmocean('tempo')
caxis([3 5])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('4P4Z Z3+Z4')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(lz_mean_4p4z))
cmocean('tempo')
caxis([3 5])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('4P4Z Z4')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[pp 'Map_FOSI_vs_4P4Z_zoo_biom.png'])

%% HP loss
figure(2)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(mzloss_mean_fosi))
cmocean('tempo')
caxis([-4 -1])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('FOSI LZ loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(mzloss_mean_4p4z))
cmocean('tempo')
caxis([-4 -1])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('4P4Z Z3 loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(zloss_mean_4p4z))
cmocean('tempo')
caxis([-4 -1])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('4P4Z Z3+Z4 loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(lzloss_mean_4p4z))
cmocean('tempo')
caxis([-4 -1])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('4P4Z Z4 loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[pp 'Map_FOSI_vs_4P4Z_zoo_loss.png'])

%% Det
figure(3)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(det_mean_fosi))
cmocean('tempo')
caxis([-4.5 -0.5])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('FOSI Det')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(det_mean_4p4z))
cmocean('tempo')
caxis([-4.5 -0.5])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('4P4Z Det')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[pp 'Map_FOSI_vs_4P4Z_det.png'])

%% HP loss min
figure(4)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(mzloss_min_fosi))
cmocean('tempo')
caxis([0 1e-4])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('FOSI LZ min loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(mzloss_min_4p4z))
cmocean('tempo')
caxis([0 1e-4])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('4P4Z Z3 min loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(zloss_min_4p4z))
cmocean('tempo')
caxis([0 1e-4])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('4P4Z Z3+Z4 min loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(lzloss_min_4p4z))
cmocean('tempo')
caxis([0 1e-4])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('4P4Z Z4 min loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[pp 'Map_FOSI_vs_4P4Z_min_zoo_loss.png'])

%% HP loss min = 0
figure(5)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,fid2)
colormap('turbo')
caxis([0 1])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('FOSI LZ min loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,id32)
colormap('turbo')
caxis([0 1])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('4P4Z Z3 min loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,zid2)
colormap('turbo')
caxis([0 1])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('4P4Z Z3+Z4 min loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,id42)
colormap('turbo')
caxis([0 1])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('4P4Z Z4 min loss')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%print('-dpng',[pp 'Map_FOSI_vs_4P4Z_min_zoo_loss.png'])


