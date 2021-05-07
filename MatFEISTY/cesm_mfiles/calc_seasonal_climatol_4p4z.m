% calc climatol b/c zoo loss too low

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/4P4Z/';
% FOSI
fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
% 4P4Z
zpath='/Volumes/MIP/GCM_DATA/CESM/4P4Z/';

%% 4P4Z 
load([zpath 'g.e22a06.G1850ECOIAF_JRA_PHYS_DEV.TL319_g17.4p4z.004.FIESTY-forcing.mat']);

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

%% reshape
[ni,nj,nt] = size(TEMP_150m);
Tp = double(reshape(TEMP_150m,ni*nj,nt));
Tb = double(reshape(TEMP_bottom,ni*nj,nt));
mz = double(reshape(zoo3C_150m,ni*nj,nt));
lz = double(reshape(zoo4C_150m,ni*nj,nt));
mzhp = double(reshape(zoo3_loss_150m,ni*nj,nt));
lzhp = double(reshape(zoo4_loss_150m,ni*nj,nt));
det = double(reshape(POC_FLUX_IN_bottom,ni*nj,nt));

z = mz+lz;
hp = mzhp + lzhp;

%% climatol
yid = (time/365); %days since 0001-01-01 
nmo = length(yid);
nyr = length(yid)/12;

mz_clim = nan*ones(ni*nj,12);
lz_clim = mz_clim;
mhp_clim = mz_clim;
lhp_clim = mz_clim;
det_clim = mz_clim;
tp_clim = mz_clim;
tb_clim = mz_clim;
z_clim = mz_clim;
hp_clim = mz_clim;
for m = 1:12
    mo = m:12:nmo;
    tp_clim(:,m) = nanmean(Tp(:,mo),2);
    tb_clim(:,m) = nanmean(Tb(:,mo),2);
    mz_clim(:,m) = nanmean(mz(:,mo),2);
    lz_clim(:,m) = nanmean(lz(:,mo),2);
    mhp_clim(:,m) = nanmean(mzhp(:,mo),2);
    lhp_clim(:,m) = nanmean(lzhp(:,mo),2);
    det_clim(:,m) = nanmean(det(:,mo),2);
    z_clim(:,m) = nanmean(z(:,mo),2);
    hp_clim(:,m) = nanmean(hp(:,mo),2);
end

%% take time means first
tp_mclim = nanmean(tp_clim);
tb_mclim = nanmean(tb_clim);
mz_mclim = nanmean(mz_clim);
lz_mclim = nanmean(lz_clim);
z_mclim = nanmean(z_clim);
hp_mclim = nanmean(hp_clim);
mzhp_mclim = nanmean(mhp_clim);
lzhp_mclim = nanmean(lhp_clim);
det_mclim = nanmean(det_clim);

hp_nclim = min(hp_clim);
mzhp_nclim = min(mhp_clim);
lzhp_nclim = min(lhp_clim);

%% reshape back
tp_clim = reshape(tp_clim,ni,nj,12);
tb_clim = reshape(tb_clim,ni,nj,12);
mz_clim = reshape(mz_clim,ni,nj,12);
lz_clim = reshape(lz_clim,ni,nj,12);
z_clim = reshape(z_clim,ni,nj,12);
hp_clim = reshape(hp_clim,ni,nj,12);
mzhp_clim = reshape(mhp_clim,ni,nj,12);
lzhp_clim = reshape(lhp_clim,ni,nj,12);
det_clim = reshape(det_clim,ni,nj,12);

%%
clatlim=[-90 90];
clonlim=[-280 80];

LAT = double(TLAT);
LON = double(TLONG);

%% all plots
for t=1:12
f1=figure(1);
subplot(4,3,t)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,squeeze(tp_clim(:,:,t)))
cmocean('thermal')
caxis([-2 35])
title(['Tp mo ' num2str(t)])
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
if t==12
    colorbar('Position',[0.25 0.05 0.4 0.03],'orientation','horizontal')
end
 
f2=figure(2);
subplot(4,3,t)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,squeeze(tb_clim(:,:,t)))
cmocean('thermal')
caxis([-2 15])
title(['Tb mo ' num2str(t)])
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
if t==12
    colorbar('Position',[0.25 0.05 0.4 0.03],'orientation','horizontal')
end

f3=figure(3);
subplot(4,3,t)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(squeeze(mz_clim(:,:,t))))
cmocean('tempo')
caxis([3 5])
title(['Z3 mo ' num2str(t)])
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
if t==12
    colorbar('Position',[0.25 0.05 0.4 0.03],'orientation','horizontal')
end

f4=figure(4);
subplot(4,3,t)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(squeeze(lz_clim(:,:,t))))
cmocean('tempo')
caxis([3 5])
title(['Z4 mo ' num2str(t)])
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
if t==12
    colorbar('Position',[0.25 0.05 0.4 0.03],'orientation','horizontal')
end

f5=figure(5);
subplot(4,3,t)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(squeeze(z_clim(:,:,t))))
cmocean('tempo')
caxis([3 5])
title(['Z3+Z4 mo ' num2str(t)])
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
if t==12
    colorbar('Position',[0.25 0.05 0.4 0.03],'orientation','horizontal')
end

f6=figure(6);
subplot(4,3,t)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(squeeze(mzhp_clim(:,:,t))))
cmocean('tempo')
caxis([-4 -1])
title(['Z3loss mo ' num2str(t)])
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
if t==12
    colorbar('Position',[0.25 0.05 0.4 0.03],'orientation','horizontal')
end

f7=figure(7);
subplot(4,3,t)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(squeeze(lzhp_clim(:,:,t))))
cmocean('tempo')
caxis([-4 -1])
title(['Z4loss mo ' num2str(t)])
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
if t==12
    colorbar('Position',[0.25 0.05 0.4 0.03],'orientation','horizontal')
end

f8=figure(8);
subplot(4,3,t)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(squeeze(hp_clim(:,:,t))))
cmocean('tempo')
caxis([-4 -1])
title(['Zloss mo ' num2str(t)])
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
if t==12
    colorbar('Position',[0.25 0.05 0.4 0.03],'orientation','horizontal')
end

f9=figure(9);
subplot(4,3,t)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,log10(squeeze(det_clim(:,:,t))))
cmocean('tempo')
caxis([-4 -1])
title(['Det mo ' num2str(t)])
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
if t==12
    colorbar('Position',[0.25 0.05 0.4 0.03],'orientation','horizontal')
end
end

%%
print(f1,'-dpng',[pp 'Map_4P4Z_Tp_clim.png'])
print(f2,'-dpng',[pp 'Map_4P4Z_Tb_clim.png'])
print(f3,'-dpng',[pp 'Map_4P4Z_MZ_clim.png'])
print(f4,'-dpng',[pp 'Map_4P4Z_LZ_clim.png'])
print(f5,'-dpng',[pp 'Map_4P4Z_Z_clim.png'])
print(f6,'-dpng',[pp 'Map_4P4Z_MZloss_clim.png'])
print(f7,'-dpng',[pp 'Map_4P4Z_LZloss_clim.png'])
print(f8,'-dpng',[pp 'Map_4P4Z_Zloss_clim.png'])
print(f9,'-dpng',[pp 'Map_4P4Z_Det_clim.png'])

%% HP loss min
figure(10)
subplot(3,3,1)
plot(1:12,tp_mclim,'k')
xlim([1 12])
title('TP')

subplot(3,3,2)
plot(1:12,tb_mclim,'k')
xlim([1 12])
title('TB')

subplot(3,3,3)
plot(1:12,det_mclim*1e-9 * 1e4 * 12.01 * 9.0*60*60*24,'k')
xlim([1 12])
title('Det')

subplot(3,3,4)
plot(1:12,mz_mclim* 1e-9 * 1e4 * 12.01 * 9.0,'k')
xlim([1 12])
title('MZ')

subplot(3,3,5)
plot(1:12,lz_mclim* 1e-9 * 1e4 * 12.01 * 9.0,'k')
xlim([1 12])
title('LZ')

subplot(3,3,6)
plot(1:12,z_mclim* 1e-9 * 1e4 * 12.01 * 9.0,'k')
xlim([1 12])
title('Z')

subplot(3,3,7)
plot(1:12,mzhp_mclim*1e-9 * 1e4 * 12.01 * 9.0*60*60*24,'k')
xlim([1 12])
title('MZloss')

subplot(3,3,8)
plot(1:12,lzhp_mclim*1e-9 * 1e4 * 12.01 * 9.0*60*60*24,'k')
xlim([1 12])
title('LZloss')

subplot(3,3,9)
plot(1:12,hp_mclim*1e-9 * 1e4 * 12.01 * 9.0*60*60*24,'k')
xlim([1 12])
title('Zloss')
print('-dpng',[pp 'ts_4P4Z_all_clim_gm2d.png'])

