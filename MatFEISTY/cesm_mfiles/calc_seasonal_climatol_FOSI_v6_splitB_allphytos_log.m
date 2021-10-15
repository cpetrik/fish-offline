% Calc seasonal climatology of CESM FOSI

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

%% Inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'TLAT','TLONG','TAREA','time','yr');
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo_totloss_allphytoC.mat'],...
    'LzooC_150m','Lzoo_loss_150m');
load([fpath 'gridspec_POP_gx1v6_noSeas.mat'],'mask');

%% nans & zeros
TEMP_150m = double(TEMP_150m);
TEMP_bottom = double(TEMP_bottom);
POC_FLUX_IN_bottom = double(POC_FLUX_IN_bottom);

TEMP_bottom(TEMP_bottom >= 9.9e+36) = nan;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;
LzooC_150m(TEMP_bottom >= 9.9e+36) = nan;
Lzoo_loss_150m(TEMP_bottom >= 9.9e+36) = nan;

LzooC_150m(LzooC_150m<0) = 0.0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0.0;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom<0) = 0.0;

%% Units
%poc flux: mmol/m^3 cm/s
%zoo loss: mmol/m^3/s cm
%zoo: mmolC/m^3 cm
%tp: degC
%tb: degC

% From nmolC cm-2 s-1 to g(WW) m-2 d-1
% 1e9 nmol in 1 mol C
% 1e4 cm2 in 1 m2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
% * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

% From nmolC cm-2 s-1 to g(WW) m-2
% * 1e-9 * 1e4 * 12.01 * 9.0;

mz = LzooC_150m * 1e-9 * 1e4 * 12.01 * 9.0;
mz_hp = Lzoo_loss_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
det = POC_FLUX_IN_bottom * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

mz(mz(:)<=0) = eps;
mz_hp(mz_hp(:)<=0) = eps;
det(det(:)<=0) = eps;

%% By hand
[ni,nj,nmo] = size(det);
nyr = nmo/12;

%%
mz_clim = nan*ones(ni,nj,12);
mhp_clim = mz_clim;
det_clim = mz_clim;
tp_clim = mz_clim;
tb_clim = mz_clim;
for m = 1:12
    mo = m:12:nyr;
    tp_clim(:,:,m) = nanmean(TEMP_150m(:,:,mo),3);
    tb_clim(:,:,m) = nanmean(TEMP_bottom(:,:,mo),3);
%     mz_clim(:,:,m) = nanmean((mz(:,:,mo)),3);
%     mhp_clim(:,:,m) = nanmean((mz_hp(:,:,mo)),3);
%     det_clim(:,:,m) = nanmean((det(:,:,mo)),3);
    mz_clim(:,:,m) = nanmean(log(mz(:,:,mo)),3);
    mhp_clim(:,:,m) = nanmean(log(mz_hp(:,:,mo)),3);
    det_clim(:,:,m) = nanmean(log(det(:,:,mo)),3);
end

%% Save climatologies
save([fpath 'cesm_fosi_v6_climatol_68yrs_natlog.mat'],'TLAT','TLONG',...
    'tp_clim','tb_clim','mz_clim','mhp_clim','det_clim','nmo','nyr');

%% Quick check
tp1=reshape(tp_clim,ni*nj,12);
tb1=reshape(tb_clim,ni*nj,12);
mz1=reshape(mz_clim,ni*nj,12);
hp1=reshape(mhp_clim,ni*nj,12);
det1=reshape(det_clim,ni*nj,12);

% take time means first
tp_mclim = nanmean(tp1);
tb_mclim = nanmean(tb1);
mz_mclim = nanmean(mz1);
mhp_mclim = nanmean(hp1);
det_mclim = nanmean(det1);

figure
subplot(3,2,1)
plot(nanmean(tp1))
title('TP')

subplot(3,2,2)
plot(nanmean(tb1))
title('TB')

subplot(3,2,3)
plot(nanmean(mz1))
title('MZ')

subplot(3,2,4)
plot(nanmean(hp1))
title('MZloss')

subplot(3,2,5)
plot(nanmean(det1))
title('Det')

%%
clatlim=[-90 90];
clonlim=[-280 80];

LAT = double(TLAT);
LON = double(TLONG);

%% all plots
for t=1:12
f3=figure(3);
subplot(4,3,t)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(squeeze(mz_clim(:,:,t))))
cmocean('tempo')
caxis([-1 3])
title(['MZ mo ' num2str(t)])
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
if t==12
    colorbar('Position',[0.25 0.05 0.4 0.03],'orientation','horizontal')
end

f6=figure(6);
subplot(4,3,t)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(squeeze(mhp_clim(:,:,t))))
cmocean('tempo')
caxis([-2 1])
title(['MZloss mo ' num2str(t)])
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
if t==12
    colorbar('Position',[0.25 0.05 0.4 0.03],'orientation','horizontal')
end

f9=figure(9);
subplot(4,3,t)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,(squeeze(det_clim(:,:,t))))
cmocean('tempo')
caxis([-5 -2])
title(['Det mo ' num2str(t)])
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
if t==12
    colorbar('Position',[0.25 0.05 0.4 0.03],'orientation','horizontal')
end
end

%%
print(f3,'-dpng',[pp 'Map_FOSI_v6_MZ_clim_natlog.png'])
print(f6,'-dpng',[pp 'Map_FOSI_v6_MZloss_clim_natlog.png'])
print(f9,'-dpng',[pp 'Map_FOSI_v6_Det_clim_natlog.png'])

%% HP loss min
figure(10)
subplot(2,3,1)
plot(1:12,tp_mclim,'k')
xlim([1 12])
title('TP')

subplot(2,3,2)
plot(1:12,exp(mz_mclim),'k')
xlim([1 12])
title('MZ')

subplot(2,3,3)
plot(1:12,exp(mhp_mclim),'k')
xlim([1 12])
title('MZloss')

subplot(2,3,4)
plot(1:12,tb_mclim,'k')
xlim([1 12])
title('TB')

subplot(2,3,5)
plot(1:12,exp(det_mclim),'k')
xlim([1 12])
title('Det')
print('-dpng',[pp 'ts_FOSI_v6_all_clim_natlog_gm2d.png'])


