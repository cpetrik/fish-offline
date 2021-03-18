% Calc seasonal climatology of CESM FOSI

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';

%% Inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'TLAT','TLONG','TAREA','time','yr');
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo.mat'],...
    'LzooC_150m','Lzoo_loss_150m');
load([fpath 'gridspec_POP_gx1v6.mat'],'mask');

%% nans
TEMP_150m = double(TEMP_150m);
TEMP_bottom = double(TEMP_bottom);
POC_FLUX_IN_bottom = double(POC_FLUX_IN_bottom);

TEMP_bottom(TEMP_bottom >= 9.9e+36) = nan;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;

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
    mz_clim(:,:,m) = nanmean(log(mz(:,:,mo)),3);
    mhp_clim(:,:,m) = nanmean(log(mz_hp(:,:,mo)),3);
    det_clim(:,:,m) = nanmean(log(det(:,:,mo)),3);
end

%%
lat = double(TLAT);
lon = double(TLONG);

% Save climatologies
save([fpath 'cesm_fosi_climatol_68yrs_ln.mat'],'lat','lon',...
    'tp_clim','tb_clim','mz_clim','mhp_clim','det_clim','nmo','nyr');

%% Quick check
tp1=reshape(tp_clim,ni*nj,12);
tb1=reshape(tb_clim,ni*nj,12);
mz1=reshape(mz_clim,ni*nj,12);
hp1=reshape(mhp_clim,ni*nj,12);
det1=reshape(det_clim,ni*nj,12);

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


