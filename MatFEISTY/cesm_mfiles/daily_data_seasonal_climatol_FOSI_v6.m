% Use seasonal climatologies of phys & BGC to create
% forcing file 

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';

%% Climatologies
% non-temp vars are ln-trans
load([fpath 'cesm_fosi_v6_climatol_68yrs_natlog.mat']);

%% Transform back before interpolating
det_clim = exp(det_clim);
mz_clim = exp(mz_clim);
mhp_clim = exp(mhp_clim);

%% Just ocean cells
load([fpath 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');

[ni,nj,nt] = size(tp_clim);
tp_clim_vec = reshape(tp_clim,ni*nj,nt);
tb_clim_vec = reshape(tb_clim,ni*nj,nt);
det_clim_vec = reshape(det_clim,ni*nj,nt);
mz_clim_vec = reshape(mz_clim,ni*nj,nt);
mhp_clim_vec = reshape(mhp_clim,ni*nj,nt);

tp_clim_vec = tp_clim_vec(GRD.ID,:);
tb_clim_vec = tb_clim_vec(GRD.ID,:);
det_clim_vec = det_clim_vec(GRD.ID,:);
mz_clim_vec = mz_clim_vec(GRD.ID,:);
mhp_clim_vec = mhp_clim_vec(GRD.ID,:);

%% Quick check
figure
subplot(3,2,1)
plot(nanmean(tp_clim_vec))
title('TP')

subplot(3,2,2)
plot(nanmean(tb_clim_vec))
title('TB')

subplot(3,2,3)
plot(nanmean(mz_clim_vec))
title('MZ')

subplot(3,2,4)
plot(nanmean(mhp_clim_vec))
title('MZloss')

subplot(3,2,5)
plot(nanmean(det_clim_vec))
title('Det')

%% Interpolate onto an annual cycle with steps equivalent to the number of
%model time steps in a year.
stepsperyear = 365;
interpx=linspace(0,12,stepsperyear);

%Interpolate variables and scale factors onto timesteps
tp_clim_spl = spline(0:1:12,[tp_clim_vec(:,12),tp_clim_vec],interpx);
tb_clim_spl = spline(0:1:12,[tb_clim_vec(:,12),tb_clim_vec],interpx);
det_clim_spl = spline(0:1:12,[det_clim_vec(:,12),det_clim_vec],interpx);
mz_clim_spl = spline(0:1:12,[mz_clim_vec(:,12),mz_clim_vec],interpx);
mhp_clim_spl = spline(0:1:12,[mhp_clim_vec(:,12),mhp_clim_vec],interpx);

%% Quick check
figure
subplot(3,2,1)
plot(nanmean(tp_clim_spl))
title('TP')

subplot(3,2,2)
plot(nanmean(tb_clim_spl))
title('TB')

subplot(3,2,3)
plot(nanmean(mz_clim_spl))
title('MZ')

subplot(3,2,4)
plot(nanmean(mhp_clim_spl))
title('MZloss')

subplot(3,2,5)
plot(nanmean(det_clim_spl))
title('Det')

%% All vars
ESM.Tp = tp_clim_spl;
ESM.Tb = tb_clim_spl;
ESM.det = det_clim_spl;
ESM.Zm = mz_clim_spl;
ESM.dZm = mhp_clim_spl;

save([fpath 'Data_cesm_fosi_v6_daily_climtol_1yr.mat'], 'ESM');

%% Temp vars
CESM.Tp = tp_clim_spl;
CESM.Tb = tb_clim_spl;

save([fpath 'Data_cesm_fosi_v6_daily_climtol_temp_1yr.mat'], 'CESM');

%% Food vars
clear CESM
CESM.det = det_clim_spl;
CESM.Zm = mz_clim_spl;
CESM.dZm = mhp_clim_spl;

save([fpath 'Data_cesm_fosi_v6_daily_climtol_food_1yr.mat'], 'CESM');
