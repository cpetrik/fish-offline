% Calc seasonal climatology - Step 1 of powerspec slope
% ln-transform abund
% No transform environ (temp)
% Then create anomaly ts - Step 2 of powerspec slope

% Steps
% 1. Calc seasonal climatology
% 2. Create anomaly time series by removing climatol
% 3. Calc Theil-Sen linear trend
% 4. Remove trend if significant
% 5. Calc spectral slopes

%clear all
clear average_DT average_T1 average_T2 btm_temp e20 fndet_btm geolat_c
clear geolon_c jhploss_nlgz_100 jhploss_nmdz_100 nlgz_100 nmdz_100 s20
clear time tp_100
close all

fpath='/Volumes/MIP/GCM_DATA/ESM2M_PI/';

%% output
load([fpath 'ocean_cobalt_temp_1001_1400.mat']);
load([fpath 'ocean_cobalt_det_1001_1400.mat']);
load([fpath 'ocean_cobalt_zoo_1001_1400.mat']);
load([fpath 'ocean_cobalt_hploss_1001_1400.mat']);

% already converted units to be same as fish (g WW m-2)

%% By hand
yid = (time_all/365)+1; %days since 0001-01-01 
nmo = length(yid);
nyr = length(yid)/12;

%%
[ni,nj] = size(geolon_t);
mz_clim = nan*ones(ni,nj,12);
lz_clim = mz_clim;
mhp_clim = mz_clim;
lhp_clim = mz_clim;
det_clim = mz_clim;
tp_clim = mz_clim;
tb_clim = mz_clim;
for m = 1:12
    mo = m:12:nmo;
    tp_clim(:,:,m) = nanmean(Tp(:,:,mo),3);
    tb_clim(:,:,m) = nanmean(Tb(:,:,mo),3);
    mz_clim(:,:,m) = nanmean(log(mz(:,:,mo)),3);
    lz_clim(:,:,m) = nanmean(log(lz(:,:,mo)),3);
    mhp_clim(:,:,m) = nanmean(log(mzhp(:,:,mo)),3);
    lhp_clim(:,:,m) = nanmean(log(lzhp(:,:,mo)),3);
    det_clim(:,:,m) = nanmean(log(det(:,:,mo)),3);
end

%% quick plot
figure
subplot(2,2,1)
plot(squeeze(mz_clim(84,5,:)))
subplot(2,2,2)
plot(squeeze(mz_clim(84,20,:)))
subplot(2,2,3)
plot(squeeze(mz_clim(84,50,:)))
subplot(2,2,4)
plot(squeeze(mz_clim(84,160,:)))

figure
subplot(2,2,1)
plot(squeeze(tp_clim(84,50,:)))
subplot(2,2,2)
plot(squeeze(tb_clim(84,50,:)))
subplot(2,2,3)
plot(squeeze(mz_clim(84,50,:)))
subplot(2,2,4)
plot(squeeze(lz_clim(84,50,:)))

%% 
lat = double(geolat_t);
lon = double(geolon_t);

%% Calc anomaly ts
tp_anom = nan*ones(ni,nj,nmo);
tb_anom = tp_anom;
mz_anom = tp_anom;
lz_anom = tp_anom;
det_anom = tp_anom;
hpmz_anom = tp_anom;
hplz_anom = tp_anom;

for m = 1:12
    mo = m:12:nmo;
    tp_anom(:,:,mo) = Tp(:,:,mo) - tp_clim(:,:,m);
    tb_anom(:,:,mo) = Tb(:,:,mo) - tb_clim(:,:,m);
    mz_anom(:,:,mo) = log(mz(:,:,mo)) - mz_clim(:,:,m);
    lz_anom(:,:,mo) = log(lz(:,:,mo)) - lz_clim(:,:,m);
    det_anom(:,:,mo) = log(det(:,:,mo)) - det_clim(:,:,m);
    hpmz_anom(:,:,mo) = log(mzhp(:,:,mo)) - mhp_clim(:,:,m);
    hplz_anom(:,:,mo) = log(lzhp(:,:,mo)) - lhp_clim(:,:,m);
    
end

%%
figure
subplot(2,2,1)
plot(yid,squeeze(mz_anom(84,5,:)))
subplot(2,2,2)
plot(yid,squeeze(mz_anom(84,20,:)))
subplot(2,2,3)
plot(yid,squeeze(mz_anom(84,50,:)))
subplot(2,2,4)
plot(yid,squeeze(mz_anom(84,160,:)))

figure
subplot(2,2,1)
plot(yid,squeeze(tp_anom(84,50,:)))
subplot(2,2,2)
plot(yid,squeeze(tb_anom(84,50,:)))
subplot(2,2,3)
plot(yid,squeeze(mz_anom(84,50,:)))
subplot(2,2,4)
plot(yid,squeeze(lz_anom(84,50,:)))

%% Save climatologies
save([fpath 'cobalt_pi400_climatol_ln.mat'],...
    'lat','lon','yid',...
    'tp_clim','tb_clim','mz_clim','lz_clim','mhp_clim','lhp_clim','det_clim');

% Save anom
save([fpath 'cobalt_pi400_temp_anom.mat'],...
    'lat','lon','yid','tp_anom','tb_anom','-v7.3');
save([fpath 'cobalt_pi400_det_anom_ln.mat'],...
    'lat','lon','yid','det_anom','-v7.3');
save([fpath 'cobalt_pi400_zoo_anom_ln.mat'],...
    'lat','lon','yid','mz_anom','lz_anom','-v7.3');
save([fpath 'cobalt_pi400_hploss_anom_ln.mat'],...
    'lat','lon','yid','hpmz_anom','hplz_anom','-v7.3');


