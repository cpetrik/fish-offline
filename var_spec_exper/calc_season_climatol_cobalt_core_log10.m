% Calc seasonal climatology - Step 1 of powerspec slope
% log10-transform abund
% No transform environ (temp)
% Then create anomaly ts - Step 2 of powerspec slope

% Steps
% 1. Calc seasonal climatology
% 2. Create anomaly time series by removing climatol
% 3. Calc Theil-Sen linear trend
% 4. Remove trend if significant
% 5. Calc spectral slopes

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

%% CORE-forced output
load([fpath 'ocean_cobalt_temp100_monthly_1950_2007.mat'],'tp_100');
load([fpath 'ocean_cobalt_temp_btm_monthly_1950_2007.mat'],'tb');
load([fpath 'ocean_cobalt_mz100_monthly_1950_2007.mat'],'mz_100','units_vint');
load([fpath 'ocean_cobalt_lz100_monthly_1950_2007.mat'],'lz_100');
load([fpath 'ocean_cobalt_hploss_mz100_monthly_1950_2007.mat'],'hploss_mz_100');
load([fpath 'ocean_cobalt_hploss_lz100_monthly_1950_2007.mat'],'hploss_lz_100');
load([fpath 'ocean_cobalt_fndet_btm_monthly_1950_2007.mat']); %,'det_btm'

tp_100 = double(tp_100);
tb = double(tb);
mz_100 = double(mz_100);
lz_100 = double(lz_100);
det_btm = double(det_btm);
hploss_mz_100 = double(hploss_mz_100);
hploss_lz_100 = double(hploss_lz_100);

%% Convert Units to be same as fish (g WW m-2)
%det btm: mol N m-2 s-1
%zoo: mol N kg-1 integrated to mol N m-2 (1 kg of water = 1 m3)
%zoo loss: mol N m-2 s-1
%tp: degC
%tb: degC

mz = mz_100 * (106.0/16.0) * 12.01 * 9.0;
lz = lz_100 * (106.0/16.0) * 12.01 * 9.0;
mz_hp = hploss_mz_100 * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
lz_hp = hploss_lz_100 * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
det = det_btm * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;

%% By hand
yid = yr(runs);
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
    mo = m:12:nyr;
    tp_clim(:,:,m) = nanmean(tp_100(:,:,mo),3);
    tb_clim(:,:,m) = nanmean(tb(:,:,mo),3);
    mz_clim(:,:,m) = nanmean(log10(mz(:,:,mo)),3);
    lz_clim(:,:,m) = nanmean(log10(lz(:,:,mo)),3);
    mhp_clim(:,:,m) = nanmean(log10(mz_hp(:,:,mo)),3);
    lhp_clim(:,:,m) = nanmean(log10(lz_hp(:,:,mo)),3);
    det_clim(:,:,m) = nanmean(log10(det(:,:,mo)),3);
end

%% 
lat = double(geolat_t);
lon = double(geolon_t);

% Save climatologies
save([fpath 'cobalt_core_climatol_1950_2007_log10.mat'],...
    'units_vint','lat','lon','yid',...
    'tp_clim','tb_clim','mz_clim','lz_clim','mhp_clim','lhp_clim','det_clim');

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
    tp_anom(:,:,mo) = tp_100(:,:,mo) - tp_clim(:,:,m);
    tb_anom(:,:,mo) = tb(:,:,mo) - tb_clim(:,:,m);
    mz_anom(:,:,mo) = log10(mz_100(:,:,mo)) - mz_clim(:,:,m);
    lz_anom(:,:,mo) = log10(lz_100(:,:,mo)) - lz_clim(:,:,m);
    det_anom(:,:,mo) = log10(det_btm(:,:,mo)) - det_clim(:,:,m);
    hpmz_anom(:,:,mo) = log10(hploss_mz_100(:,:,mo)) - mhp_clim(:,:,m);
    hplz_anom(:,:,mo) = log10(hploss_lz_100(:,:,mo)) - lhp_clim(:,:,m);
end
 
%%
% mz_clim = repmat(mz_clim,1,1,(696/12));
% lz_clim = repmat(lz_clim,1,1,(696/12));
% mhp_clim = repmat(mhp_clim,1,1,(696/12));
% lhp_clim = repmat(lhp_clim,1,1,(696/12));
% det_clim = repmat(det_clim,1,1,(696/12));
% tp_clim = repmat(tp_clim,1,1,(696/12));
% tb_clim = repmat(tb_clim,1,1,(696/12));
% 
% %
% tp_anom = tp_100 - tp_clim;
% tb_anom = tb - tb_clim;
% mz_anom = sqrt(mz_100) - mz_clim;
% lz_anom = sqrt(lz_100) - lz_clim;
% det_anom = sqrt(det_btm) - det_clim;
% hpmz_anom = sqrt(hploss_mz_100) - mhp_clim;
% hplz_anom = sqrt(hploss_lz_100) - lhp_clim;

%% Save anom
save([fpath 'cobalt_core_anom_1950_2007_log10.mat'],...
    'units_vint','lat','lon','yid',...
    'tp_anom','tb_anom','mz_anom','lz_anom','hpmz_anom','hplz_anom','det_anom');
