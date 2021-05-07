% Calc seasonal climatology - Step 1 of powerspec slope
% Log-transform abund
% No transform environ (temp)
% Then create anomaly ts - Step 2 of powerspec slope
% For 4 biomes and each hemisphere separately

% Steps
% 1. Calc seasonal climatology
% 2. Create anomaly time series by removing climatol
% 3. Calc Theil-Sen linear trend
% 4. Remove trend if significant
% 5. Calc spectral slopes

% OLD VERSION, NOT CORRECT METHODS

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

% Biomes
load([cpath 'COBALT_hist_biomes_1950_2005.mat']);
% Grid info
load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

%% Sep hemis
nid = find(geolat_t(:)>0);
sid = find(geolat_t(:)<0);
vmask = lmask;
vmask(nid) = vmask(nid)*2;
nhem = find(vmask==2);
shem = find(vmask==1);

biome8_hist = biome4_hist;
biome8_hist(nid) = biome8_hist(nid) +4;
%1 & 5: LC
%2 & 6: ECCS
%3 & 7: ECSS
%4 & 8: Coastal
pcolor(biome8_hist)
shading flat

%% CORE-forced output
load([fpath 'ocean_cobalt_temp100_monthly_1950_2007.mat'],'tp_100');
load([fpath 'ocean_cobalt_temp_btm_monthly_1950_2007.mat'],'tb');
load([fpath 'ocean_cobalt_mz100_monthly_1950_2007.mat'],'mz_100','units_vint');
load([fpath 'ocean_cobalt_lz100_monthly_1950_2007.mat'],'lz_100');
load([fpath 'ocean_cobalt_hploss_mz100_monthly_1950_2007.mat'],'hploss_mz_100');
load([fpath 'ocean_cobalt_hploss_lz100_monthly_1950_2007.mat'],'hploss_lz_100');
load([fpath 'ocean_cobalt_fndet_btm_monthly_1950_2007.mat']); %,'det_btm'

[ni,nj,nt] = size(tb);

tp_100 = double(reshape(tp_100,ni*nj,nt));
tb = double(reshape(tb,ni*nj,nt));
mz_100 = double(reshape(mz_100,ni*nj,nt));
lz_100 = double(reshape(lz_100,ni*nj,nt));
det_btm = double(reshape(det_btm,ni*nj,nt));
hploss_mz_100 = double(reshape(hploss_mz_100,ni*nj,nt));
hploss_lz_100 = double(reshape(hploss_lz_100,ni*nj,nt));

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

%% Take means over each biome

for L=1:8
    lid = find(biome8_hist==L);
    Btp(L,:) = nanmean(tp_100(lid,:));
    Btb(L,:) = nanmean(tb(lid,:));
    Bdet(L,:) = nanmean(det(lid,:));
    Bmz(L,:) = nanmean(mz(lid,:));
    Blz(L,:) = nanmean(lz(lid,:));
    Bhpmz(L,:) = nanmean(mz_hp(lid,:));
    Bhplz(L,:) = nanmean(lz_hp(lid,:));
end

%% By hand
yid = yr(runs);
nmo = length(yid);
nyr = length(yid)/12;

%%
mz_clim = nan*ones(8,12);
lz_clim = mz_clim;
mhp_clim = mz_clim;
lhp_clim = mz_clim;
det_clim = mz_clim;
tp_clim = mz_clim;
tb_clim = mz_clim;
for m = 1:12
    mo = m:12:nyr;
    tp_clim(:,m) = nanmean(Btp(:,mo),2);
    tb_clim(:,m) = nanmean(Btb(:,mo),2);
    mz_clim(:,m) = nanmean(log10(Bmz(:,mo)),2);
    lz_clim(:,m) = nanmean(log10(Blz(:,mo)),2);
    mhp_clim(:,m) = nanmean(log10(Bhpmz(:,mo)),2);
    lhp_clim(:,m) = nanmean(log10(Bhplz(:,mo)),2);
    det_clim(:,m) = nanmean(log10(Bdet(:,mo)),2);
end

%% Look at seasonal diffs
figure(2)
plot(1:12,mz_clim);
legend('SLC','SECCS','SECSS','SCoast','NLC','NECCS','NECSS','NCoast')
legend('location','eastoutside')

%% Calc anomaly ts
tp_anom = nan*ones(8,nmo);
tb_anom = tp_anom;
mz_anom = tp_anom;
lz_anom = tp_anom;
det_anom = tp_anom;
hpmz_anom = tp_anom;
hplz_anom = tp_anom;

for m = 1:12
    mo = m:12:nmo;
    tp_anom(:,mo) = Btp(:,mo) - tp_clim(:,m);
    tb_anom(:,mo) = Btb(:,mo) - tb_clim(:,m);
    mz_anom(:,mo) = log10(Bmz(:,mo)) - mz_clim(:,m);
    lz_anom(:,mo) = log10(Blz(:,mo)) - lz_clim(:,m);
    det_anom(:,mo) = log10(Bdet(:,mo)) - det_clim(:,m);
    hpmz_anom(:,mo) = log10(Bhpmz(:,mo)) - mhp_clim(:,m);
    hplz_anom(:,mo) = log10(Bhplz(:,mo)) - lhp_clim(:,m);
    
end

