% Calc seasonal climatology of each biome
% Calc seasonal climatology of std deviations

% Steps - A Barton methods
% 1. Take sqrt of ts
% 2. Calc monthly climatology of sqrt ts
% 3. Calc anomaly ts (sqrt ts - clim ts)
% 4. Calc monthly climatology of std dev of anom ts

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

%% 1. Square root of everything
% Not temperature, because there are negatives
tp_sq = (tp_100);
tb_sq = (tb);
det_sq = sqrt(det);
mz_sq = sqrt(mz);
lz_sq = sqrt(lz);
mhp_sq = sqrt(mz_hp);
lhp_sq = sqrt(lz_hp);

%% Take means over each biome
for L=1:8
    lid = find(biome8_hist==L);
    Btp(L,:) = nanmean(tp_sq(lid,:));
    Btb(L,:) = nanmean(tb_sq(lid,:));
    Bdet(L,:) = nanmean(det_sq(lid,:));
    Bmz(L,:) = nanmean(mz_sq(lid,:));
    Blz(L,:) = nanmean(lz_sq(lid,:));
    Bhpmz(L,:) = nanmean(mhp_sq(lid,:));
    Bhplz(L,:) = nanmean(lhp_sq(lid,:));
end

%% 2. Monthly climatol by hand
yid = yr(runs);
nmo = length(yid);
nyr = length(yid)/12;

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
    mz_clim(:,m) = nanmean((Bmz(:,mo)),2);
    lz_clim(:,m) = nanmean((Blz(:,mo)),2);
    mhp_clim(:,m) = nanmean((Bhpmz(:,mo)),2);
    lhp_clim(:,m) = nanmean((Bhplz(:,mo)),2);
    det_clim(:,m) = nanmean((Bdet(:,mo)),2);
end

%% Look at seasonal diffs
figure(2)
plot(1:12,mz_clim);
legend('SLC','SECCS','SECSS','SCoast','NLC','NECCS','NECSS','NCoast')
legend('location','eastoutside')

%% Save climatologies
save([fpath 'cobalt_biome_core_climatol_1950_2007.mat'],...
    'units_vint','geolat_t','geolon_t','yid',...
    'tp_clim','tb_clim','mz_clim','lz_clim','mhp_clim','lhp_clim','det_clim');

%% 3. Calc anomaly ts
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
    mz_anom(:,mo) = (Bmz(:,mo)) - mz_clim(:,m);
    lz_anom(:,mo) = (Blz(:,mo)) - lz_clim(:,m);
    det_anom(:,mo) = (Bdet(:,mo)) - det_clim(:,m);
    hpmz_anom(:,mo) = (Bhpmz(:,mo)) - mhp_clim(:,m);
    hplz_anom(:,mo) = (Bhplz(:,mo)) - lhp_clim(:,m);
    
end

%% Look at seasonal diffs
figure(3)
plot(yid(1:12),tp_anom(:,1:12));
legend('SLC','SECCS','SECSS','SCoast','NLC','NECCS','NECSS','NCoast')
legend('location','eastoutside')
 
%% 4. Scale factor climatol 
mz_sc = nan*ones(8,12);
lz_sc = mz_sc;
mhp_sc = mz_sc;
lhp_sc = mz_sc;
det_sc = mz_sc;
tp_sc = mz_sc;
tb_sc = mz_sc;
for m = 1:12
    mo = m:12:nyr;
    tp_sc(:,m) = std(tp_anom(:,mo),0,2);
    tb_sc(:,m) = std(tb_anom(:,mo),0,2);
    mz_sc(:,m) = std((mz_anom(:,mo)),0,2);
    lz_sc(:,m) = std((lz_anom(:,mo)),0,2);
    mhp_sc(:,m) = std((hpmz_anom(:,mo)),0,2);
    lhp_sc(:,m) = std((hplz_anom(:,mo)),0,2);
    det_sc(:,m) = std((det_anom(:,mo)),0,2);
end

%% Save anom & scale factor ts
save([fpath 'cobalt_biome_core_anom_1950_2007.mat'],...
    'units_vint','geolat_t','geolon_t','yid',...
    'tp_anom','tb_anom','mz_anom','lz_anom','hpmz_anom','hplz_anom','det_anom',...
    'tp_sc','tb_sc','mz_sc','lz_sc','mhp_sc','lhp_sc','det_sc');
