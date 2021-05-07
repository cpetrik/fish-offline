% Calc seasonal climatology - Step 1 of powerspec slope
% ln-transform abund
% No transform environ (temp)
% Spatially average across biomes
% Then create anomaly ts - Step 2 of powerspec slope

% Steps
% 1. Calc seasonal climatology
% 2. Create anomaly time series by removing climatol
% 3. Calc Theil-Sen linear trend
% 4. Remove trend if significant
% 5. Calc spectral slopes

clear all
close all

%% output
fpath='/Volumes/MIP/GCM_DATA/ESM2M_PI/';
load([fpath 'ocean_cobalt_temp_1001_1400.mat']);
load([fpath 'ocean_cobalt_det_1001_1400.mat']);
load([fpath 'ocean_cobalt_zoo_1001_1400.mat']);
load([fpath 'ocean_cobalt_hploss_1001_1400.mat']);

% already converted units to be same as fish (g WW m-2)

%% Biomes 
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([cpath 'COBALT_hist_biomes_1950_2005.mat']);

biome = {'LC','ECCS','ECSS','Coast'};

lid = find(biome4_hist==1);
uid = find(biome4_hist==2);
sid = find(biome4_hist==3);
cid = find(biome4_hist==4);

%% Biome means
[ni,nj,nt] = size(det);
Tp = reshape(Tp,ni*nj,nt);
Tb = reshape(Tb,ni*nj,nt);
mz = reshape(mz,ni*nj,nt);
lz = reshape(lz,ni*nj,nt);
mzhp = reshape(mzhp,ni*nj,nt);
lzhp = reshape(lzhp,ni*nj,nt);
det = reshape(det,ni*nj,nt);

%%
tp_biome(1,:) = nanmean(Tp(lid,:),1);
tp_biome(2,:) = nanmean(Tp(uid,:),1);
tp_biome(3,:) = nanmean(Tp(sid,:),1);
tp_biome(4,:) = nanmean(Tp(cid,:),1);

tb_biome(1,:) = nanmean(Tb(lid,:),1);
tb_biome(2,:) = nanmean(Tb(uid,:),1);
tb_biome(3,:) = nanmean(Tb(sid,:),1);
tb_biome(4,:) = nanmean(Tb(cid,:),1);

mz_biome(1,:) = log(nanmean(mz(lid,:),1));
mz_biome(2,:) = log(nanmean(mz(uid,:),1));
mz_biome(3,:) = log(nanmean(mz(sid,:),1));
mz_biome(4,:) = log(nanmean(mz(cid,:),1));

lz_biome(1,:) = log(nanmean(lz(lid,:),1));
lz_biome(2,:) = log(nanmean(lz(uid,:),1));
lz_biome(3,:) = log(nanmean(lz(sid,:),1));
lz_biome(4,:) = log(nanmean(lz(cid,:),1));

mhp_biome(1,:) = log(nanmean(mzhp(lid,:),1));
mhp_biome(2,:) = log(nanmean(mzhp(uid,:),1));
mhp_biome(3,:) = log(nanmean(mzhp(sid,:),1));
mhp_biome(4,:) = log(nanmean(mzhp(cid,:),1));

lhp_biome(1,:) = log(nanmean(lzhp(lid,:),1));
lhp_biome(2,:) = log(nanmean(lzhp(uid,:),1));
lhp_biome(3,:) = log(nanmean(lzhp(sid,:),1));
lhp_biome(4,:) = log(nanmean(lzhp(cid,:),1));

det_biome(1,:) = log(nanmean(det(lid,:),1));
det_biome(2,:) = log(nanmean(det(uid,:),1));
det_biome(3,:) = log(nanmean(det(sid,:),1));
det_biome(4,:) = log(nanmean(det(cid,:),1));

z_biome(1,:) = log(nanmean(mz(lid,:)+lz(lid,:),1));
z_biome(2,:) = log(nanmean(mz(uid,:)+lz(uid,:),1));
z_biome(3,:) = log(nanmean(mz(sid,:)+lz(sid,:),1));
z_biome(4,:) = log(nanmean(mz(cid,:)+lz(cid,:),1));

hp_biome(1,:) = log(nanmean(mzhp(lid,:)+lzhp(lid,:),1));
hp_biome(2,:) = log(nanmean(mzhp(uid,:)+lzhp(uid,:),1));
hp_biome(3,:) = log(nanmean(mzhp(sid,:)+lzhp(sid,:),1));
hp_biome(4,:) = log(nanmean(mzhp(cid,:)+lzhp(cid,:),1));

%% By hand
yid = (time_all/365)+1; %days since 0001-01-01 
nmo = length(yid);
nyr = length(yid)/12;

%%
mz_clim = nan*ones(4,12);
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
    tp_clim(:,m) = nanmean(tp_biome(:,mo),2);
    tb_clim(:,m) = nanmean(tb_biome(:,mo),2);
    mz_clim(:,m) = nanmean(mz_biome(:,mo),2);
    lz_clim(:,m) = nanmean(lz_biome(:,mo),2);
    mhp_clim(:,m) = nanmean(mhp_biome(:,mo),2);
    lhp_clim(:,m) = nanmean(lhp_biome(:,mo),2);
    det_clim(:,m) = nanmean(det_biome(:,mo),2);
    z_clim(:,m) = nanmean(z_biome(:,mo),2);
    hp_clim(:,m) = nanmean(hp_biome(:,mo),2);
end

%% quick plot
figure
subplot(2,2,1)
plot(1:12,tp_clim)
subplot(2,2,2)
plot(1:12,tb_clim)
subplot(2,2,3)
plot(1:12,mz_clim)
subplot(2,2,4)
plot(1:12,lz_clim)

%% Calc anomaly ts
tp_anom = nan*ones(4,nmo);
tb_anom = tp_anom;
mz_anom = tp_anom;
lz_anom = tp_anom;
det_anom = tp_anom;
hpmz_anom = tp_anom;
hplz_anom = tp_anom;
z_anom = tp_anom;
hp_anom = tp_anom;

for m = 1:12
    mo = m:12:nmo;
    tp_anom(:,mo) = tp_biome(:,mo) - tp_clim(:,m);
    tb_anom(:,mo) = tb_biome(:,mo) - tb_clim(:,m);
    mz_anom(:,mo) = mz_biome(:,mo) - mz_clim(:,m);
    lz_anom(:,mo) = lz_biome(:,mo) - lz_clim(:,m);
    det_anom(:,mo) = det_biome(:,mo) - det_clim(:,m);
    hpmz_anom(:,mo) = mhp_biome(:,mo) - mhp_clim(:,m);
    hplz_anom(:,mo) = lhp_biome(:,mo) - lhp_clim(:,m);
    z_anom(:,mo) = z_biome(:,mo) - z_clim(:,m);
    hp_anom(:,mo) = hp_biome(:,mo) - hp_clim(:,m);
    
end

%%
figure
subplot(2,2,1)
plot(yid,mz_anom(1,:))
subplot(2,2,2)
plot(yid,mz_anom(2,:))
subplot(2,2,3)
plot(yid,mz_anom(3,:))
subplot(2,2,4)
plot(yid,mz_anom(4,:))

figure
subplot(2,2,1)
plot(yid,tp_anom(3,:))
subplot(2,2,2)
plot(yid,tb_anom(3,:))
subplot(2,2,3)
plot(yid,mz_anom(3,:))
subplot(2,2,4)
plot(yid,lz_anom(3,:))

%% Save climatologies
save([fpath 'cobalt_pi400_biomes_climatol_ln.mat'],...
    'yid','biome4_hist','biome','z_clim','hp_clim',...
    'tp_clim','tb_clim','mz_clim','lz_clim','mhp_clim','lhp_clim','det_clim');

% Save anom
save([fpath 'cobalt_pi400_biomes_temp_anom.mat'],...
    'biome4_hist','biome','yid','tp_anom','tb_anom','-v7.3');
save([fpath 'cobalt_pi400_biomes_det_anom_ln.mat'],...
    'biome4_hist','biome','yid','det_anom','-v7.3');
save([fpath 'cobalt_pi400_biomes_zoo_anom_ln.mat'],...
    'biome4_hist','biome','yid','mz_anom','lz_anom','z_anom','-v7.3');
save([fpath 'cobalt_pi400_biomes_hploss_anom_ln.mat'],...
    'biome4_hist','biome','yid','hpmz_anom','hplz_anom','hp_anom','-v7.3');


