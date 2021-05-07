% Calc temporal variance of CORE-forced Cobalt inputs

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

%% Just ocean cells
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat','GRD');

[ni,nj,nt] = size(tp_100);
tp_100_vec = reshape(tp_100,ni*nj,nt);
tb_vec = reshape(tb,ni*nj,nt);
det_vec = reshape(det,ni*nj,nt);
mz_vec = reshape(mz,ni*nj,nt);
lz_vec = reshape(lz,ni*nj,nt);
mz_hp_vec = reshape(mz_hp,ni*nj,nt);
lz_hp_vec = reshape(lz_hp,ni*nj,nt);

tp_100_vec = tp_100_vec(GRD.ID,:);
tb_vec = tb_vec(GRD.ID,:);
det_vec = det_vec(GRD.ID,:);
mz_vec = mz_vec(GRD.ID,:);
lz_vec = lz_vec(GRD.ID,:);
mz_hp_vec = mz_hp_vec(GRD.ID,:);
lz_hp_vec = lz_hp_vec(GRD.ID,:);

%% variance
tp_100_var = var(tp_100_vec,0,2);
tb_var = var(tb_vec,0,2);
det_var = var(det_vec,0,2);
mz_var = var(mz_vec,0,2);
lz_var = var(lz_vec,0,2);
mz_hp_var = var(mz_hp_vec,0,2);
lz_hp_var = var(lz_hp_vec,0,2);

tp_100_mean = mean(tp_100_vec,2);

%%
save([fpath 'cobalt_core_variance_1950_2007.mat'],...
    'geolat_t','geolon_t','GRD',...
    'tp_100_var','tb_var','det_var','mz_var','lz_var','mz_hp_var',...
    'lz_hp_var','tp_100_mean');

