% Use seasonal climatologies of phys & BGC to create 
% experimental time-series
% Log-transform abund
% No transform environ (temp)
% These are before adding noise of known spectral color


clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([cpath 'COBALT_hist_biomes_1950_2005.mat']);

%% Climatologies
load([fpath 'cobalt_biome_core_climatol_1950_2007.mat']);

% convert abund to non-log-transformed (same as fish (g WW m-2 or g WW m-2 d-1))
mz_clim = 10.^(mz_clim);
lz_clim = 10.^(lz_clim);
mhp_clim = 10.^(mhp_clim);
lhp_clim = 10.^(lhp_clim);
det_clim = 10.^(det_clim);

%% METHOD 3 - interp so start and end the same
Tdays=1:365;
Time=Tdays(15:30:end);

tp_sp = csape(Time,tp_clim,'periodic');
tp_day = fnval(tp_sp,1:365);

tb_sp = csape(Time,tb_clim,'periodic');
tb_day = fnval(tb_sp,1:365);

mz_sp = csape(Time,mz_clim,'periodic');
mz_day = fnval(mz_sp,1:365);

lz_sp = csape(Time,lz_clim,'periodic');
lz_day = fnval(lz_sp,1:365);

det_sp = csape(Time,det_clim,'periodic');
det_day = fnval(det_sp,1:365);

mhp_sp = csape(Time,mhp_clim,'periodic');
mhp_day = fnval(mhp_sp,1:365);

lhp_sp = csape(Time,lhp_clim,'periodic');
lhp_day = fnval(lhp_sp,1:365);

%%
figure(1)
subplot(3,3,1)
plot(Time,tp_clim,'o'); hold on
plot(1:365,tp_day)

subplot(3,3,2)
plot(Time,tb_clim,'o'); hold on
plot(1:365,tb_day)

subplot(3,3,3)
plot(Time,det_clim,'o'); hold on
plot(1:365,det_day)

subplot(3,3,4)
plot(Time,mz_clim,'o'); hold on
plot(1:365,mz_day)

subplot(3,3,5)
plot(Time,lz_clim,'o'); hold on
plot(1:365,lz_day)

subplot(3,3,7)
plot(Time,mhp_clim,'o'); hold on
plot(1:365,mhp_day)

subplot(3,3,8)
plot(Time,lhp_clim,'o'); hold on
plot(1:365,lhp_day)