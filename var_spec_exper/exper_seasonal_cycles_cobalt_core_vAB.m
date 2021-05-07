% Use seasonal climatologies of phys & BGC to create
% experimental time-series
% sqrt-transform abund
% No transform environ (temp)

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([cpath 'COBALT_hist_biomes_1950_2005.mat']);

%% Climatologies
% of sqrt variables
load([fpath 'cobalt_biome_core_climatol_1950_2007.mat']);

%% Scale factors
%of anomaly of sqrt variables
load([fpath 'cobalt_biome_core_anom_1950_2007.mat'],...
    'tp_sc','tb_sc','mz_sc','lz_sc','mhp_sc','lhp_sc','det_sc');

%% 5. Interpolate onto an annual cycle with steps equivalent to the number of
%model time steps in a year.
stepsperyear = 365;
interpx=linspace(0,12,stepsperyear);

%Interpolate variables and scale factors onto timesteps
tp_clim_spl = spline(0:1:12,[tp_clim(:,12),tp_clim],interpx);
tb_clim_spl = spline(0:1:12,[tb_clim(:,12),tb_clim],interpx);
det_clim_spl = spline(0:1:12,[det_clim(:,12),det_clim],interpx);
mz_clim_spl = spline(0:1:12,[mz_clim(:,12),mz_clim],interpx);
lz_clim_spl = spline(0:1:12,[lz_clim(:,12),lz_clim],interpx);
mhp_clim_spl = spline(0:1:12,[mhp_clim(:,12),mhp_clim],interpx);
lhp_clim_spl = spline(0:1:12,[lhp_clim(:,12),lhp_clim],interpx);

%Make sure variance is seasonal
tp_sc_spl = spline(0:1:12,[tp_sc(:,12),tp_sc],interpx);
tb_sc_spl = spline(0:1:12,[tb_sc(:,12),tb_sc],interpx);
det_sc_spl = spline(0:1:12,[det_sc(:,12),det_sc],interpx);
mz_sc_spl = spline(0:1:12,[mz_sc(:,12),mz_sc],interpx);
lz_sc_spl = spline(0:1:12,[lz_sc(:,12),lz_sc],interpx);
mhp_sc_spl = spline(0:1:12,[mhp_sc(:,12),mhp_sc],interpx);
lhp_sc_spl = spline(0:1:12,[lhp_sc(:,12),lhp_sc],interpx);

%%
Time = 1:12;

figure(1)
subplot(3,3,1)
plot(Time,tp_clim,'o'); hold on
plot(interpx,tp_clim_spl)

subplot(3,3,2)
plot(Time,tb_clim,'o'); hold on
plot(interpx,tb_clim_spl)

subplot(3,3,3)
plot(Time,det_clim,'o'); hold on
plot(interpx,det_clim_spl)

subplot(3,3,4)
plot(Time,mz_clim,'o'); hold on
plot(interpx,mz_clim_spl)

subplot(3,3,5)
plot(Time,lz_clim,'o'); hold on
plot(interpx,lz_clim_spl)

subplot(3,3,7)
plot(Time,mhp_clim,'o'); hold on
plot(interpx,mhp_clim_spl)

subplot(3,3,8)
plot(Time,lhp_clim,'o'); hold on
plot(interpx,lhp_clim_spl)

%% Use these as control forcing
% Repeat each biome 220 times = 20 random draws x 11 noise slopes

tp_clim_spl_rep=repelem(tp_clim_spl,220,1);
tb_clim_spl_rep=repelem(tb_clim_spl,220,1);
det_clim_spl_rep=(repelem(det_clim_spl,220,1)).^2;
mz_clim_spl_rep=(repelem(mz_clim_spl,220,1)).^2;
lz_clim_spl_rep=(repelem(lz_clim_spl,220,1)).^2;
mhp_clim_spl_rep=(repelem(mhp_clim_spl,220,1)).^2;
lhp_clim_spl_rep=(repelem(lhp_clim_spl,220,1)).^2;

%%
COBALT.Tp = tp_clim_spl_rep;
COBALT.Tb = tb_clim_spl_rep;
COBALT.det = det_clim_spl_rep;
COBALT.Zm = mz_clim_spl_rep;
COBALT.Zl = lz_clim_spl_rep;
COBALT.dZm = mhp_clim_spl_rep;
COBALT.dZl = lhp_clim_spl_rep;

save([fpath 'Data_core_cobalt_biome_climatol_daily_1yr.mat'], 'COBALT');

