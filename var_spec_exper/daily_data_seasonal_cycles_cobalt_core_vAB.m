% Use seasonal climatologies of phys & BGC to create
% forcing file separate from perturbations

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Climatologies
% non-temp vars are ln-trans
load([fpath 'cobalt_core_climatol_1950_2007_ln.mat']);

%% Transform back before interpolating
det_clim = exp(det_clim);
mz_clim = exp(mz_clim);
lz_clim = exp(lz_clim);
mhp_clim = exp(mhp_clim);
lhp_clim = exp(lhp_clim);

%% Just ocean cells
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat','GRD');

[ni,nj,nt] = size(tp_clim);
tp_clim_vec = reshape(tp_clim,ni*nj,nt);
tb_clim_vec = reshape(tb_clim,ni*nj,nt);
det_clim_vec = reshape(det_clim,ni*nj,nt);
mz_clim_vec = reshape(mz_clim,ni*nj,nt);
lz_clim_vec = reshape(lz_clim,ni*nj,nt);
mhp_clim_vec = reshape(mhp_clim,ni*nj,nt);
lhp_clim_vec = reshape(lhp_clim,ni*nj,nt);

tp_clim_vec = tp_clim_vec(GRD.ID,:);
tb_clim_vec = tb_clim_vec(GRD.ID,:);
det_clim_vec = det_clim_vec(GRD.ID,:);
mz_clim_vec = mz_clim_vec(GRD.ID,:);
lz_clim_vec = lz_clim_vec(GRD.ID,:);
mhp_clim_vec = mhp_clim_vec(GRD.ID,:);
lhp_clim_vec = lhp_clim_vec(GRD.ID,:);

%% Interpolate onto an annual cycle with steps equivalent to the number of
%model time steps in a year.
stepsperyear = 365;
interpx=linspace(0,12,stepsperyear);

%Interpolate variables and scale factors onto timesteps
tp_clim_spl = spline(0:1:12,[tp_clim_vec(:,12),tp_clim_vec],interpx);
tb_clim_spl = spline(0:1:12,[tb_clim_vec(:,12),tb_clim_vec],interpx);
det_clim_spl = spline(0:1:12,[det_clim_vec(:,12),det_clim_vec],interpx);
mz_clim_spl = spline(0:1:12,[mz_clim_vec(:,12),mz_clim_vec],interpx);
lz_clim_spl = spline(0:1:12,[lz_clim_vec(:,12),lz_clim_vec],interpx);
mhp_clim_spl = spline(0:1:12,[mhp_clim_vec(:,12),mhp_clim_vec],interpx);
lhp_clim_spl = spline(0:1:12,[lhp_clim_vec(:,12),lhp_clim_vec],interpx);


%%
COBALT.Tp = tp_clim_spl;
COBALT.Tb = tb_clim_spl;
COBALT.det = det_clim_spl;
COBALT.Zm = mz_clim_spl;
COBALT.Zl = lz_clim_spl;
COBALT.dZm = mhp_clim_spl;
COBALT.dZl = lhp_clim_spl;

save([fpath 'Data_CORE_cobalt_climatol_daily_1yr.mat'], 'COBALT');
