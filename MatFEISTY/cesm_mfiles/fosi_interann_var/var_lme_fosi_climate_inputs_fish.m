% Calc var of climate indices, inputs, and fish biomass of FEISTY by LME
% CESM FOSI mod = v15_All_fish03

clear all
close all

%% Climate anomalies
apath = '/Users/cpetrik/Dropbox/Princeton/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
%load([apath 'climate_anomalies.mat'])
load([apath 'Climate_anomalies_annual_means.mat'])

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%% FOSI input forcing
%load([cpath 'lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'])
spath='/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/FOSI/';
load([spath 'LME_fosi_input_anomalies_annual.mat']);

% put inputs in matrix
manom = nan*ones(5,66,68);
manom(1,:,:) = lme_tpa;
manom(2,:,:) = lme_tba;
manom(3,:,:) = lme_deta;
manom(4,:,:) = lme_mza;
manom(5,:,:) = lme_losa;

yanom = 1948:2015;

canom = {'Tp','Tb','Det','LZbiom','LZloss'};

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

dpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

sims = {'v15_All_fish03_';'v15_climatol_';'v15_varFood_';'v15_varTemp_'};
mod = sims{1};

%load([dpath 'LME_fosi_fished_',mod,cfile '.mat']);
fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
load([fpath 'LME_fosi_fished_',mod,'anomalies_annual.mat'])


%% mean & std by lme
% tlme = double(lme_mask);
% tlme(tlme<0) = nan;

lme_sf_std = std(lme_msfa,0,2,'omitnan');
lme_sp_std = std(lme_mspa,0,2,'omitnan');
lme_sd_std = std(lme_msda,0,2,'omitnan');
lme_mf_std = std(lme_mmfa,0,2,'omitnan');
lme_mp_std = std(lme_mmpa,0,2,'omitnan');
lme_md_std = std(lme_mmda,0,2,'omitnan');
lme_lp_std = std(lme_mlpa,0,2,'omitnan');
lme_ld_std = std(lme_mlda,0,2,'omitnan');
lme_b_std = std(lme_mba,0,2,'omitnan');
lme_a_std = std(lme_mAa,0,2,'omitnan');
lme_s_std = std(lme_mSa,0,2,'omitnan');
lme_m_std = std(lme_mMa,0,2,'omitnan');
lme_l_std = std(lme_mLa,0,2,'omitnan');
lme_f_std = std(lme_mFa,0,2,'omitnan');
lme_p_std = std(lme_mPa,0,2,'omitnan');
lme_d_std = std(lme_mDa,0,2,'omitnan');

lme_sf_mean = nanmean(lme_msfa,2);
lme_sp_mean = nanmean(lme_mspa,2);
lme_sd_mean = nanmean(lme_msda,2);
lme_mf_mean = nanmean(lme_mmfa,2);
lme_mp_mean = nanmean(lme_mmpa,2);
lme_md_mean = nanmean(lme_mmda,2);
lme_lp_mean = nanmean(lme_mlpa,2);
lme_ld_mean = nanmean(lme_mlda,2);
lme_b_mean = nanmean(lme_mba,2);
lme_a_mean = nanmean(lme_mAa,2);
lme_s_mean = nanmean(lme_mSa,2);
lme_m_mean = nanmean(lme_mMa,2);
lme_l_mean = nanmean(lme_mLa,2);
lme_f_mean = nanmean(lme_mFa,2);
lme_p_mean = nanmean(lme_mPa,2);
lme_d_mean = nanmean(lme_mDa,2);

lme_Tp_std = std(lme_tpa,0,2,'omitnan');
lme_Tb_std = std(lme_tba,0,2,'omitnan');
lme_Det_std = std(lme_deta,0,2,'omitnan');
lme_MZ_std = std(lme_mza,0,2,'omitnan');
lme_MZl_std = std(lme_losa,0,2,'omitnan');

lme_Tp_mean = nanmean(lme_tpa,2);
lme_Tb_mean = nanmean(lme_tba,2);
lme_Det_mean = nanmean(lme_deta,2);
lme_MZ_mean = nanmean(lme_mza,2);
lme_MZl_mean = nanmean(lme_losa,2);

%% Coefficient of variance
lme_sf_cv = lme_sf_std ./ lme_sf_mean;
lme_sp_cv = lme_sp_std ./ lme_sp_mean;
lme_sd_cv = lme_sd_std ./ lme_sd_mean;
lme_mf_cv = lme_mf_std ./ lme_mf_mean;
lme_mp_cv = lme_mp_std ./ lme_mp_mean;
lme_md_cv = lme_md_std ./ lme_md_mean;
lme_lp_cv = lme_lp_std ./ lme_lp_mean;
lme_ld_cv = lme_ld_std ./ lme_ld_mean;
lme_b_cv = lme_b_std ./ lme_b_mean;
lme_a_cv = lme_a_std ./ lme_a_mean;
lme_s_cv = lme_s_std ./ lme_s_mean;
lme_m_cv = lme_m_std ./ lme_m_mean;
lme_l_cv = lme_l_std ./ lme_l_mean;
lme_f_cv = lme_f_std ./ lme_f_mean;
lme_p_cv = lme_p_std ./ lme_p_mean;
lme_d_cv = lme_d_std ./ lme_d_mean;
% lme_pf_cv = lme_pf_std ./ lme_pf_mean;
% lme_pd_cv = lme_pd_std ./ lme_pd_mean;
% lme_lm_cv = lme_lm_std ./ lme_lm_mean;
