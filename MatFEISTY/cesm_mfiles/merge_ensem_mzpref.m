% FEISTY forced by FOSI
% changed mzpref on S and M 
% initialized with spinup pref=1
% mean of all 68 yrs saved

clear all
close all

cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6.mat']);
load([cpath 'Data_grid_POP_gx1v6.mat']);

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
spath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% orig results with mz = 1
mod = 'v12_All_fish03_';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
dpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];

load([dpath 'LME_fosi_fished_',mod,cfile '.mat'],'lme_mcatch','lme_area');

[r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch,lme_area);

%% param ensemble results with mz = 0.1:0.1:0.9
load([dpath 'FOSI_v12_ensemble_mzpref_equal.mat'],'rmse_all','mis_all','mz',...
    'r_all');

mis_all0 = NaN*ones(10,45,5);
mis_all0(1,:,:) = mis;
mis_all0(2:10,:,:) = mis_all;
r_all0 = [r',r_all];
rmse_all0 = [rmse',rmse_all];

fx_all{1} = '1';
for n=1:length(mz)
    fx_all{n+1} = num2str(mz(n));
end

clear rmse_all mis_all mz r_all

%% param ensemble results with mz = 0.1:0.1:0.9 for M only
load([dpath 'FOSI_v12_ensemble_S1_Mmzpref.mat'],'rmse_all','mis_all','mz',...
    'r_all');
% load([dpath 'FOSI_v13_ensemble_mzpref_S1_Mmzpref.mat'],'rmse_all','mis_all','mz',...
%     'r_all');

mis_all0 = [mis_all0; mis_all];
r_all0 = [r_all0,r_all];
rmse_all0 = [rmse_all0,rmse_all];

for n=1:length(mz)
    fx_all{n+10} = ['M' num2str(mz(n))];
end

%% param ensemble results with mz = 0.1:0.1:0.9 for S only
load([dpath 'FOSI_v12_ensemble_M1_Smzpref.mat'],'rmse_all','mis_all','mz',...
    'r_all');
% load([dpath 'FOSI_v13_ensemble_mzpref_M1_Smzpref.mat'],'rmse_all','mis_all','mz',...
%     'r_all');

mis_all0 = [mis_all0; mis_all];
r_all0 = [r_all0,r_all];
rmse_all0 = [rmse_all0,rmse_all];

for n=1:length(mz)
    fx_all{n+19} = ['S' num2str(mz(n))];
end

%% param ensemble results with mz = 0.1:0.1:0.9 w/ M=0.5S
load([dpath 'FOSI_v12_ensemble_mzpref_MhalfS.mat'],'rmse_all','mis_all','mz',...
    'r_all');

mis_all0 = [mis_all0; mis_all];
r_all0 = [r_all0,r_all];
rmse_all0 = [rmse_all0,rmse_all];

for n=1:length(mz)
    fx_all{n+28} = ['S' num2str(mz(n)) '_M' num2str(mz(n)/2)];
end

%% param ensemble results with mz = 0.1:0.1:0.9 w/ S=0.5M
load([dpath 'FOSI_v12_ensemble_mzpref_ShalfM.mat'],'rmse_all','mis_all','mz',...
    'r_all');

mis_all0 = [mis_all0; mis_all];
r_all0 = [r_all0,r_all];
rmse_all0 = [rmse_all0,rmse_all];

for n=1:length(mz)
    fx_all{n+37} = ['S' num2str(mz(n)) '_M' num2str(mz(n)*2)];
end

%% param ensemble results with mz = 0.1:0.1:0.9 w/ S+M=1
load([dpath 'FOSI_v12_ensemble_mzpref_SMsum1.mat'],'rmse_all','mis_all','mz',...
    'r_all');

mis_all0 = [mis_all0; mis_all];
r_all0 = [r_all0,r_all];
rmse_all0 = [rmse_all0,rmse_all];

for n=1:length(mz)
    fx_all{n+46} = ['S' num2str(mz(n)) '_M' num2str(1-mz(n))];
end

%% save total ensemble to merge with other ensemble
save([dpath 'FOSI_v12_ensemble_mzpref_all.mat'],'rmse_all0','mis_all0',...
    'r_all0','fx_all');
