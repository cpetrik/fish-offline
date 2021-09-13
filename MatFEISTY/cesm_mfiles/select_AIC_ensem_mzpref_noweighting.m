% FEISTY forced by FOSI
% changed mzpref on S and M equally
% initialized with spinup pref=1
% mean of all 68 yrs saved
% AIC using SAU obs

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

mis_all0 = mis_all;
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

mis_all0 = [mis_all0; mis_all];
r_all0 = [r_all0,r_all];
rmse_all0 = [rmse_all0,rmse_all];

for n=1:length(mz)
    fx_all{n+10} = ['M' num2str(mz(n))];
end

%% param ensemble results with mz = 0.1:0.1:0.9 for S only
load([dpath 'FOSI_v12_ensemble_M1_Smzpref.mat'],'rmse_all','mis_all','mz',...
    'r_all');

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

%% SAU comparison
% SAUP data
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/SAUP_Stock_top10.mat');
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

l10all = [l10sF(keep);l10sP(keep);l10sD(keep)];
%variance of catch observations
sig = var(l10all);
%num of observations
n = length(l10all);

%% Outputs
%mean square error
mse_all = rmse_all.^2;
mse = rmse.^2;

%% Just look at F, P, and D ind
mis_all_F = mis_all0(:,:,2);
mis_all_F2 = mis_all_F;
mis_all_P = mis_all0(:,:,3);
mis_all_P2 = mis_all_P;
mis_all_D = mis_all0(:,:,4);

%put residuals of all fn types in one vector
mis_combo = [mis_all_F2,mis_all_P2,mis_all_D];

%% Do same to orig
mis_fn = mis(:,2:4);
mis_fn = reshape(mis_fn,45*3,1);

%% Classic AIC 
% AIC = -2*log(L) + 2*K
% log(L) = (-n/2) * log(2*pi*var) - (1/(2*var)) * sum(resid^2)

%logLike
LL = -n/2 * log(2*pi*sig) - (1/(2*sig)) * sum(mis_combo.^2,2);
LL0 = -n/2 * log(2*pi*sig) - (1/(2*sig)) * sum(mis_fn.^2);
LL_all = [LL0;LL];

caic_all = -2 * LL_all;
caic = -2 * LL0;

[caic_srt,idc] = sort(caic_all);
idc2 = find(caic_srt<caic);

cdel = caic_srt - caic_srt(1);
cw = exp(-0.5*cdel) ./ sum( exp(-0.5*cdel) );

caicv(:,1) = idc;
caicv(:,2) = caic_srt;
caicv(:,3) = cdel;
caicv(:,4) = cw;
cT = array2table(caicv,'VariableNames',{'ParamSet','AIC','delta','weight'});
writetable(cT,[dpath 'AIC_FOSI_v12_ensemble_mzpref.csv'])

%% AICs <= AIC(orig) + 2
pset(:,1) = nan*ones(length(idc),1);
pset(:,2:5) = caicv;
pset(:,6:10) = rmse_all0(:,idc)';
pset(:,11:15) = r_all0(:,idc)';

pT = array2table(pset,'VariableNames',{'MZpref','ParamSet','AIC','dAIC',...
    'wAIC','rmseAll','rmseF','rmseP','rmseD','rmsePD','rAll','rF','rP','rD','rPD'});
% test = char(fx_all);
% test2 = test(idc,:);
% pT(:,1) = test2;
pT(:,16) = fx_all(idc)';

writetable(pT,[dpath 'bestAIC_FOSI_ensemble_mzpref.csv'])

