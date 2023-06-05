% FEISTY forced by FOSI
% changed mzpref on S and M 
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
dpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% split P results
modP = 'v12_All_fish03_';
load([dpath 'FOSI_v12_ensemble_mzpref_all.mat']);

r_allP = r_all0;
mis_allP = mis_all0;
rmse_allP = rmse_all0;

for n=1:length(fx_all)
    fx_allP{n} = ['P' fx_all{n}];
end

clear rmse_all0 mis_all0 r_all0 fx_all

%% split B results
modB = 'v13_All_fish03_';
load([dpath 'FOSI_v13_ensemble_mzpref_all.mat']);

r_allB = r_all0;
mis_allB = mis_all0;
rmse_allB = rmse_all0;

for n=1:length(fx_all)
    fx_allB{n} = ['B' fx_all{n}];
end

clear rmse_all0 mis_all0 r_all0 fx_all

%% combine
mis_all0 = [mis_allP; mis_allB];
r_all0 = [r_allP,r_allB];
rmse_all0 = [rmse_allP,rmse_allB];
fx_all = [fx_allP,fx_allB];

%% save total ensemble to merge with other ensemble
save([dpath 'FOSI_v12_v13_ensemble_mzpref_all.mat'],'rmse_all0','mis_all0',...
    'r_all0','fx_all');

%% SAU comparison
% SAUP w/squid as forage, use top 10 catch yrs
load([spath 'forage_with_squids/SAUP_LME_Catch_top10_Stock_newF.mat']);
sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

l10all = [l10sF(keep);l10sP(keep);l10sD(keep)];
%variance of catch observations
sig = var(l10all);
%num of observations
n = length(l10all);

%% Outputs
%mean square error
mse_all = rmse_all0.^2;

%% Just look at F, P, and D ind
mis_all_F = mis_all0(:,:,2);
mis_all_F2 = mis_all_F;
mis_all_P = mis_all0(:,:,3);
mis_all_P2 = mis_all_P;
mis_all_D = mis_all0(:,:,4);

%put residuals of all fn types in one vector
mis_combo = [mis_all_F2,mis_all_P2,mis_all_D];

%% Classic AIC 
% AIC = -2*log(L) + 2*K
% log(L) = (-n/2) * log(2*pi*var) - (1/(2*var)) * sum(resid^2)

%logLike
LL = -n/2 * log(2*pi*sig) - (1/(2*sig)) * sum(mis_combo.^2,2);
caic_all = -2 * LL;

[caic_srt,idc] = sort(caic_all);
cdel = caic_srt - caic_srt(1);
cw = exp(-0.5*cdel) ./ sum( exp(-0.5*cdel) );

caicv(:,1) = idc;
caicv(:,2) = caic_srt;
caicv(:,3) = cdel;
caicv(:,4) = cw;
cT = array2table(caicv,'VariableNames',{'ParamSet','AIC','delta','weight'});
writetable(cT,[dpath 'AIC_FOSI_v12_v13_ensemble_mzpref.csv'])

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

writetable(pT,[dpath 'bestAIC_FOSI_v12_v13_ensemble_mzpref.csv'])

