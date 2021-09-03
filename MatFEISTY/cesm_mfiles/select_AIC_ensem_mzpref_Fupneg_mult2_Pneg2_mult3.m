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

%%
% %% param ensemble results with kt=0.0885
% dp  = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';
% load([dp 'Climatol_ensemble_param5_mid5.mat'],'rmse_all','mis_all','sim',...
%     'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch','r_all');
% r_all1 = r_all;
% rmse_all1 = rmse_all;
% mis_all1 = mis_all;
% sim1 = sim;
% lme_Fmcatch1 = lme_Fmcatch;
% lme_Pmcatch1 = lme_Pmcatch;
% lme_Dmcatch1 = lme_Dmcatch;
% clear rmse_all mis_all sim lme_Fmcatch lme_Pmcatch lme_Dmcatch r_all

%% orig results with mz = 1
mod = 'v12_All_fish03_';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
dpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];

load([dpath 'LME_fosi_fished_',mod,cfile '.mat'],'lme_mcatch','lme_area');

[r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch,lme_area);

%% param ensemble results with mz = 0.1:0.1:0.9
load([dpath 'FOSI_ensemble_mzpref_equal.mat'],'rmse_all','mis_all','mz',...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch','r_all');

fx_all = [1,mz];
r_all0 = [r',r_all];
rmse_all0 = [rmse',rmse_all];

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

%% Multiply the neg F upwelling LME misfits so they weigh more
up = [3;11;13;27;28;29];
[uboth,uid,kid]=intersect(up,keep);
mis_all_F = mis_all(:,:,2);
negF = mis_all_F(:,kid) < 0;
negF2 = mis_all_F;
negF3 = double(negF);
negF3(negF3==1) = 2;
negF3(negF3==0) = 1;
mis_all_F2 = mis_all_F;
mis_all_F2(:,kid) = mis_all_F(:,kid) .* negF3;

%% Multiply the P misfits < - log10(2) so they weigh more
% used to be log10(5) so they weigh more
mis_all_P = mis_all(:,:,3);
negP = mis_all_P < (-1*log10(2));
negP2 = mis_all_P;
negP3 = double(negP);
negP3(negP3==1) = 3;
negP3(negP3==0) = 1;
mis_all_P2 = mis_all_P .* negP3;

%% D ensem
mis_all_D = mis_all(:,:,4);
%put residuals of all fn types in one vector
mis_combo = [mis_all_F2,mis_all_P2,mis_all_D];

%% Do same to orig
mis_fn = mis(:,2:4);
nid = find(mis_fn(:,1) < 0);
unid = intersect(kid,nid);
mis_fn(unid,1) = mis_fn(unid,1) .* 2;

pid = find(mis_fn(:,2) < (-1*log10(2)));
mis_fn(pid,2) = mis_fn(pid,2) .* 3;

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
writetable(cT,[dpath 'AIC_Fupneg_mult2_Pneg2_mult3_FOSI_ensemble_mzpref_equal.csv'])

%% AICs <= AIC(orig) + 2
pset(:,1) = fx_all(idc);
pset(:,2:5) = caicv;
pset(:,6:10) = rmse_all0(:,idc)';
pset(:,11:15) = r_all0(:,idc)';

pT = array2table(pset,'VariableNames',{'eqMZpref','ParamSet','AIC','dAIC',...
    'wAIC','rmseAll','rmseF','rmseP','rmseD','rmsePD','rAll','rF','rP','rD','rPD'});
writetable(pT,[dpath 'bestAIC_Fupneg_mult2_Pneg2_mult3_FOSI_ensemble_mzpref_equal.csv'])

