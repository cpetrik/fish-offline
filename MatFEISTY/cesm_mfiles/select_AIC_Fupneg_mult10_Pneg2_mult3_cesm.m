% FEISTY diff parameter sets with CESM-4P2Z

clear 
close all

%% SAUP data
spath = '/Volumes/petrik-lab/Feisty/Obs_data/SAUP/';
load([spath 'Climatol_ms/Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

% SAUP w/squid as forage, use top 10 catch yrs
load([spath 'forage_with_squids/SAUP_LME_Catch_top10_Stock_newF.mat']);
sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);
l10all = [l10sF(keep);l10sP(keep);l10sD(keep)];

%variance of catch observations
sig = var(l10all);
%num of observations
n = length(l10all);

%% climatol baseline parameter set (bmet=0.175, A=0.5, MZ=0.9)
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
harv = 'All_fish03_1deg';
dpath = ['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/4P2Z/'];
load([dpath 'LME_4P2Z_',harv,'_' cfile '.mat'],'lme_mcatch');

% SAU comparison of baseline params
[r,rmse,ss,mis] = lme_saup_newF_corr_stock_ensem(lme_mcatch);
mis_all_F(1,:,:) = mis(:,2);
mis_all_P(1,:,:) = mis(:,3);
mis_all_D(1,:,:) = mis(:,4);

param(1,1) = 0.175;
param(1,2) = 0.50;
param(1,3) = 0.9;

%% bmet=0.175, A=0.75, MZ=0.9
clear lme_mcatch mis
cfile3 = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A075_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
npath3 = ['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile3 '/4P2Z/'];
load([npath3 'LME_4P2Z_',harv,'_' cfile3 '.mat'],'lme_mcatch');
[r,rmse,ss,mis] = lme_saup_newF_corr_stock_ensem(lme_mcatch);

mis_all_F(2,:,:) = mis(:,2);
mis_all_P(2,:,:) = mis(:,3);
mis_all_D(2,:,:) = mis(:,4);

param(2,1) = 0.175;
param(2,2) = 0.75;
param(2,3) = 0.9;

%% bmet=0.185, A=0.5, MZ=0.9
clear lme_mcatch mis
cfile4 = 'Dc_Lam700_enc70-b200_m400-b185-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
npath4 = ['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile4 '/4P2Z/'];
load([npath4 'LME_4P2Z_',harv,'_' cfile4 '.mat'],'lme_mcatch');
[r,rmse,ss,mis] = lme_saup_newF_corr_stock_ensem(lme_mcatch);

mis_all_F(3,:,:) = mis(:,2);
mis_all_P(3,:,:) = mis(:,3);
mis_all_D(3,:,:) = mis(:,4);

param(3,1) = 0.185;
param(3,2) = 0.50;
param(3,3) = 0.9;

%% bmet=0.185, A=0.65, MZ=0.9
clear lme_mcatch mis
cfile5 = 'Dc_Lam700_enc70-b200_m400-b185-k086_c20-b250_D075_A065_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
npath5 = ['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile5 '/4P2Z/'];
load([npath5 'LME_4P2Z_',harv,'_' cfile5 '.mat'],'lme_mcatch');
[r,rmse,ss,mis] = lme_saup_newF_corr_stock_ensem(lme_mcatch);

mis_all_F(4,:,:) = mis(:,2);
mis_all_P(4,:,:) = mis(:,3);
mis_all_D(4,:,:) = mis(:,4);

param(4,1) = 0.185;
param(4,2) = 0.65;
param(4,3) = 0.9;

%% bmet=0.195, A=0.5, MZ=0.9
clear lme_mcatch mis
cfile6 = 'Dc_Lam700_enc70-b200_m400-b195-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
npath6 = ['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile6 '/4P2Z/'];
load([npath6 'LME_4P2Z_',harv,'_' cfile6 '.mat'],'lme_mcatch');
[r,rmse,ss,mis] = lme_saup_newF_corr_stock_ensem(lme_mcatch);

mis_all_F(5,:,:) = mis(:,2);
mis_all_P(5,:,:) = mis(:,3);
mis_all_D(5,:,:) = mis(:,4);

param(5,1) = 0.195;
param(5,2) = 0.50;
param(5,3) = 0.9;

%% bmet=0.200, A=0.65, MZ=0.9
clear lme_mcatch mis
cfile8 = 'Dc_Lam700_enc70-b200_m400-b200-k086_c20-b250_D075_A065_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
npath8 = ['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile8 '/4P2Z/'];
load([npath8 'LME_4P2Z_',harv,'_' cfile8 '.mat'],'lme_mcatch');
[r,rmse,ss,mis] = lme_saup_newF_corr_stock_ensem(lme_mcatch);

mis_all_F(6,:,:) = mis(:,2);
mis_all_P(6,:,:) = mis(:,3);
mis_all_D(6,:,:) = mis(:,4);

param(6,1) = 0.200;
param(6,2) = 0.65;
param(6,3) = 0.9;

%% bmet=0.2125, A=0.5, MZ=0.9
clear lme_mcatch mis
cfile9 = 'Dc_Lam700_enc70-b213_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
npath9 = ['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile9 '/4P2Z/'];
load([npath9 'LME_4P2Z_',harv,'_' cfile9 '.mat'],'lme_mcatch');
[r,rmse,ss,mis] = lme_saup_newF_corr_stock_ensem(lme_mcatch);

mis_all_F(7,:,:) = mis(:,2);
mis_all_P(7,:,:) = mis(:,3);
mis_all_D(7,:,:) = mis(:,4);

param(7,1) = 0.2125;
param(7,2) = 0.50;
param(7,3) = 0.9;

%% bmet=0.2125, A=0.75, MZ=0.9
clear lme_mcatch mis
cfile10 = 'Dc_Lam700_enc70-b213_m400-b175-k086_c20-b250_D075_A075_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
npath10 = ['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile10 '/4P2Z/'];
load([npath10 'LME_4P2Z_',harv,'_' cfile10 '.mat'],'lme_mcatch');
[r,rmse,ss,mis] = lme_saup_newF_corr_stock_ensem(lme_mcatch);

mis_all_F(8,:,:) = mis(:,2);
mis_all_P(8,:,:) = mis(:,3);
mis_all_D(8,:,:) = mis(:,4);

param(8,1) = 0.2125;
param(8,2) = 0.75;
param(8,3) = 0.9;

% %% bmet=0.175, A=0.5, MZ=0.8  %CATCH NOT RUN
% clear lme_mcatch mis
% cfile1 = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ080_mMZ040_nmort1_BE08_CC80_RE00100';
% npath1 = ['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile1 '/4P2Z/'];
% load([npath1 'LME_4P2Z_',harv,'_' cfile1 '.mat'],'lme_mcatch');
% [r,rmse,ss,mis] = lme_saup_newF_corr_stock_ensem(lme_mcatch);
% 
% mis_all_F(9,:,:) = mis(:,2);
% mis_all_P(9,:,:) = mis(:,3);
% mis_all_D(9,:,:) = mis(:,4);
% 
% param(9,1) = 0.175;
% param(9,2) = 0.50;
% param(9,3) = 0.8;
% 
% %% bmet=0.175, A=0.5, MZ=0.7  %CATCH NOT RUN
% clear lme_mcatch mis
% cfile2 = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ070_mMZ035_nmort1_BE08_CC80_RE00100';
% npath2 = ['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile2 '/4P2Z/'];
% load([npath2 'LME_4P2Z_',harv,'_' cfile2 '.mat'],'lme_mcatch');
% [r,rmse,ss,mis] = lme_saup_newF_corr_stock_ensem(lme_mcatch);
% 
% mis_all_F(10,:,:) = mis(:,2);
% mis_all_P(10,:,:) = mis(:,3);
% mis_all_D(10,:,:) = mis(:,4);
% 
% param(10,1) = 0.175;
% param(10,2) = 0.50;
% param(10,3) = 0.7;
% 
% %% bmet=0.200, A=0.5, MZ=0.9   %CATCH NOT RUN
% clear lme_mcatch mis
% cfile7 = 'Dc_Lam700_enc70-b200_m400-b200-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
% npath7 = ['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile7 '/4P2Z/'];
% load([npath7 'LME_4P2Z_',harv,'_' cfile7 '.mat'],'lme_mcatch');
% [r,rmse,ss,mis] = lme_saup_newF_corr_stock_ensem(lme_mcatch);
% 
% mis_all_F(11,:,:) = mis(:,2);
% mis_all_P(11,:,:) = mis(:,3);
% mis_all_D(11,:,:) = mis(:,4);
% 
% param(11,1) = 0.200;
% param(11,2) = 0.50;
% param(11,3) = 0.9;

%% Multiply the neg F upwelling LME misfits so they weigh more ----------
up = [3;11;13;27;28;29];
[uboth,uid,kid]=intersect(up,keep);
negF = mis_all_F(:,kid) < 0;
negF2 = mis_all_F;
negF3 = double(negF);
negF3(negF3==1) = 10;
negF3(negF3==0) = 1;
mis_all_F2 = mis_all_F;
mis_all_F2(:,kid) = mis_all_F(:,kid) .* negF3;

%% Multiply the P misfits < - log10(2) so they weigh more
% used to be log10(5) so they weigh more
negP = mis_all_P < (-1*log10(2));
negP2 = mis_all_P;
negP3 = double(negP);
negP3(negP3==1) = 3;
negP3(negP3==0) = 1;
mis_all_P2 = mis_all_P .* negP3;

%% put residuals of all fn types in one vector
mis_combo = [mis_all_F2,mis_all_P2,mis_all_D];

%% Classic AIC 
% AIC = -2*log(L) + 2*K
% log(L) = (-n/2) * log(2*pi*var) - (1/(2*var)) * sum(resid^2)

%logLike
LL = -n/2 * log(2*pi*sig) - (1/(2*sig)) * sum(mis_combo.^2,2);

caic_all = -2 * LL;

[caic_srt,idc] = sort(caic_all);
%id_srt = id(idc);

%%
cdel = caic_srt - caic_srt(1);
cw = exp(-0.5*cdel) ./ sum( exp(-0.5*cdel) );

caicv(:,1) = idc;
caicv(:,2) = caic_srt;
caicv(:,3) = cdel;
caicv(:,4) = cw;
caicv(:,5) = param(idc,1);
caicv(:,6) = param(idc,2);
caicv(:,7) = param(idc,3);
%cT = array2table(caicv,'VariableNames',{'ParamSet','AIC','delta','weight'});
cT = array2table(caicv,'VariableNames',{'ID','AIC','delta','weight','bmet','A','MZ'});
writetable(cT,'/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/4P2Z_AIC_Fupneg_mult10_Pneg2_mult3.csv')

