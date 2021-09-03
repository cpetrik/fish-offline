% All correlations
clear all
close all

%% CESM FOSI grid
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6.mat']);
load([cpath 'Data_grid_POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;
hID = ID;

load([cpath 'LME-mask-POP_gx1v6.mat']);
load([cpath 'lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'],...
    'lme_tp_fosi','lme_area');

hlme_ptemp = lme_tp_fosi;

hlme_area = lme_area;
hlme_area_km2 = lme_area * 1e-6;
clear lme_area

hlme = lme_mask;

%% FEISTY Output
cfile2 = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
mod = 'v12_All_fish03_';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
ppath = [pp cfile2 '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

% CESM
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile2 '/'];
load([fpath 'LME_fosi_fished_',mod,cfile2 '.mat'],'lme_mcatch',...
    'lme_mbio');

%%
hlme_mcatch = nansum(lme_mcatch,2) * 1e-6;
hlme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
hlme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
hlme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
hlme_Bmbio = lme_mbio(:,9) * 1e-6;
hlme_Mmbio = (lme_mbio(:,4) + lme_mbio(:,5) + lme_mbio(:,6)) * 1e-6;
hlme_Lmbio = (lme_mbio(:,7) + lme_mbio(:,8)) * 1e-6;

% MT/km2
hlme_mcatch = hlme_mcatch ./ hlme_area_km2;
hlme_Fmcatch = hlme_Fmcatch ./ hlme_area_km2;
hlme_Pmcatch = hlme_Pmcatch ./ hlme_area_km2;
hlme_Dmcatch = hlme_Dmcatch ./ hlme_area_km2;
hlme_Bmbio = hlme_Bmbio ./ hlme_area_km2;
hlme_Mmbio = hlme_Mmbio ./ hlme_area_km2;
hlme_Lmbio = hlme_Lmbio ./ hlme_area_km2;

pFracPD = hlme_Pmcatch ./ (hlme_Pmcatch + hlme_Dmcatch);
pFracPF = hlme_Pmcatch ./ (hlme_Pmcatch + hlme_Fmcatch);
pFracLM = hlme_Lmbio ./ (hlme_Lmbio + hlme_Mmbio);

l10p=log10(hlme_mcatch);
l10pF=log10(hlme_Fmcatch);
l10pP=log10(hlme_Pmcatch);
l10pD=log10(hlme_Dmcatch);
l10pB=log10(hlme_Bmbio);

% l10ATL = log10(lme_te(:,2));
% l10HTL = log10(lme_te(:,3));
% l10LTL = log10(lme_te(:,4));

%% SAUP: All, F, P, D, Frac P
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%% Stats
%r
[rall,pall]=corr(l10s(keep),l10p(keep));
[rF,pF]=corr(l10sF(keep),l10pF(keep));
[rP,pP]=corr(l10sP(keep),l10pP(keep));
[rD,pD]=corr(l10sD(keep),l10pD(keep));
[rPD,pPD]=corr(sFracPD(keep),pFracPD(keep));

%root mean square error
o=l10s(keep);
p=l10p(keep);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=l10sF(keep);
p=l10pF(keep);
n = length(o);
num=nansum((p-o).^2);
rmseF = sqrt(num/n);

o=l10sP(keep);
p=l10pP(keep);
n = length(o);
num=nansum((p-o).^2);
rmseP = sqrt(num/n);

o=l10sD(keep);
p=l10pD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseD = sqrt(num/n);

o=sFracPD(keep);
p=pFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmsePD = sqrt(num/n);

%Fmed
Fall=10^(median(l10s(keep)-l10p(keep)));
FF=10^(median(l10sF(keep)-l10pF(keep)));
FP=10^(median(l10sP(keep)-l10pP(keep)));
FD=10^(median(l10sD(keep)-l10pD(keep)));
FPD=10^(median(sFracPD(keep)-pFracPD(keep)));

% Bias (FOSI minus SAUP)
%average error = bias
%skill(3,c) = nansum(p-o) / n;
p=l10p(keep);
o=l10s(keep);
n = length(o);
bias = nansum(o-p) / n;

p=l10pF(keep);
o=l10sF(keep);
n = length(o);
biasF = nansum(o-p) / n;

p=l10pP(keep);
o=l10sP(keep);
n = length(o);
biasP = nansum(o-p) / n;

p=l10pD(keep);
o=l10sD(keep);
n = length(o);
biasD = nansum(o-p) / n;

p=pFracPD(keep);
o=sFracPD(keep);
n = length(o);
biasPD = nansum(o-p) / n;

% o=l10pATL(keep);
% p=l10sATL(keep);
% n = length(o);
% biasATL = nansum(o-p) / n;
% 
% o=l10pHTL(keep);
% p=l10sHTL(keep);
% n = length(o);
% biasHTL = nansum(o-p) / n;
% 
% o=l10pLTL(keep);
% p=l10sLTL(keep);
% n = length(o);
% biasLTL = nansum(o-p) / n;

%% DvD: Frac P
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat')

did=[1:61,63];
did2 = notLELC(notLELC<=63);

%r
[rDvD,pDvD]=corr(FracLP(did),pFracPD(did));
[rDvD2,pDvD2]=corr(FracLP(did2),pFracPD(did2));

%root mean square error
o=FracLP(did);
p=pFracPD(did);
n = length(o);
num=nansum((p-o).^2);
rmseDvD = sqrt(num/n);

o=FracLP(did2);
p=pFracPD(did2);
n = length(o);
num=nansum((p-o).^2);
rmseDvD2 = sqrt(num/n);

%Fmed
FDvD=10^(median(FracLP(did)-pFracPD(did)));
FDvD2=10^(median(FracLP(did2)-pFracPD(did2)));

%Bias
o=FracLP(did);
p=pFracPD(did);
n = length(o);
biasDvD = nansum(o-p) / n;

o=FracLP(did2);
p=pFracPD(did2);
n = length(o);
biasDvD2 = nansum(o-p) / n;

%% Stock: All
%Reconciling Fisheries Catch and Ocean Productivity
%***TEMPLATE FOR FEEDBACK, PENDING FINAL CHECKS***
%model: 4
%alpha: 0.14
%TE_p: 0.13
%TE_b: 0.40
%f_T: 0.74
%T_100,warm: 19.99
%All fluxes in g C m-2 yr-1, Temperature in degrees celsius
%cols = LME  NPP   MESOZP  FDET   TLeq     T  modcatch SAUcatch
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'])

% POEM in gC/m2
%gWW
wlme_mcatch = nansum(lme_mcatch,2);
%gWW/m2
wmlme_mcatch = wlme_mcatch ./ hlme_area;
%gC/m2
glme_mcatch = wmlme_mcatch / 9;

%r
[rPNAS,pPNAS]=corr(StockPNAS(:,7),glme_mcatch(keep));

%root mean square error
o=StockPNAS(:,7);
p=glme_mcatch(keep);
n = length(o);
num=nansum((p-o).^2);
rmsePNAS = sqrt(num/n);

%Fmed
FPNAS=10^(median(StockPNAS(:,7)-glme_mcatch(keep)));

o=StockPNAS(:,7);
p=glme_mcatch(keep);
n = length(o);
biasPNAS = nansum(o-p) / n;

x=-5:0.5:5;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);
figure(1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(StockPNAS(:,7),glme_mcatch(keep),20,lme_tp_fosi(keep,1),'filled'); hold on;
cmocean('thermal');
text(-0.25,1.4,['r = ' sprintf('%2.2f',rPNAS) ' (p = ' sprintf('%2.2f',pPNAS) ')'])
text(-0.25,1.3,['RMSE = ' sprintf('%2.2f',rmsePNAS)])
axis([-0.5 1.5 -0.5 1.5])
xlabel('Stock PNAS')
ylabel('FOSI CESM')
title('All fishes')
stamp([harv '_' cfile2])
print('-dpng',[ppath 'FOSI_StockPNAS_',mod,'comp_temp.png'])

%% Mauread TE eff
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/';
load([spath 'Maureaud_etal_2017_s002_ECI.mat']);

% ECI for clim years (1991-1995?)
mECI = mean(ECI(:,2:6),2);
mid = ECI(:,1);

% FEISTY LME TEeffs
pECI = lme_te(mid,3);

Lma = log10(mECI); %log10(mECI)
Lpo = log10(pECI);

%r
[rL,pL]=corr(Lma,Lpo);

%root mean square error
o=Lma;
p=Lpo;
n = length(o);
num=nansum((p-o).^2);
rmseL = sqrt(num/n);

%Fmed
FL=10^(median(Lma-Lpo));

%% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rF;
fish_stat(3,1) = rP;
fish_stat(4,1) = rD;
fish_stat(5,1) = rPD;
fish_stat(6,1) = rDvD2;
fish_stat(7,1) = rPNAS;
fish_stat(1,2) = pall;
fish_stat(2,2) = pF;
fish_stat(3,2) = pP;
fish_stat(4,2) = pD;
fish_stat(5,2) = pPD;
fish_stat(6,2) = pDvD2;
fish_stat(7,2) = pPNAS;
fish_stat(1,3) = rmse;
fish_stat(2,3) = rmseF;
fish_stat(3,3) = rmseP;
fish_stat(4,3) = rmseD;
fish_stat(5,3) = rmsePD;
fish_stat(6,3) = rmseDvD2;
fish_stat(7,3) = rmsePNAS;
fish_stat(1,4) = bias;
fish_stat(2,4) = biasF;
fish_stat(3,4) = biasP;
fish_stat(4,4) = biasD;
fish_stat(5,4) = biasPD;
fish_stat(6,4) = biasDvD2;
fish_stat(7,4) = biasPNAS;

% save
Fstat = array2table(fish_stat,'VariableNames',{'r','p','RMSE','Bias'},...
    'RowNames',{'SAU All Fish','SAU F','SAU P','SAU D','SAU Frac Pelagic',...
    'vanD Frac Pelagic','Stock All Fish'});
writetable(Fstat,[fpath 'FOSI_',mod,'obs_LME_all_ms_stats_' cfile2 '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([fpath 'FOSI_',mod,'obs_LME_all_ms_stats_' cfile2 '.mat'],'fish_stat')


