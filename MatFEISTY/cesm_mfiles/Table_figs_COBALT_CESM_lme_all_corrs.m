% All correlations between COBALT-Hist and CESM-FOSI

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
%keep = notLELC;
keep=1:66;

%% COBALT Hindcast grid
hpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([hpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t','AREA_OCN'); 
grid = csvread([hpath 'grid_csv.csv']); 
[hi,hj]=size(geolon_t);

load([hpath 'lme_mask_esm2m.mat']);
load([bpath 'LME_hist9095_temp_zoop_det.mat'],'lme_ptemp','lme_area');

hlme_ptemp = lme_ptemp;
clear lme_ptemp

hlme_area_km2 = lme_area * 1e-6;
clear lme_area

hlme = lme_mask_esm2m';
clear lme_mask_esm2m

hID = grid(:,1);

%% CESM FOSI grid
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6.mat']);
load([cpath 'Data_grid_POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

load([cpath 'LME-mask-POP_gx1v6.mat']);
load([cpath 'lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'],...
    'lme_tp_fosi','lme_area');

clme_ptemp = lme_tp_fosi;

clme_area_km2 = lme_area * 1e-6;
clear lme_area

clme = lme_mask;

%% FEISTY Output
cfile2 = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
mod = 'All_fish03';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile2 '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% COBALT
cfile1 = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
gpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile1 '/Historic_ESM2M/'];
load([gpath 'LME_hist_',harv,'_' cfile1 '.mat'],'lme_mcatch',...
    'lme_mbio','lme_sbio');

% hlme_mcatch = lme_mcatch;
% hlme_mbio = lme_mbio;
% hlme_sbio = lme_sbio;

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

hFracPD = hlme_Pmcatch ./ (hlme_Pmcatch + hlme_Dmcatch);
hFracPF = hlme_Pmcatch ./ (hlme_Pmcatch + hlme_Fmcatch);
hFracLM = hlme_Lmbio ./ (hlme_Lmbio + hlme_Mmbio);

l10h=log10(hlme_mcatch);
l10hF=log10(hlme_Fmcatch);
l10hP=log10(hlme_Pmcatch);
l10hD=log10(hlme_Dmcatch);
l10hB=log10(hlme_Bmbio);

% l10hATL = log10(lme_te(:,2));
% l10hHTL = log10(lme_te(:,4));
% l10hLTL = log10(lme_te(:,6));

clear lme_mcatch lme_mbio lme_sbio

%% CESM
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile2 '/'];
load([fpath 'LME_fosi_fished_',mod,'_' cfile2 '.mat'],'lme_mcatch',...
    'lme_mbio','lme_sbio');

% clme_mcatch = lme_mcatch;
% clme_mbio = lme_mbio;
% clme_sbio = lme_sbio;

clme_mcatch = nansum(lme_mcatch,2) * 1e-6;
clme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
clme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
clme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
clme_Bmbio = lme_mbio(:,9) * 1e-6;
clme_Mmbio = (lme_mbio(:,4) + lme_mbio(:,5) + lme_mbio(:,6)) * 1e-6;
clme_Lmbio = (lme_mbio(:,7) + lme_mbio(:,8)) * 1e-6;
% MT/km2
clme_mcatch = clme_mcatch ./ clme_area_km2;
clme_Fmcatch = clme_Fmcatch ./ clme_area_km2;
clme_Pmcatch = clme_Pmcatch ./ clme_area_km2;
clme_Dmcatch = clme_Dmcatch ./ clme_area_km2;
clme_Bmbio = clme_Bmbio ./ clme_area_km2;
clme_Mmbio = clme_Mmbio ./ clme_area_km2;
clme_Lmbio = clme_Lmbio ./ clme_area_km2;

cFracPD = clme_Pmcatch ./ (clme_Pmcatch + clme_Dmcatch);
cFracPF = clme_Pmcatch ./ (clme_Pmcatch + clme_Fmcatch);
cFracLM = clme_Lmbio ./ (clme_Lmbio + clme_Mmbio);

l10c=log10(clme_mcatch);
l10cF=log10(clme_Fmcatch);
l10cP=log10(clme_Pmcatch);
l10cD=log10(clme_Dmcatch);
l10cB=log10(clme_Bmbio);

% FEISTY LME TEeffs
%     lme_te(L,2) = nanmean(TEeff_L(lid));
%     lme_te(L,4) = nanmean(TEeff_HTLd(lid));
%     lme_te(L,6) = nanmean(TEeff_LTLd(lid));
% l10cATL = log10(lme_te(:,2));
% l10cHTL = log10(lme_te(:,4));
% l10cLTL = log10(lme_te(:,6));

clear lme_mcatch lme_mbio lme_sbio %lme_te


%% Stats
%r %  ALL CORRS LOOK VERY WRONG
[rall,pall]=corr(l10h(keep),l10c(keep));
[rF,pF]=corr(l10hF(keep),l10cF(keep));
[rP,pP]=corr(l10hP(keep),l10cP(keep));
[rD,pD]=corr(l10hD(keep),l10cD(keep));
[rB,pB]=corr(l10hB(keep),l10cB(keep));
[rPD,pPD]=corr(hFracPD(keep),cFracPD(keep));
[rPF,pPF]=corr(hFracPF(keep),cFracPF(keep));
[rLM,pLM]=corr(hFracLM(keep),cFracLM(keep));
% [rATL,pATL]=corr(l10hATL(keep),l10cATL(keep));
% [rHTL,pHTL]=corr(l10hHTL(keep),l10cHTL(keep));
% [rLTL,pLTL]=corr(l10hLTL(keep),l10cLTL(keep));

%root mean square error
o=l10h(keep);
p=l10c(keep);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=l10hF(keep);
p=l10cF(keep);
n = length(o);
num=nansum((p-o).^2);
rmseF = sqrt(num/n);

o=l10hP(keep);
p=l10cP(keep);
n = length(o);
num=nansum((p-o).^2);
rmseP = sqrt(num/n);

o=l10hD(keep);
p=l10cD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseD = sqrt(num/n);

o=l10hB(keep);
p=l10cB(keep);
n = length(o);
num=nansum((p-o).^2);
rmseB = sqrt(num/n);

o=hFracPD(keep);
p=cFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmsePD = sqrt(num/n);

o=hFracPF(keep);
p=cFracPF(keep);
n = length(o);
num=nansum((p-o).^2);
rmsePF = sqrt(num/n);

o=hFracLM(keep);
p=cFracLM(keep);
n = length(o);
num=nansum((p-o).^2);
rmseLM = sqrt(num/n);

% o=l10hATL(keep);
% p=l10cATL(keep);
% n = length(o);
% num=nansum((p-o).^2);
% rmseATL = sqrt(num/n);
% 
% o=l10hHTL(keep);
% p=l10cHTL(keep);
% n = length(o);
% num=nansum((p-o).^2);
% rmseHTL = sqrt(num/n);
% 
% o=l10hLTL(keep);
% p=l10cLTL(keep);
% n = length(o);
% num=nansum((p-o).^2);
% rmseLTL = sqrt(num/n);

% Bias (Historic minus FOSI)
%average error = bias
%skill(3,c) = nansum(p-o) / n;
o=l10h(keep);
p=l10c(keep);
n = length(o);
bias = nansum(o-p) / n;

o=l10hF(keep);
p=l10cF(keep);
n = length(o);
biasF = nansum(o-p) / n;

o=l10hP(keep);
p=l10cP(keep);
n = length(o);
biasP = nansum(o-p) / n;

o=l10hD(keep);
p=l10cD(keep);
n = length(o);
biasD = nansum(o-p) / n;

o=l10hB(keep);
p=l10cB(keep);
n = length(o);
biasB = nansum(o-p) / n;

o=hFracPD(keep);
p=cFracPD(keep);
n = length(o);
biasPD = nansum(o-p) / n;

o=hFracPF(keep);
p=cFracPF(keep);
n = length(o);
biasPF = nansum(o-p) / n;

o=hFracLM(keep);
p=cFracLM(keep);
n = length(o);
biasLM = nansum(o-p) / n;

% o=l10hATL(keep);
% p=l10cATL(keep);
% n = length(o);
% biasATL = nansum(o-p) / n;
% 
% o=l10hHTL(keep);
% p=l10cHTL(keep);
% n = length(o);
% biasHTL = nansum(o-p) / n;
% 
% o=l10hLTL(keep);
% p=l10cLTL(keep);
% n = length(o);
% biasLTL = nansum(o-p) / n;

%% Table
fish_stat(1,1) = rall; %  ALL CORRS LOOK VERY WRONG - could be very low biomasses with FOSI
fish_stat(2,1) = rF;
fish_stat(3,1) = rP;
fish_stat(4,1) = rD;
fish_stat(5,1) = rB;
fish_stat(6,1) = rPD;
fish_stat(7,1) = rPF;
fish_stat(8,1) = rLM;
% fish_stat(9,1) = rATL;
% fish_stat(10,1) = rHTL;
% fish_stat(11,1) = rLTL;

fish_stat(1,2) = pall;
fish_stat(2,2) = pF;
fish_stat(3,2) = pP;
fish_stat(4,2) = pD;
fish_stat(5,2) = pB;
fish_stat(6,2) = pPD;
fish_stat(7,2) = pPF;
fish_stat(8,2) = pLM;
% fish_stat(9,2) = pATL;
% fish_stat(10,2) = pHTL;
% fish_stat(11,2) = pLTL;

fish_stat(1,3) = rmse;
fish_stat(2,3) = rmseF;
fish_stat(3,3) = rmseP;
fish_stat(4,3) = rmseD;
fish_stat(5,3) = rmseB;
fish_stat(6,3) = rmsePD;
fish_stat(7,3) = rmsePF;
fish_stat(8,3) = rmseLM;
% fish_stat(9,3) = rmseATL;
% fish_stat(10,3) = rmseHTL;
% fish_stat(11,3) = rmseLTL;

fish_stat(1,4) = bias;
fish_stat(2,4) = biasF;
fish_stat(3,4) = biasP;
fish_stat(4,4) = biasD;
fish_stat(5,4) = biasB;
fish_stat(6,4) = biasPD;
fish_stat(7,4) = biasPF;
fish_stat(8,4) = biasLM;
% fish_stat(9,4) = biasATL;
% fish_stat(10,4) = biasHTL;
% fish_stat(11,4) = biasLTL;

% save
Fstat = array2table(fish_stat,'VariableNames',{'r','p','RMSE','Bias'},...
    'RowNames',{'All Fish','F','P','D','B','Frac Pel-Dem','Frac Pel-Forage',...
    'Frac Large-Med'});%,'TEeffL','TEeffHTL','TEeffLTL'});
writetable(Fstat,[fpath 'LME_fosi_hist_stats_' cfile2 '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([fpath 'LME_fosi_hist_stats_' cfile2 '.mat'],'fish_stat')

%% Figures
x=-8:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

figure(1)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hF(keep),l10cF(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
% text(-2.75,0.75,['r = ' sprintf('%2.2f',rF) ' (p = ' sprintf('%2.2f',pF) ')'])
% text(-2.75,0.5,['RMSE = ' sprintf('%2.2f',rmseF)])
axis([-3 1 -3 1])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('Forage Fishes')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hP(keep),l10cP(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
% text(-5.5,0.5,['r = ' sprintf('%2.2f',rP) ' (p = ' sprintf('%2.2f',pP) ')'])
% text(-5.5,0.1,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-6 1 -6 1])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('Large Pelagics')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hD(keep),l10cD(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
% text(-1.75,0.8,['r = ' sprintf('%2.2f',rD) ' (p = ' sprintf('%2.2f',pD) ')'])
% text(-1.75,0.5,['RMSE = ' sprintf('%2.2f',rmseD)])
axis([-2 1 -2 1])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('Demersals')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10h(keep),l10c(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
% text(-1.75,1.4,['r = ' sprintf('%2.2f',rall) ' (p = ' sprintf('%2.2f',pall) ')'])
% text(-1.75,1.1,['RMSE = ' sprintf('%2.2f',rmse)])
axis([-2 1.5 -2 1.5])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('All fishes')
stamp([harv '_' cfile2])
print('-dpng',[ppath 'FOSI_Hist_',harv,'_comp_types_temp.png'])

%%
figure(10)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
scatter(l10hF(keep),l10cF(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
% text(-2.75,0.75,['r = ' sprintf('%2.2f',rF) ' (p = ' sprintf('%2.2f',pF) ')'])
% text(-2.75,0.5,['RMSE = ' sprintf('%2.2f',rmseF)])
axis([-3 1 -300 1])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('Forage Fishes')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
scatter(l10hP(keep),l10cP(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
% text(-5.5,0.5,['r = ' sprintf('%2.2f',rP) ' (p = ' sprintf('%2.2f',pP) ')'])
% text(-5.5,0.1,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-15 1 -300 1])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('Large Pelagics')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
scatter(l10hD(keep),l10cD(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
% text(-1.75,0.8,['r = ' sprintf('%2.2f',rD) ' (p = ' sprintf('%2.2f',pD) ')'])
% text(-1.75,0.5,['RMSE = ' sprintf('%2.2f',rmseD)])
axis([-2 1 -300 1])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('Demersals')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
scatter(l10h(keep),l10c(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
% text(-1.75,1.4,['r = ' sprintf('%2.2f',rall) ' (p = ' sprintf('%2.2f',pall) ')'])
% text(-1.75,1.1,['RMSE = ' sprintf('%2.2f',rmse)])
axis([-2 1.5 -300 1.5])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('All fishes')
stamp([harv '_' cfile2])
print('-dpng',[ppath 'FOSI_Hist_',harv,'_comp_types_temp_outliers.png'])

%% Fractions
figure(2)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(hFracPF(keep),cFracPF(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
text(-0.1,1.3,['r = ' sprintf('%2.2f',rPF) ' (p = ' sprintf('%2.2f',pPF) ')'])
text(-0.1,1.2,['RMSE = ' sprintf('%2.2f',rmsePF)])
axis([-0.2 1.4 -0.2 1.4])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('P / (P+F)')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(hFracPD(keep),cFracPD(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-0.1,1.3,['r = ' sprintf('%2.2f',rPD) ' (p = ' sprintf('%2.2f',pPD) ')'])
text(-0.1,1.2,['RMSE = ' sprintf('%2.2f',rmsePD)])
axis([-0.2 1.4 -0.2 1.4])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('P / (P+D)')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(hFracLM(keep),cFracLM(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-0.1,1.3,['r = ' sprintf('%2.2f',rLM) ' (p = ' sprintf('%2.2f',pLM) ')'])
text(-0.1,1.2,['RMSE = ' sprintf('%2.2f',rmseLM)])
axis([-0.2 1.4 -0.2 1.4])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('L / (L+M)')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hB(keep),l10cB(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
% text(-3.75,-0.3,['r = ' sprintf('%2.2f',rB) ' (p = ' sprintf('%2.2f',pB) ')'])
% text(-3.75,-0.7,['RMSE = ' sprintf('%2.2f',rmseB)])
axis([-4 0 -4 0])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('Benthos')
stamp([harv '_' cfile2])
print('-dpng',[ppath 'FOSI_Hist_',harv,'_comp_fracs_temp.png'])


%% TEeffs
figure(3)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hATL(keep),l10cATL(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
text(-4.75,-1.25,['r = ' sprintf('%2.2f',rATL) ' (p = ' sprintf('%2.2f',pATL) ')'])
text(-4.75,-1.5,['RMSE = ' sprintf('%2.2f',rmseATL)])
axis([-5 -1 -5 -1])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('log_1_0 TEeff ATL')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hHTL(keep),l10cHTL(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-3.25,-0.75,['r = ' sprintf('%2.2f',rHTL) ' (p = ' sprintf('%2.2f',pHTL) ')'])
text(-3.25,-1.0,['RMSE = ' sprintf('%2.2f',rmseHTL)])
axis([-3.5 -0.5 -3.5 -0.5])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('log_1_0 TEeff HTL')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10hLTL(keep),l10cLTL(keep),20,clme_ptemp(keep,1),'filled'); hold on;
cmocean('thermal');
text(-2.4,-0.65,['r = ' sprintf('%2.2f',rLTL) ' (p = ' sprintf('%2.2f',pLTL) ')'])
text(-2.4,-0.8,['RMSE = ' sprintf('%2.2f',rmseLTL)])
axis([-2.5 -0.5 -2.5 -0.5])
xlabel('Hist COBALT')
ylabel('FOSI CESM')
title('log_1_0 TEeff LTL')

stamp([harv '_' cfile])
print('-dpng',[ppath 'FOSI_Hist_',harv,'_comp_TEeffs_temp.png'])

