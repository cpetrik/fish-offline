%CESM-FEISTY catch vs. SAUP catch by LME
%Use same methods as Stock et al. 2017 to reduce SAUP dataset

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6.mat']);
load([cpath 'Data_grid_POP_gx1v6.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);
load([fpath 'lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'],'lme_tp_fosi')

%TAREA units 'cm^2'
AREA_OCN = TAREA * 1e-4;

%% use weighted catches
% load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
%     'Flme_wcatch','Dlme_wcatch','Plme_wcatch');

load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
keep = notLELC;

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

%% FEISTY file info
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
mod = 'v12_All_fish03_';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
dpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

load([dpath 'LME_fosi_fished_',mod,cfile '.mat']);
lme_area_km2 = lme_area * 1e-6;

%% plot info
[ni,nj]=size(TLONG);
geolon_t = double(TLONG);
geolat_t = double(TLAT);
plotminlat=-90; 
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; 

load coastlines

% Assign a color to each LME based on temp
% tmap=colormap(jet(66));
% lme_ptemp(:,2)=1:length(lme_ptemp);
% [B,I] = sort(lme_ptemp(:,1));
% I(:,2)=1:length(lme_ptemp);
% [B2,I2] = sort(I(:,1));
% tid = I(I2,:);
% close all

x=-8:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% SAUP
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);
%all P and D
Llme_mcatch10 = Plme_mcatch10 + Dlme_mcatch10; 
l10sL = log10(Llme_mcatch10+eps); 

%% FEISTY LME biomass in MT
plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;

% MT/km2
plme_mcatch = plme_mcatch ./ lme_area_km2;
plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);

l10p=log10(plme_mcatch+eps);
l10pF=log10(plme_Fmcatch+eps);
l10pP=log10(plme_Pmcatch+eps);
l10pD=log10(plme_Dmcatch+eps);
%all P and D
plme_Lmcatch = plme_Pmcatch + plme_Dmcatch; 
l10pL = log10(plme_Lmcatch+eps);

%% on grid
tlme = double(lme_mask);
tlme(tlme<0) = nan;
pFracPD_grid = NaN*ones(ni,nj);
sFracPD_grid = NaN*ones(ni,nj);
l10sF_grid = NaN*ones(ni,nj);
l10pF_grid = NaN*ones(ni,nj);
l10sP_grid = NaN*ones(ni,nj);
l10pP_grid = NaN*ones(ni,nj);
for L=1:66
    lid = find(tlme==L);
    pFracPD_grid(lid) = pFracPD(L);
    sFracPD_grid(lid) = sFracPD(L);
%     l10sF_grid(lid) = l10sF(L);
%     l10pF_grid(lid) = l10pF(L);
%     l10sP_grid(lid) = l10sP(L);
%     l10pP_grid(lid) = l10pP(L);
end

for L=1:length(keep)
    lme=keep(L);
    lid = find(tlme==lme);
    l10sF_grid(lid) = l10sF(lme);
    l10pF_grid(lid) = l10pF(lme);
    l10sP_grid(lid) = l10sP(lme);
    l10pP_grid(lid) = l10pP(lme);
end

diffPD = pFracPD_grid - sFracPD_grid;
diffF = (l10pF_grid - l10sF_grid);
diffP = (l10pP_grid - l10sP_grid);

%% Drop Arctic, Antarctic, Hawaii, Australia -------------------------
% Stats
%r
[rall,pall]=corr(l10s(keep),l10p(keep));
[rF,pF]=corr(l10sF(keep),l10pF(keep));
[rP,pP]=corr(l10sP(keep),l10pP(keep));
[rD,pD]=corr(l10sD(keep),l10pD(keep));
[rL,pL]=corr(l10sL(keep),l10pL(keep));
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

o=l10sL(keep);
p=l10pL(keep);
n = length(o);
num=nansum((p-o).^2);
rmseL = sqrt(num/n);

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
FL=10^(median(l10sL(keep)-l10pL(keep)));
FPD=10^(median(sFracPD(keep)-pFracPD(keep)));

% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rF;
fish_stat(3,1) = rP;
fish_stat(4,1) = rD;
fish_stat(5,1) = rL;
fish_stat(6,1) = rPD;
fish_stat(1,2) = rmse;
fish_stat(2,2) = rmseF;
fish_stat(3,2) = rmseP;
fish_stat(4,2) = rmseD;
fish_stat(5,2) = rmseL;
fish_stat(6,2) = rmsePD;
fish_stat(1,3) = Fall;
fish_stat(2,3) = FF;
fish_stat(3,3) = FP;
fish_stat(4,3) = FD;
fish_stat(5,3) = FL;
fish_stat(6,3) = FPD;

Fstat = array2table(fish_stat,'VariableNames',{'r','RMSE','Fmed'},...
    'RowNames',{'All','F','P','D','L','FracP'});
writetable(Fstat,[dpath 'LME_SAUP_stats_' mod cfile '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([dpath 'LME_SAUP_stats_' mod cfile '.mat'],'fish_stat')

%% Plots - same as ms
figure(1)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
% scatter(l10sF(keep),l10pF(keep),20,'k','filled'); hold on;
scatter(l10sF(keep),l10pF(keep),20,lme_tp_fosi(keep,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.375 0.5 0.3 0.025],'orientation','horizontal')
text(-5.5,1.5,['r = ' sprintf('%2.2f',rF) ' (p = ' sprintf('%2.2f',pF) ')'])
text(-5.5,1.0,['RMSE = ' sprintf('%2.2f',rmseF)])
axis([-6 2 -6 2])
xlabel('SAU')
ylabel('FEISTY ')
title('Forage Fishes')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
% scatter(l10sP(keep),l10pP(keep),20,'k','filled'); hold on;
scatter(l10sP(keep),l10pP(keep),20,lme_tp_fosi(keep,1),'filled'); hold on;
cmocean('thermal');
text(-4.5,1.5,['r = ' sprintf('%2.2f',rP) ' (p = ' sprintf('%2.2f',pP) ')'])
text(-4.5,1.0,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-5 2 -5 2])
xlabel('SAU')
ylabel('FEISTY ')
title('Large Pelagics')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
% scatter(l10sD(keep),l10pD(keep),20,'k','filled'); hold on;
scatter(l10sD(keep),l10pD(keep),20,lme_tp_fosi(keep,1),'filled'); hold on;
cmocean('thermal');
text(-1.75,1.7,['r = ' sprintf('%2.2f',rD) ' (p = ' sprintf('%2.2f',pD) ')'])
text(-1.75,1.4,['RMSE = ' sprintf('%2.2f',rmseD)])
axis([-2 2 -2 2])
xlabel('SAU')
ylabel('FEISTY ')
title('Demersals')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
% scatter(l10s(keep),l10p(keep),20,'k','filled'); hold on;
scatter(l10s(keep),l10p(keep),20,lme_tp_fosi(keep,1),'filled'); hold on;
cmocean('thermal');
text(-1.75,1.7,['r = ' sprintf('%2.2f',rall) ' (p = ' sprintf('%2.2f',pall) ')'])
text(-1.75,1.4,['RMSE = ' sprintf('%2.2f',rmse)])
axis([-2 2 -2 2])
xlabel('SAU')
ylabel('FEISTY ')
title('All fishes')
print('-dpng',[ppath 'FOSI_',mod,'_SAUP_comp_types_temp_Stock_LELC.png'])

%% 4plot - no Forage
figure(10)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10sP(keep),l10pP(keep),20,lme_tp_fosi(keep,1),'filled'); hold on;
cmocean('thermal');
text(-4.5,1.5,['r = ' sprintf('%2.2f',rP) ' (p = ' sprintf('%2.2f',pP) ')'])
text(-4.5,1.0,['RMSE = ' sprintf('%2.2f',rmseP)])
axis([-5 2 -5 2])
xlabel('SAU')
ylabel('FEISTY ')
title('Large Pelagics')

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10sD(keep),l10pD(keep),20,lme_tp_fosi(keep,1),'filled'); hold on;
cmocean('thermal');
text(-1.75,1.7,['r = ' sprintf('%2.2f',rD) ' (p = ' sprintf('%2.2f',pD) ')'])
text(-1.75,1.4,['RMSE = ' sprintf('%2.2f',rmseD)])
axis([-2 2 -2 2])
xlabel('SAU')
ylabel('FEISTY ')
title('Demersals')

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
scatter(l10sL(keep),l10pL(keep),20,lme_tp_fosi(keep,1),'filled'); hold on;
cmocean('thermal');
text(-1.75,1.7,['r = ' sprintf('%2.2f',rL) ' (p = ' sprintf('%2.2f',pL) ')'])
text(-1.75,1.4,['RMSE = ' sprintf('%2.2f',rmseL)])
axis([-2 2 -2 2])
xlabel('SAU')
ylabel('FEISTY ')
title('Large pelagics + Demersals')

subplot(2,2,4)
plot(x,x,'--k'); hold on;
% scatter(sFracPD(keep),pFracPD(keep),20,'k','filled'); hold on;
scatter(sFracPD(keep),pFracPD(keep),20,lme_tp_fosi(keep,1),'filled'); hold on;
cmocean('thermal');
axis([0 1.05 0 1.05])
text(0.05,0.99,['r = ' sprintf('%2.2f',rPD) ' (p = ' sprintf('%2.2f',pPD) ')'])
text(0.05,0.93,['RMSE = ' sprintf('%2.2f',rmsePD)])
xlabel('SAU')
ylabel('FEISTY ')
title('Large Pelagics vs Demersals')
print('-dpng',[ppath 'FOSI_',mod,'_SAUP_comp_types_temp_Stock_LELC_noF.png'])

%% P outliers
figure(2)
plot(x,x,'--k'); hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
% scatter(l10sP(keep),l10pP(keep),20,'k','filled'); hold on;
scatter(l10sP(keep),l10pP(keep),20,lme_tp_fosi(keep,1),'filled'); hold on;
cmocean('thermal');
for i=1:length(keep)
    if (l10pP(keep(i))<-3)
        text(l10sP(keep(i)),l10pP(keep(i)),num2str(keep(i)))
    end
end
axis([-3 2 -65 2])
xlabel('SAU')
ylabel('FEISTY ')
title('Large Pelagics')
print('-dpng',[ppath 'FOSI_',mod,'_SAUP_comp_P_Stock_LELC_temp.png'])

%% without outliers
keep2 = keep([1:18,21:42,45]);
[rP2,pP2]=corr(l10sP(keep2),l10pP(keep2));

%root mean square error
o=l10sP(keep2);
p=l10pP(keep2);
n = length(o);
num=nansum((p-o).^2);
rmseP2 = sqrt(num/n);

%% P:D ratio
figure(3)
plot(x,x,'--k'); hold on;
% scatter(sFracPD(keep),pFracPD(keep),20,'k','filled'); hold on;
scatter(sFracPD(keep),pFracPD(keep),20,lme_tp_fosi(keep,1),'filled'); hold on;
cmocean('thermal');
axis([0 1 0 1])
text(0.05,0.96,['r = ' sprintf('%2.2f',rPD) ' (p = ' sprintf('%2.2f',pPD) ')'])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmsePD)])
xlabel('SAU')
ylabel('FEISTY ')
title('Large Pelagics vs Demersals')
%print('-dpng',[ppath 'FOSI_',mod,'_SAUP_comp_PDratio_Stock_LELC_temp.png'])

