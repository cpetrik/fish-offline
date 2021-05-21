% P:D ratio by LME 
% CESM FOSI
% Compare to Daniel's model results &SAUP

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6.mat']);
load([cpath 'Data_grid_POP_gx1v6.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

%TAREA units 'cm^2'
AREA_OCN = TAREA * 1e-4;
tlme = lme_mask_esm2m';

%% FEISTY 
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
mod = 'All_fish03';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
dpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

load([dpath 'LME_fosi_fished_',mod,'_' cfile '.mat']);

lme_area_km2 = lme_area * 1e-6;

% FEISTY  LME biomass in MT
plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
% MT/km2
plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;

plme_rPDcatch = plme_Pmcatch ./ (plme_Pmcatch+plme_Dmcatch);

%% DvD on grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/DanielVD_PelDem/Colleen_modeledfish_LME.mat')
dlme_Pfrac = NaN*ones(360,200);
for L=1:63
    lid = find(tlme==L);
    dlme_Pfrac(lid) = FracLP(L);
end

%% SAUP
load([spath 'SAUP_LME_Catch_top10_Stock.mat']);
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')

sFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

l10s=log10(slme_mcatch10+eps);
l10sF=log10(Flme_mcatch10+eps);
l10sP=log10(Plme_mcatch10+eps);
l10sD=log10(Dlme_mcatch10+eps);

%on grid
sFracPD_grid = NaN*ones(360,200);
for L=1:66
    lid = find(tlme==L);
    sFracPD_grid(lid) = sFracPD(L);
end

%% Comparison stats
did=[1:61,63];
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
did2 = notLELC(notLELC<=63);

diffD = rPD_catch - dlme_Pfrac;
diffS = rPD_catch - sFracPD_grid;

%r
rall=corr(FracLP(did),plme_rPDcatch(did));
rall2=corr(FracLP(did2),plme_rPDcatch(did2));
rPD=corr(sFracPD(notLELC),plme_rPDcatch(notLELC));

%root mean square error
o=FracLP(did);
p=plme_rPDcatch(did);
n = length(o);
num=nansum((p-o).^2);
rmse = sqrt(num/n);

o=FracLP(did2);
p=plme_rPDcatch(did2);
n = length(o);
num=nansum((p-o).^2);
rmse2 = sqrt(num/n);

o=sFracPD(notLELC);
p=plme_rPDcatch(notLELC);
n = length(o);
num=nansum((p-o).^2);
rmsePD = sqrt(num/n);

%Fmed
Fall=10^(median(FracLP(did)-plme_rPDcatch(did)));
Fall2=10^(median(FracLP(did2)-plme_rPDcatch(did2)));
FPD=10^(median(sFracPD(notLELC)-plme_rPDcatch(notLELC)));


% Table
fish_stat(1,1) = rall;
fish_stat(2,1) = rmse;
fish_stat(3,1) = Fall;
fish_stat(1,2) = rall2;
fish_stat(2,2) = rmse2;
fish_stat(3,2) = Fall2;
fish_stat(1,3) = rPD;
fish_stat(2,3) = rmsePD;
fish_stat(3,3) = FPD;

Fstat = array2table(fish_stat,'RowNames',{'r','RMSE','Fmed'},...
    'VariableNames',{'DvDAllLMEs','DvDnoLELC','SAUnoLELC'});
writetable(Fstat,[dpath 'Hist_LME_DvD_SAU_stats_' cfile '.csv'],'Delimiter',',','WriteRowNames',true)
save([dpath 'FOSI_LME_DvD_SAU_stats_' cfile '.mat'],'fish_stat')

%% Plot info
[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

cmYOR=cbrewer('seq','YlOrRd',28);
cmRP=cbrewer('seq','RdPu',28);
cmPR=cbrewer('seq','PuRd',28);

x=0:0.1:1;

%% Subplot with maps and corr
figure(1)
%SAU
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffS)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.525 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('FEISTY - SAU difference')

%vanD
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('FEISTY - vanD difference')

%SAU corr
subplot('Position',[0.075 0.075 0.4 0.4])
plot(x,x,'--k');hold on;
%scatter(sFracPD(notLELC),plme_rPDcatch(notLELC),20,'filled'); hold on;
scatter(sFracPD(notLELC),plme_rPDcatch(notLELC),20,lme_ptemp(notLELC,1),'filled'); hold on;
cmocean('thermal');
text(0.75,0.55,['r = ' sprintf('%2.2f',rPD)])
text(0.75,0.5,['RMSE = ' sprintf('%2.2f',rmsePD)])
text(0.75,0.45,['Fmed = ' sprintf('%2.2f',FPD)])
axis([0 1.05 0 1.05])
xlabel('SAU')
ylabel('FEISTY')
%title('Fraction Large Pelagics')

%vanD Corr
subplot('Position',[0.55 0.075 0.4 0.4])
plot(x,x,'--k');hold on;
%scatter(FracLP(did),plme_rPDcatch(did),20,'filled'); hold on;
scatter(FracLP(did),plme_rPDcatch(did),20,lme_ptemp(did,1),'filled'); hold on;
cmocean('thermal');
text(0.75,0.55,['r = ' sprintf('%2.2f',rall)])
text(0.75,0.5,['RMSE = ' sprintf('%2.2f',rmse)])
text(0.75,0.45,['Fmed = ' sprintf('%2.2f',Fall)])
axis([0 1.05 0 1.05])
xlabel('vanD')
ylabel('FEISTY')
%title('Fraction Large Pelagics')
%stamp(cfile)
print('-dpng',[ppath 'FOSI_' harv '_LME_fracPD_catch_SAUP_DvD_comp_subplot_Fmed.png'])

%% Subplot with maps and corr no Fmed
figure(2)
%SAU
subplot('Position',[0 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffS)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.56 0.5 0.025],'orientation','horizontal')
set(gcf,'renderer','painters')
title('FEISTY - SAU difference')

%DvD
subplot('Position',[0.5 0.53 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffD)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('FEISTY - vanD difference')

%SAU corr
subplot('Position',[0.1 0.16 0.35 0.35])
plot(x,x,'--k');hold on;
%scatter(sFracPD(notLELC),plme_rPDcatch(notLELC),20,'filled'); hold on;
scatter(sFracPD(notLELC),plme_rPDcatch(notLELC),20,lme_ptemp(notLELC,1),'filled'); hold on;
cmocean('thermal');
text(0.725,0.55,['r = ' sprintf('%2.2f',rPD)])
text(0.725,0.49,['RMSE = ' sprintf('%2.2f',rmsePD)])
axis([0 1.05 0 1.05])
xlabel('SAU')
ylabel('FEISTY')
%title('Fraction Large Pelagics')

%DvD Corr
subplot('Position',[0.575 0.16 0.35 0.35])
plot(x,x,'--k');hold on;
%scatter(FracLP(did),plme_rPDcatch(did),20,'filled'); hold on;
scatter(FracLP(did),plme_rPDcatch(did),20,lme_ptemp(did,1),'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.25 0.05 0.5 0.025],'orientation','horizontal')
text(0.725,0.55,['r = ' sprintf('%2.2f',rall)])
text(0.725,0.49,['RMSE = ' sprintf('%2.2f',rmse)])
axis([0 1.05 0 1.05])
xlabel('vanD')
ylabel('FEISTY')
%title('Fraction Large Pelagics')
%stamp(cfile)
print('-dpng',[ppath 'FOSI_' harv '_LME_fracPD_catch_SAUP_DvD_comp_subplot.png'])


