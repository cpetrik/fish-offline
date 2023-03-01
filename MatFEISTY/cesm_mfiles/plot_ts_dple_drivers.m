% Visualize input forcing of FEISTY
% CESM DPLE
% Time series plots and maps

clear
close all

%% FOSI
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

ID = GRD.ID;

tlme = double(lme_mask);
tlme(tlme<0) = nan;
lme_grid = tlme(ID);

%% DPLE
fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/DPLE/';

% GET DPLE time info ----------------------------
ncdisp([fpath 'DPLE-FIESTY-forcing_zooC_150m.nc'])

% Pelagic temperature
ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_TEMP_150m.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end

Y = double(Y);

%%
[ni,nj] = size(TLONG);
mos = length(L);

% FOSI is 1948-2015
% DPLE is 1954-2017
% leadtimes = 1:24;
firstyear = 1954; %of DPLE initializations
lastyear  = 2015; %of DPLE initializations
startmonth = 11;  %of DPLE initializations

%% Select year and member
yr = 1; %1:length(Y)

load([fpath 'Time_Means_DPLE_LME_drivers_Y',num2str(Y(yr)),'.mat']);

%% Units
% meso zoo: nmolC cm-2 to g(WW) m-2
% 1e9 nmol in 1 mol C
% 1e4 cm2 in 1 m2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
tMZ = tMZ * 1e-9 * 1e4 * 12.01 * 9.0;
lme_mMZ = lme_mMZ * 1e-9 * 1e4 * 12.01 * 9.0;

% meso zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
tMZl = tMZl * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
lme_mMZl = lme_mMZl * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

% detrital flux to benthos: nmolC cm-2 s-1 to g(WW) m-2 d-1
tDet = tDet * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
lme_mDet = lme_mDet * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

%% Plots
y=(Y(yr)+(11/12)):(1/12):(Y(yr)+11);

mmz = max(tMZ');
%mid = find(~isinf(mmz));
mid = find(mmz<100);

md = max(tDet');
did = ~isinf(md);

%% global
figure(1)
subplot(2,2,1)
plot(y,tTP,'r','Linewidth',0.5); hold on;
plot(y,nanmean(tTP,1),'color',[0.75 0 0],'Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Temperature (^oC)')
title('T Pelagic')

subplot(2,2,2)
plot(y,tTB,'b','Linewidth',0.5); hold on;
plot(y,nanmean(tTB,1),'color',[0 0 0.65],'Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Temperature (^oC)')
title('T Benthic')

subplot(2,2,3)
plot(y,tMZ(mid,:),'color',[0 0.75 0.5],'Linewidth',0.5); hold on;
plot(y,nanmean(tMZ(mid,:),1),'color',[0 0.5 0.25],'Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('biomass (g m^-^2)')
title('Mesozooplankton')

subplot(2,2,4)
plot(y,tDet,'color',[0.6 0.6 0.6],'Linewidth',0.5); hold on;
plot(y,nanmean(tDet,1),'color',[0.3 0.3 0.3],'Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('flux (g m^-^2 d^-^1)')
title('Detritus')

%print('-dpng',[ppath 'DPLE_v14_Y2015_All_fish03_CCLME_all_types.png'])

%% CCE
figure(2)
subplot(2,2,1)
plot(y,squeeze(lme_mTP(3,:,:)),'r','Linewidth',0.5); hold on;
plot(y,nanmean(squeeze(lme_mTP(3,:,:)),2),'color',[0.75 0 0],'Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Temperature (^oC)')
title('T Pelagic')

subplot(2,2,2)
plot(y,squeeze(lme_mTB(3,:,:)),'b','Linewidth',0.5); hold on;
plot(y,nanmean(squeeze(lme_mTB(3,:,:)),2),'color',[0 0 0.65],'Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Temperature (^oC)')
title('T Benthic')

subplot(2,2,3)
plot(y,squeeze(lme_mMZ(3,:,:)),'color',[0 0.75 0.5],'Linewidth',0.5); hold on;
plot(y,nanmean(squeeze(lme_mMZ(3,:,:)),2),'color',[0 0.5 0.25],'Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('biomass (g m^-^2)')
title('Mesozooplankton')

subplot(2,2,4)
plot(y,squeeze(lme_mDet(3,:,:)),'color',[0.6 0.6 0.6],'Linewidth',0.5); hold on;
plot(y,nanmean(squeeze(lme_mDet(3,:,:)),2),'color',[0.3 0.3 0.3],'Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('flux (g m^-^2 d^-^1)')
title('Detritus')

%print('-dpng',[ppath 'DPLE_v14_Y2015_All_fish03_CCLME_all_types.png'])



