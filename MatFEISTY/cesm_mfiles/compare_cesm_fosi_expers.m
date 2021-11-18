% Compare FEISTY runs w/ CESM FOSI experiments
% Time series plots and maps

clear all
close all

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI//';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v14_All_fish03';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% Full
load([fpath 'Plot_Means_FOSI_v14_All_fish03_' cfile '.mat']);

FulltF = F;
FulltP = P;
FulltD = D;
FulltB = B;
FullsF = AllF;
FullsP = AllP;
FullsD = AllD;
FullsS = AllS;
FullsM = AllM;
FullsL = AllL;

FulltA = F+P+D;
FullsA = AllF+AllP+AllD;

clear F P D B AllF AllP AllD AllS AllM AllL

%%
%exper = varFood, varTemp, climatol
load([fpath 'Plot_Means_FOSI_v14_climatol_' cfile '.mat']);

ClimtF = F;
ClimtP = P;
ClimtD = D;
ClimtB = B;
ClimsF = AllF;
ClimsP = AllP;
ClimsD = AllD;
ClimsS = AllS;
ClimsM = AllM;
ClimsL = AllL;

ClimtA = F+P+D;
ClimsA = AllF+AllP+AllD;

clear F P D B AllF AllP AllD AllS AllM AllL

%%
load([fpath 'Plot_Means_FOSI_v14_varTemp_' cfile '.mat']);

TemptF = F;
TemptP = P;
TemptD = D;
TemptB = B;
TempsF = AllF;
TempsP = AllP;
TempsD = AllD;
TempsS = AllS;
TempsM = AllM;
TempsL = AllL;

TemptA = F+P+D;
TempsA = AllF+AllP+AllD;

clear F P D B AllF AllP AllD AllS AllM AllL

%%
load([fpath 'Plot_Means_FOSI_v14_varFood_' cfile '.mat']);

FoodtF = F;
FoodtP = P;
FoodtD = D;
FoodtB = B;
FoodsF = AllF;
FoodsP = AllP;
FoodsD = AllD;
FoodsS = AllS;
FoodsM = AllM;
FoodsL = AllL;

FoodtA = F+P+D;
FoodsA = AllF+AllP+AllD;

clear F P D B AllF AllP AllD AllS AllM AllL

%% Diffs from climatol
dFoodtF = FoodtF - ClimtF;
dFoodtP = FoodtP - ClimtP;
dFoodtD = FoodtD - ClimtD;
dFoodtB = FoodtB - ClimtB;
dFoodtA = FoodtA - ClimtA;
dFoodsF = FoodsF - ClimsF;
dFoodsP = FoodsP - ClimsP;
dFoodsD = FoodsD - ClimsD;
dFoodsS = FoodsS - ClimsS;
dFoodsM = FoodsM - ClimsM;
dFoodsL = FoodsL - ClimsL;
dFoodsA = FoodsA - ClimsA;

dTemptF = TemptF - ClimtF;
dTemptP = TemptP - ClimtP;
dTemptD = TemptD - ClimtD;
dTemptB = TemptB - ClimtB;
dTemptA = TemptA - ClimtA;
dTempsF = TempsF - ClimsF;
dTempsP = TempsP - ClimsP;
dTempsD = TempsD - ClimsD;
dTempsS = TempsS - ClimsS;
dTempsM = TempsM - ClimsM;
dTempsL = TempsL - ClimsL;
dTempsA = TempsA - ClimsA;

dFulltF = FulltF - ClimtF;
dFulltP = FulltP - ClimtP;
dFulltD = FulltD - ClimtD;
dFulltB = FulltB - ClimtB;
dFulltA = FulltA - ClimtA;
dFullsF = FullsF - ClimsF;
dFullsP = FullsP - ClimsP;
dFullsD = FullsD - ClimsD;
dFullsS = FullsS - ClimsS;
dFullsM = FullsM - ClimsM;
dFullsL = FullsL - ClimsL;
dFullsA = FullsA - ClimsA;

%%
[ni,nj]=size(TLONG);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;     

%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black

set(groot,'defaultAxesColorOrder',cm10);

cmBP50=cbrewer('seq','BuPu',50,'PCHIP');

%% Plots in time
t = 1:length(ClimtF); %time;
y = t/12;

% Benthos
figure(1)
plot(y,log10(ClimtB),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(TemptB),'r','Linewidth',2); hold on;
plot(y,log10(FoodtB),'b','Linewidth',2); hold on;
plot(y,log10(FulltB),'k','Linewidth',2); hold on;
legend('clim','temp','prey','full')
legend('location','west')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_comp_exper_ts_B.png'])

%% Forage
figure(2)
plot(y,log10(ClimtF),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(TemptF),'r','Linewidth',2); hold on;
plot(y,log10(FoodtF),'b','Linewidth',2); hold on;
plot(y,log10(FulltF),'k','Linewidth',2); hold on;
legend('clim','temp','prey','full')
legend('location','southwest')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('Forage')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_comp_exper_ts_F.png'])

%% Lg Pel
figure(3)
plot(y,log10(ClimtP),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(TemptP),'r','Linewidth',2); hold on;
plot(y,log10(FoodtP),'b','Linewidth',2); hold on;
plot(y,log10(FulltP),'k','Linewidth',2); hold on;
legend('clim','temp','prey','full')
legend('location','southwest')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('Large pelagic')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_comp_exper_ts_P.png'])

%% Dem
figure(4)
plot(y,log10(ClimtD),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(TemptD),'r','Linewidth',2); hold on;
plot(y,log10(FoodtD),'b','Linewidth',2); hold on;
plot(y,log10(FulltD),'k','Linewidth',2); hold on;
legend('clim','temp','prey','full')
legend('location','northwest')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('Demersal')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_comp_exper_ts_D.png'])

%% All
figure(5)
plot(y,log10(ClimtA),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(TemptA),'r','Linewidth',2); hold on;
plot(y,log10(FoodtA),'b','Linewidth',2); hold on;
plot(y,log10(FulltA),'k','Linewidth',2); hold on;
legend('clim','temp','prey','full')
legend('location','southwest')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('All')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_comp_exper_ts_All.png'])

%% 4 subplots - log10
figure(6)
subplot(2,2,1)
plot(y,log10(ClimtF),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
plot(y,log10(TemptF),'r','Linewidth',1); hold on;
plot(y,log10(FoodtF),'b','Linewidth',1); hold on;
plot(y,log10(FulltF),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('Forage')

% Lg Pel
subplot(2,2,2)
plot(y,log10(ClimtP),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
plot(y,log10(TemptP),'r','Linewidth',1); hold on;
plot(y,log10(FoodtP),'b','Linewidth',1); hold on;
plot(y,log10(FulltP),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('Large pelagic')

% Dem
subplot(2,2,3)
plot(y,log10(ClimtD),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
plot(y,log10(TemptD),'r','Linewidth',1); hold on;
plot(y,log10(FoodtD),'b','Linewidth',1); hold on;
plot(y,log10(FulltD),'k','Linewidth',1); hold on;
legend('clim','temp','prey','full')
legend('location','northwest')
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('Demersal')

% All
subplot(2,2,4)
plot(y,log10(ClimtA),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
plot(y,log10(TemptA),'r','Linewidth',1); hold on;
plot(y,log10(FoodtA),'b','Linewidth',1); hold on;
plot(y,log10(FulltA),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('All')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_comp_exper_ts_All_subplot_log.png'])

%% 4 subplots
figure(7)
subplot(2,2,1)
plot(y,(ClimtF),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
plot(y,(TemptF),'r','Linewidth',1); hold on;
plot(y,(FoodtF),'b','Linewidth',1); hold on;
plot(y,(FulltF),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title('Forage')

% Lg Pel
subplot(2,2,2)
plot(y,(ClimtP),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
plot(y,(TemptP),'r','Linewidth',1); hold on;
plot(y,(FoodtP),'b','Linewidth',1); hold on;
plot(y,(FulltP),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title('Large pelagic')

% Dem
subplot(2,2,3)
plot(y,(ClimtD),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
plot(y,(TemptD),'r','Linewidth',1); hold on;
plot(y,(FoodtD),'b','Linewidth',1); hold on;
plot(y,(FulltD),'k','Linewidth',1); hold on;
legend('clim','temp','prey','full')
legend('location','northwest')
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title('Demersal')

% All
subplot(2,2,4)
plot(y,(ClimtA),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
plot(y,(TemptA),'r','Linewidth',1); hold on;
plot(y,(FoodtA),'b','Linewidth',1); hold on;
plot(y,(FulltA),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title('All')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_comp_exper_ts_All_subplot.png'])

%% 4 subplots
figure(8)
subplot(2,2,1)
%plot(y,(ClimtF),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
plot(y,(dTemptF),'r','Linewidth',1); hold on;
plot(y,(dFoodtF),'b','Linewidth',1); hold on;
plot(y,(dFulltF),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title('Forage')

% Lg Pel
subplot(2,2,2)
%plot(y,(ClimtP),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
plot(y,(dTemptP),'r','Linewidth',1); hold on;
plot(y,(dFoodtP),'b','Linewidth',1); hold on;
plot(y,(dFulltP),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title('Large pelagic')

% Dem
subplot(2,2,3)
%plot(y,(ClimtD),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
plot(y,(dTemptD),'r','Linewidth',1); hold on;
plot(y,(dFoodtD),'b','Linewidth',1); hold on;
plot(y,(dFulltD),'k','Linewidth',1); hold on;
legend('temp','prey','full')
legend('location','northwest')
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title('Demersal')

% All
subplot(2,2,4)
%plot(y,(ClimtA),'color',[0.5 0.5 0.5],'Linewidth',1); hold on;
plot(y,(dTemptA),'r','Linewidth',1); hold on;
plot(y,(dFoodtA),'b','Linewidth',1); hold on;
plot(y,(dFulltA),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g m^-^2)')
title('All')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_comp_exper_ts_diff_subplot.png'])
 
%% Plots in space

% Forage 4 on subplots
figure(9)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(ClimsF))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
cb = colorbar('Position',[0.05 0.5 0.4 0.025],'orientation','horizontal');%,'AxisLocation','in');
title(cb,'log10 g m^-^2')
set(gcf,'renderer','painters')
title('F Climatol')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(dTempsF))
cmocean('balance')   
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
cb2 = colorbar('Position',[0.55 0.5 0.4 0.025],'orientation','horizontal','AxisLocation','in');
title(cb2,'g m^-^2')
set(gcf,'renderer','painters')
title('diff F var Temp')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(dFullsF))
cmocean('balance')   
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('diff F Full')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(dFoodsF))
cmocean('balance')               
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('diff F var Prey')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_global_exper_diffF_subplot.png'])

%% LgPel 4 on subplots
figure(10)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(ClimsP))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
cb = colorbar('Position',[0.05 0.5 0.4 0.025],'orientation','horizontal');%,'AxisLocation','in');
title(cb,'log10 g m^-^2')
set(gcf,'renderer','painters')
title('P Climatol')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(dTempsP))
cmocean('balance')   
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
cb2 = colorbar('Position',[0.55 0.5 0.4 0.025],'orientation','horizontal','AxisLocation','in');
title(cb2,'g m^-^2')
set(gcf,'renderer','painters')
title('diff P var Temp')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(dFullsP))
cmocean('balance')   
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('diff P Full')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(dFoodsP))
cmocean('balance')               
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('diff P var Prey')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_global_exper_diffP_subplot.png'])

%% Dem 4 on subplots
figure(11)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(ClimsD))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
cb = colorbar('Position',[0.05 0.5 0.4 0.025],'orientation','horizontal');%,'AxisLocation','in');
title(cb,'log10 g m^-^2')
set(gcf,'renderer','painters')
title('D Climatol')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(dTempsD))
cmocean('balance')   
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
cb2 = colorbar('Position',[0.55 0.5 0.4 0.025],'orientation','horizontal','AxisLocation','in');
title(cb2,'g m^-^2')
set(gcf,'renderer','painters')
title('diff D var Temp')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(dFullsD))
cmocean('balance')   
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('diff D Full')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(dFoodsD))
cmocean('balance')               
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 2]);
set(gcf,'renderer','painters')
title('diff D var Prey')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_global_exper_diffD_subplot.png'])

%% All 4 on subplots
figure(12)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(ClimsA))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2.5]);
cb = colorbar('Position',[0.05 0.5 0.4 0.025],'orientation','horizontal');%,'AxisLocation','in');
title(cb,'log10 g m^-^2')
set(gcf,'renderer','painters')
title('All Climatol')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(dTempsA))
cmocean('balance')   
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 3]);
cb2 = colorbar('Position',[0.55 0.5 0.4 0.025],'orientation','horizontal','AxisLocation','in');
title(cb2,'g m^-^2')
set(gcf,'renderer','painters')
title('diff All var Temp')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(dFullsA))
cmocean('balance')   
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 3]);
set(gcf,'renderer','painters')
title('diff All Full')

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(dFoodsA))
cmocean('balance')               
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 3]);
set(gcf,'renderer','painters')
title('diff All var Prey')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_global_exper_diffAll_subplot.png'])

%% Bent 4 on subplots
% figure(13)
% % all F
% subplot('Position',[0 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(TLAT,TLONG,log10(ClimsB))
% colormap(cmBP50)
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-2 2]);
% cb = colorbar('Position',[0.05 0.5 0.4 0.025],'orientation','horizontal');%,'AxisLocation','in');
% title(cb,'log10 g m^-^2')
% set(gcf,'renderer','painters')
% title('B Climatol')
% 
% % all D
% subplot('Position',[0 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(TLAT,TLONG,(dTempsB))
% cmocean('balance')   
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-2 2]);
% cb2 = colorbar('Position',[0.55 0.5 0.4 0.025],'orientation','horizontal','AxisLocation','in');
% title(cb2,'g m^-^2')
% set(gcf,'renderer','painters')
% title('diff B var Temp')
% 
% % All P
% subplot('Position',[0.5 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(TLAT,TLONG,(dFullsB))
% cmocean('balance')   
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-2 2]);
% set(gcf,'renderer','painters')
% title('diff B Full')
% 
% % All
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(TLAT,TLONG,(dFoodsB))
% cmocean('balance')               
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-2 2]);
% set(gcf,'renderer','painters')
% title('diff B var Prey')
% stamp('')
% print('-dpng',[ppath 'FOSI_',mod,'_global_exper_diffB_subplot.png'])

