% Visualize output of FEISTY
% CESM 4P4Z
% Time series plots and maps

clear all
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_Sm100_nmort1_BE08_noCC_RE00100';
mod = 'Spinup_All_fish03';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/4P4Z/';
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
% load([fpath 'Time_Means_4P4Z_' cfile '.mat'],'mz_tmfrac','mz_ttf',...
%     'lz_tmfrac','lz_ttf');
% load([fpath 'Space_Means_4P4Z_' cfile '.mat'],'mz_smfrac','mz_stf',...
%     'lz_smfrac','lz_stf');
% load([fpath 'Annual_Means_4P4Z_' cfile '.mat'],'mz_mtf','lz_mtf');
load([fpath 'Means_4P4Z_Spinup_' cfile '.mat']);

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/4P4Z/';
load([cpath 'gridspec_POP_gx1v6_4p4z.mat']);
load([cpath 'Data_grid_POP_gx1v6_4p4z.mat']);

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

%% 
t = 1:length(mz_tmfrac); %time;
y = t/12;

% MZ loss plots
nx = length(mz_stf);
nt = length(mz_ttf);
%Mean fraction
Cmz_smfrac = mz_smfrac;
Clz_smfrac = lz_smfrac;
%Total times it happens over time
Cmz_ttover = mz_ttf/nx;
Clz_ttover = lz_ttf/nx;
%Total times it happens in each year in space?
Cmz_stover5 = mz_stf/nt;
Clz_stover5 = lz_stf/nt;

%%
%happens whole year
test1=floor(Clz_stover5);
sum(test1)/length(test1) % = 7.90

%happens whole year
test2=floor(Cmz_stover5);
%histogram(test2)
sum(test2)/length(test2) % = 0.02

test3=round(Clz_stover5);
sum(test3)/length(test3) % = 8.38

%happens >=50% of year
test4=round(Cmz_stover5);
%histogram(test4)
sum(test4)/length(test4) % = 0.03

%% Plot in time
figure(6)
plot(y, Cmz_ttover,'k','LineWidth',2); hold on;
plot(y, Clz_ttover,'b','LineWidth',2); hold on;
legend('Z3','Z4')
xlabel('Years')
ylabel('Fraction of grid cells over-consumed')
print('-dpng',[ppath '4P4Z_',mod '_timeseries_zoop_overcon.png'])

%% Plots in space
CFmz=NaN*ones(ni,nj);
COmz=NaN*ones(ni,nj);
COmz5=NaN*ones(ni,nj);
CFmz(GRD.ID)=Cmz_smfrac;
COmz(GRD.ID)=Cmz_stover5;

CFlz=NaN*ones(ni,nj);
COlz=NaN*ones(ni,nj);
COlz5=NaN*ones(ni,nj);
CFlz(GRD.ID)=Clz_smfrac;
COlz(GRD.ID)=Clz_stover5;

% save([bpath '4P4Z_' mod '_ts_map_zoop_overcon.mat'],...
%     'Cmz_ttover','CFmz','COmz','-append');

cmBP=cbrewer('seq','BuPu',10,'PCHIP');
cmBP2 = cmBP;
cmBP2(11,:) = [0 0 0];

%% 2 plots of both maps
figure(7)
%1 - m frac
subplot('Position',[0.01 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,CFmz)
colormap(cmBP2)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1.1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean fraction Z3 hploss consumed','HorizontalAlignment','center')

%2 - m over
subplot('Position',[0.01 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,COmz)
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean times Z3 overconsumed','HorizontalAlignment','center')


%3 - l frac
subplot('Position',[0.5 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,CFlz)
colormap(cmBP2)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1.1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean fraction Z4 hploss consumed','HorizontalAlignment','center')

%4 - l over
subplot('Position',[0.5 0.01 0.45 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,COlz)
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'Mean times Z4 overconsumed','HorizontalAlignment','center')
print('-dpng',[ppath '4P4Z_',mod '_global_zoop_overcon.png'])

%%
figure(8)
%1 - m frac
subplot('Position',[0.01 0.5 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(CFmz))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([-2 12]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'log_1_0 Mean fraction Z3 hploss consumed','HorizontalAlignment','center')

subplot('Position',[0.01 0.01 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(CFlz))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 12]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'log_1_0 Mean fraction Z4 hploss consumed','HorizontalAlignment','center')
print('-dpng',[ppath '4P4Z_',mod '_global_zoop_fraccon.png'])

%%
figure(9)
%1 - m frac
subplot('Position',[0.01 0.5 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(COmz*nt))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 24]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'total times Z3 hploss overconsumed','HorizontalAlignment','center')

subplot('Position',[0.01 0.01 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(COlz*nt))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 700]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'total times Z4 hploss overconsumed','HorizontalAlignment','center')
print('-dpng',[ppath '4P4Z_',mod '_global_zoop_months_overcon.png'])

%%
figure(10)
%1 - m frac
subplot('Position',[0.01 0.5 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(COmz*nt))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 10]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'total times Z3 hploss overconsumed','HorizontalAlignment','center')

subplot('Position',[0.01 0.01 0.9 0.45])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(COlz*nt))
colormap(cmBP)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 10]);
colorbar
set(gcf,'renderer','painters')
text(0,1.6,'total times Z4 hploss overconsumed','HorizontalAlignment','center')
%print('-dpng',[ppath '4P4Z_',mod '_global_zoop_fraccon.png'])
