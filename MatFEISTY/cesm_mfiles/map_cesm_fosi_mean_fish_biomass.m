% Visualize output of FEISTY
% CESM FOSI
% Time series plots and maps

clear 
close all

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v15_All_fish03';
%mod = 'All_fish03';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([fpath 'Time_Means_FOSI_' mod '_' cfile '.mat']);
load([fpath 'Space_Means_FOSI_' mod '_' cfile '.mat']);
load([fpath 'Annual_Means_FOSI_' mod '_' cfile '.mat'],'mz_mtf');
% load([fpath 'Time_Means_FOSI_' cfile '.mat']);
% load([fpath 'Space_Means_FOSI_' cfile '.mat']);
% load([fpath 'Annual_Means_FOSI_' cfile '.mat'],'mz_mtf');

% Map data
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

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

%% Plots in space

Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

Zsf(GRD.ID)=sf_sbio;
Zsp(GRD.ID)=sp_sbio;
Zsd(GRD.ID)=sd_sbio;
Zmf(GRD.ID)=mf_sbio;
Zmp(GRD.ID)=mp_sbio;
Zmd(GRD.ID)=md_sbio;
Zlp(GRD.ID)=lp_sbio;
Zld(GRD.ID)=ld_sbio;
Zb(GRD.ID)=b_sbio;

All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;
AllS = Zsp+Zsf+Zsd;
AllM = Zmp+Zmf+Zmd;
AllL = Zlp+Zld;
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracLM = AllL ./ (AllL+AllM);

%% save outputs for comparison
save([fpath 'Plot_Means_FOSI_' mod '_' cfile '.mat'],...
    'AllF','AllP','AllD','AllS','AllM','AllL','Zb','-append');

%% load means
load([fpath 'Plot_Means_FOSI_' mod '_' cfile '.mat'],...
    'AllF','AllP','AllD','AllS','AllM','AllL','Zb');
All = AllF + AllP + AllD;

%% fix seam
[TLAT2,TLON2,AllB] = cyclic_map_seam_cesm(TLAT,TLONG,Zb);
[~,~,AllS2] = cyclic_map_seam_cesm(TLAT,TLONG,AllS);
[~,~,AllM2] = cyclic_map_seam_cesm(TLAT,TLONG,AllM);
[~,~,AllL2] = cyclic_map_seam_cesm(TLAT,TLONG,AllL);
[~,~,AllF2] = cyclic_map_seam_cesm(TLAT,TLONG,AllF);
[~,~,AllP2] = cyclic_map_seam_cesm(TLAT,TLONG,AllP);
[~,~,AllD2] = cyclic_map_seam_cesm(TLAT,TLONG,AllD);
[~,~,All2]  = cyclic_map_seam_cesm(TLAT,TLONG,All);

%% 8plot by type
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - s
subplot('Position',[0.015 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllS2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2.5])
set(gcf,'renderer','painters')
text(0,1.75,'Sm','HorizontalAlignment','center')

%B - m
subplot('Position',[0.015 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllM2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2.5])
set(gcf,'renderer','painters')
text(0,1.75,'Md','HorizontalAlignment','center')

%C - l
subplot('Position',[0.015 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllL2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2.5])
set(gcf,'renderer','painters')
text(0,1.75,'Lg','HorizontalAlignment','center')

%D - F
subplot('Position',[0.015 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(All2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2.5])
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')

%E - P
subplot('Position',[0.47 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllF2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2.5])
set(gcf,'renderer','painters')
text(0,1.75,'F','HorizontalAlignment','center')

%F - D
subplot('Position',[0.47 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllP2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2.5])
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'P','HorizontalAlignment','center')

%G - B
subplot('Position',[0.47 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllD2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2.5])
set(gcf,'renderer','painters')
text(0,1.75,'D','HorizontalAlignment','center')

%H - all
subplot('Position',[0.47 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllB))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 2.5])
set(gcf,'renderer','painters')
text(0,1.75,'B','HorizontalAlignment','center')

print('-dpng',[ppath 'Map_FEISTY_FOSI_',mod,'_mean_biomass_types.png'])
