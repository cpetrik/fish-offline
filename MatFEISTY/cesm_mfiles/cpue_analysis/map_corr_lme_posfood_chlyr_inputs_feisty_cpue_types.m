% Mapp corr coefs of ind driver combos with cpue
% Subplots together for comparison
% Restricted analysis to chl yrs

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regress_cpue/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',...
    cfile,'/corrs_cpue/'];

mod = 'v15_All_fish03';
mod2 = 'v15_obsfish';

%%
load('LMEs_corr_cpue_chlyrs_inputs_feisty_mostsiglag_posfood.mat')
%dim: LME x driver x (corr, p-val, lag)

%% Map data
cpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
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

cmRB = cbrewer('div','RdBu',21,'PCHIP');
cmRB = flipud(cmRB);

%% 8plot by driver
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - F biom
subplot('Position',[0.015 0.75 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllF2))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1.5 2.0])
set(gcf,'renderer','painters')
text(0,1.75,'F biom','HorizontalAlignment','center')

%B - P biom
subplot('Position',[0.015 0.5 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllP2))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1.5 2.0])
set(gcf,'renderer','painters')
text(0,1.75,'P biom','HorizontalAlignment','center')

%C - D biom
subplot('Position',[0.015 0.25 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllD2))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1.5 2.0])
set(gcf,'renderer','painters')
text(0,1.75,'D biom','HorizontalAlignment','center')
colorbar('Position',[0.445 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out',...
    'Ticks',-1:2,'TickLabels',-1:2)

%D - A biom
subplot('Position',[0.015 0.0 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(All2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1.5 2.0])
set(gcf,'renderer','painters')
text(0,1.75,'All biom','HorizontalAlignment','center')


%E - F prod
subplot('Position',[0.51 0.75 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(pF2+eps))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1])
%colorbar('Position',[0.95 0.765 0.025 0.225],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'F prod','HorizontalAlignment','center')

%F - P prod
subplot('Position',[0.51 0.5 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(pP2+eps))
colormap(cmBP50)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1])
colorbar('Position',[0.945 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out',...
    'Ticks',-1:1,'TickLabels',-1:1)
set(gcf,'renderer','painters')
text(0,1.75,'P prod','HorizontalAlignment','center')

%G - D prod
subplot('Position',[0.51 0.25 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(pD2+eps))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1])
set(gcf,'renderer','painters')
text(0,1.75,'D','HorizontalAlignment','center')

%H - All prod
subplot('Position',[0.51 0.0 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(pA2+eps))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1])
set(gcf,'renderer','painters')
text(0,1.75,'All prod','HorizontalAlignment','center')

%print('-dpng',[ppath 'Map_FEISTY_FOSI_',mod,'_mean_biomass_prod_types.png'])


