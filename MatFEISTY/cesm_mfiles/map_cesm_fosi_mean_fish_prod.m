% Visualize output of FEISTY prod
% CESM FOSI
% maps

clear 
close all

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v15_All_fish03';
%mod = 'All_fish03';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
ppath = [pp cfile '/FOSI/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

load([fpath 'Space_Means_FOSI_' mod '_' cfile '.mat'],'time',...
    'sf_sprod','sp_sprod','sd_sprod',...
    'mf_sprod','mp_sprod','md_sprod',...
    'lp_sprod','ld_sprod');

% Map data
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

[ni,nj]=size(TLONG);
ID=GRD.ID;

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

Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);

Psf(ID)=sf_sprod;
Psp(ID)=sp_sprod;
Psd(ID)=sd_sprod;
Pmf(ID)=mf_sprod;
Pmp(ID)=mp_sprod;
Pmd(ID)=md_sprod;
Plp(ID)=lp_sprod;
Pld(ID)=ld_sprod;

Psf(Psf(:)<0) = nan; 
Psp(Psp(:)<0) = nan;
Psd(Psd(:)<0) = nan;
Pmf(Pmf(:)<0) = nan;
Pmp(Pmp(:)<0) = nan;
Pmd(Pmd(:)<0) = nan;
Plp(Plp(:)<0) = nan;
Pld(Pld(:)<0) = nan;

All = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
% AllF = Psf+Pmf;
% AllP = Psp+Pmp+Plp;
% AllD = Psd+Pmd+Pld;
AllS = Psp+Psf+Psd;
AllM = Pmp+Pmf+Pmd;
AllL = Plp+Pld;
% does this need to be biomass-weighted mean? 
%- no, biomass already part of it

%% save outputs for comparison
save([fpath 'Plot_Prod_Means_FOSI_' mod '_' cfile '.mat'],...
    'Pmf','Plp','Pld','AllS','AllM','AllL','All');

%% load means
load([fpath 'Plot_Prod_Means_FOSI_' mod '_' cfile '.mat'],...
    'Pmf','Plp','Pld','AllS','AllM','AllL','All');

%% fix seam
[TLAT2,TLON2,AllS2] = cyclic_map_seam_cesm(TLAT,TLONG,AllS);
[~,~,AllM2] = cyclic_map_seam_cesm(TLAT,TLONG,AllM);
[~,~,AllL2] = cyclic_map_seam_cesm(TLAT,TLONG,AllL);
[~,~,AllF2] = cyclic_map_seam_cesm(TLAT,TLONG,Pmf);
[~,~,AllP2] = cyclic_map_seam_cesm(TLAT,TLONG,Plp);
[~,~,AllD2] = cyclic_map_seam_cesm(TLAT,TLONG,Pld);
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
caxis([-4 0])
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
caxis([-4 0])
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
caxis([-4 0])
set(gcf,'renderer','painters')
text(0,1.75,'Lg','HorizontalAlignment','center')

%D - All
subplot('Position',[0.2425 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(All2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0])
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')

%E - F
subplot('Position',[0.47 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllF2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0])
set(gcf,'renderer','painters')
text(0,1.75,'F','HorizontalAlignment','center')

%F - P
subplot('Position',[0.47 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllP2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0])
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'P','HorizontalAlignment','center')

%G - D
subplot('Position',[0.47 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllD2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0])
set(gcf,'renderer','painters')
text(0,1.75,'D','HorizontalAlignment','center')

%H - all
% subplot('Position',[0.47 0.0 0.44 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(TLAT2,TLON2,log10(AllB))
% colormap(cmBP50)
% %colorbar
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-4 0])
% set(gcf,'renderer','painters')
% text(0,1.75,'B','HorizontalAlignment','center')

print('-dpng',[ppath 'Map_FEISTY_FOSI_',mod,'_mean_prod_types.png'])

%% 8plot by stage
f3 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - s
subplot('Position',[0.015 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,log10(Psf))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0])
set(gcf,'renderer','painters')
text(0,1.75,'SF','HorizontalAlignment','center')

%B - m
subplot('Position',[0.015 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,log10(Psp))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0])
set(gcf,'renderer','painters')
text(0,1.75,'SP','HorizontalAlignment','center')

%C - l
subplot('Position',[0.015 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,log10(Pmp))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0])
set(gcf,'renderer','painters')
text(0,1.75,'MP','HorizontalAlignment','center')

%D - F
subplot('Position',[0.015 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,log10(Plp))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0])
set(gcf,'renderer','painters')
text(0,1.75,'LP','HorizontalAlignment','center')

%E - P
subplot('Position',[0.47 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,log10(Pmf))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0])
set(gcf,'renderer','painters')
text(0,1.75,'MF','HorizontalAlignment','center')

%F - D
subplot('Position',[0.47 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,log10(Psd))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0])
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'SD','HorizontalAlignment','center')

%G - B
subplot('Position',[0.47 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,log10(Pmd))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0])
set(gcf,'renderer','painters')
text(0,1.75,'MD','HorizontalAlignment','center')

%H - all
subplot('Position',[0.47 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,log10(Pld))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 0])
set(gcf,'renderer','painters')
text(0,1.75,'LD','HorizontalAlignment','center')

print('-dpng',[ppath 'Map_FEISTY_FOSI_',mod,'_mean_prod_stage.png'])
