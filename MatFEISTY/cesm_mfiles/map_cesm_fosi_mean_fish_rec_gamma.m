% Visualize output of FEISTY recruitment & gamma
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
    'mf_srec','lp_srec','ld_srec',...
    'mf_sgam','lp_sgam','ld_sgam');

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

Rmf=NaN*ones(ni,nj);
Rlp=NaN*ones(ni,nj);
Rld=NaN*ones(ni,nj);

Gmf=NaN*ones(ni,nj);
Glp=NaN*ones(ni,nj);
Gld=NaN*ones(ni,nj);

Rmf(ID)=mf_srec;
Rlp(ID)=lp_srec;
Rld(ID)=ld_srec;

Gmf(ID)=mf_sgam;
Glp(ID)=lp_sgam;
Gld(ID)=ld_sgam;

RA = (Rmf+Rlp+Rld)./3;
GA = (Gmf+Glp+Gld)./3;

%% save outputs for comparison
save([fpath 'Plot_rec_gamma_Means_FOSI_' mod '_' cfile '.mat'],...
    'Rmf','Rlp','Rld','RA','Gmf','Glp','Gld','GA');

%% load means
load([fpath 'Plot_rec_gamma_Means_FOSI_' mod '_' cfile '.mat'],...
    'Rmf','Rlp','Rld','RA','Gmf','Glp','Gld','GA');

%%
Rmf(Rmf(:)<0) = nan;
Rlp(Rlp(:)<0) = nan;
Rld(Rld(:)<0) = nan;
RA(RA(:)<0) = nan;

Gmf(Gmf(:)<0) = nan;
Glp(Glp(:)<0) = nan;
Gld(Gld(:)<0) = nan;
GA(GA(:)<0) = nan;

%% fix seam
[TLAT2,TLON2,AllF2] = cyclic_map_seam_cesm(TLAT,TLONG,Rmf);
[~,~,AllP2] = cyclic_map_seam_cesm(TLAT,TLONG,Rlp);
[~,~,AllD2] = cyclic_map_seam_cesm(TLAT,TLONG,Rld);
[~,~,All] = cyclic_map_seam_cesm(TLAT,TLONG,RA);

[~,~,GF2] = cyclic_map_seam_cesm(TLAT,TLONG,Gmf);
[~,~,GP2] = cyclic_map_seam_cesm(TLAT,TLONG,Glp);
[~,~,GD2] = cyclic_map_seam_cesm(TLAT,TLONG,Gld);
[~,~,Gall] = cyclic_map_seam_cesm(TLAT,TLONG,GA);

%% 8plot by type
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%E - F
subplot('Position',[0.47 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllF2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 -2])
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
caxis([-4 -2])
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
caxis([-4 -2])
set(gcf,'renderer','painters')
text(0,1.75,'D','HorizontalAlignment','center')

%All
subplot('Position',[0.47 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(All))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 -2])
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')
stamp('rec')

print('-dpng',[ppath 'Map_FEISTY_FOSI_',mod,'_mean_rec_types.png'])

%% Gamma
f3 = figure('Units','inches','Position',[1 3 6.5 8]);

%E - F
subplot('Position',[0.47 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(GF2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 -1])
set(gcf,'renderer','painters')
text(0,1.75,'F','HorizontalAlignment','center')

%F - P
subplot('Position',[0.47 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(GP2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 -1])
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'P','HorizontalAlignment','center')

%G - D
subplot('Position',[0.47 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(GD2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 -1])
set(gcf,'renderer','painters')
text(0,1.75,'D','HorizontalAlignment','center')

%All
subplot('Position',[0.47 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(Gall))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 -1])
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')
stamp('gamma')

print('-dpng',[ppath 'Map_FEISTY_FOSI_',mod,'_mean_gamma_types.png'])
