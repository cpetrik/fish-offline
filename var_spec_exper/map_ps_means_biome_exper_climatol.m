% Calc mean and std of powerspectra from COBALT ts for diff var
% By 4 biomes, N&S hemispheres separate
% Make maps of mean ps slope
% Make maps of difference from forcing ps slope

clear all
close all

%% COBALT --------------------------------------------------------
spath='/Volumes/MIP/GCM_DATA/CORE-forced/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
exper = 'Biome_exper_TP_';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/' exper];
ppath = [pp cfile '/Biome_exper/'];

% power spectra
load([fpath 'powerspec_feisty_means_stds.mat'])

%% Biomes
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([cpath 'COBALT_hist_biomes_1950_2005.mat']);
% Grid info
load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

% Sep hemis
nid = find(geolat_t(:)>0);
sid = find(geolat_t(:)<0);
vmask = lmask;
vmask(nid) = vmask(nid)*2;
nhem = find(vmask==2);
shem = find(vmask==1);

biome8_hist = biome4_hist;
biome8_hist(nid) = biome8_hist(nid) +4;
%1 & 5: LC
%2 & 6: ECCS
%3 & 7: ECSS
%4 & 8: Coastal
% pcolor(biome8_hist)
% shading flat

%% Put means on grid
[ni,nj] = size(biome8_hist);
gsf = NaN*ones(ni*nj,11);
gsp = NaN*ones(ni*nj,11);
gsd = NaN*ones(ni*nj,11);
gmf = NaN*ones(ni*nj,11);
gmp = NaN*ones(ni*nj,11);
gmd = NaN*ones(ni*nj,11);
glp = NaN*ones(ni*nj,11);
gld = NaN*ones(ni*nj,11);
gB  = NaN*ones(ni*nj,11);
gF  = NaN*ones(ni*nj,11);
gP  = NaN*ones(ni*nj,11);
gD  = NaN*ones(ni*nj,11);
gS  = NaN*ones(ni*nj,11);
gM  = NaN*ones(ni*nj,11);
gL  = NaN*ones(ni*nj,11);
gall = NaN*ones(ni*nj,11);

for i=1:8
    id = find(biome8_hist==i);
    gsf(id,:) = repmat(msf(:,i)',length(id),1);
    gsp(id,:) = repmat(msp(:,i)',length(id),1);
    gsd(id,:) = repmat(msd(:,i)',length(id),1);
    gmf(id,:) = repmat(mmf(:,i)',length(id),1);
    gmp(id,:) = repmat(mmp(:,i)',length(id),1);
    gmd(id,:) = repmat(mmd(:,i)',length(id),1);
    glp(id,:) = repmat(mlp(:,i)',length(id),1);
    gld(id,:) = repmat(mld(:,i)',length(id),1);
    gB(id,:) = repmat(mB(:,i)',length(id),1);
    gF(id,:) = repmat(mF(:,i)',length(id),1);
    gP(id,:) = repmat(mP(:,i)',length(id),1);
    gD(id,:) = repmat(mD(:,i)',length(id),1);
    gS(id,:) = repmat(mS(:,i)',length(id),1);
    gM(id,:) = repmat(mM(:,i)',length(id),1);
    gL(id,:) = repmat(mL(:,i)',length(id),1);
    gall(id,:) = repmat(mall(:,i)',length(id),1);
end

%% reshape
gsf = reshape(gsf,ni,nj,11);
gsp = reshape(gsp,ni,nj,11);
gsd = reshape(gsd,ni,nj,11);
gmf = reshape(gmf,ni,nj,11);
gmp = reshape(gmp,ni,nj,11);
gmd = reshape(gmd,ni,nj,11);
glp = reshape(glp,ni,nj,11);
gld = reshape(gld,ni,nj,11);
gB  = reshape(gB,ni,nj,11);
gF  = reshape(gF,ni,nj,11);
gP  = reshape(gP,ni,nj,11);
gD  = reshape(gD,ni,nj,11);
gS  = reshape(gS,ni,nj,11);
gM  = reshape(gM,ni,nj,11);
gL  = reshape(gL,ni,nj,11);
gall = reshape(gall,ni,nj,11);

%% Plot mean for each biome and slope
biome = {'SLC','SECCS','SECSS','SCoast','NLC','NECCS','NECSS','NCoast'};
alp = [0:-0.15:-1.5];

cm9=[0.5 0.5 0.5 ;...   %grey
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0 0.7 0;...    %green
    0 0 0];...      %black
    
set(groot,'defaultAxesColorOrder',cm9);

% plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];
% ENTER -100 TO MAP ORIGIN LONG

cmapN = cmocean('balance');
cmapN = flipud(cmapN);

%%
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,biome8_hist)
colormap(cm9(1:8,:))              
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%colorbar
colorbar('Ticks',1:8,...
         'TickLabels',biome)
caxis([1 8])                   
set(gcf,'renderer','painters')
title('Ocean Biomes CORE 1951-2007')
print('-dpng',[ppath 'biomes8_historic_1951-2007.png'])

%% B
figure(2)
%A 
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gB(:,:,1)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(1))])
%B 
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gB(:,:,3)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(3))])
%C 
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gB(:,:,5)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(5))])
%D 
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gB(:,:,7)))
colormap(cmapN)
colorbar('Position',[0.8 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(7))])
%E 
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gB(:,:,9)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(9))])
%F 
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gB(:,:,11)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(11))])
print('-dpng',[ppath exper 'map_ps_mean_B.png'])

%% F
figure(3)
%A 
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gF(:,:,1)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(1))])
%B 
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gF(:,:,3)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(3))])
%C 
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gF(:,:,5)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(5))])
%D 
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gF(:,:,7)))
colormap(cmapN)
colorbar('Position',[0.8 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(7))])
%E 
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gF(:,:,9)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(9))])
%F 
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gF(:,:,11)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(11))])
print('-dpng',[ppath exper 'map_ps_mean_F.png'])

%% P
figure(4)
%A 
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gP(:,:,1)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(1))])
%B 
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gP(:,:,3)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(3))])
%C 
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gP(:,:,5)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(5))])
%D 
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gP(:,:,7)))
colormap(cmapN)
colorbar('Position',[0.8 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(7))])
%E 
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gP(:,:,9)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(9))])
%F 
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gP(:,:,11)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Large pelagic','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(11))])
print('-dpng',[ppath exper 'map_ps_mean_P.png'])

%% D
figure(5)
%A 
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gD(:,:,1)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(1))])
%B 
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gD(:,:,3)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(3))])
%C 
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gD(:,:,5)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(5))])
%D 
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gD(:,:,7)))
colormap(cmapN)
colorbar('Position',[0.8 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(7))])
%E 
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gD(:,:,9)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(9))])
%F 
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gD(:,:,11)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(11))])
print('-dpng',[ppath exper 'map_ps_mean_D.png'])

%% S
figure(7)
%A 
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gS(:,:,1)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Small fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(1))])
%B 
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gS(:,:,3)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Small fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(3))])
%C 
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gS(:,:,5)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Small fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(5))])
%D 
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gS(:,:,7)))
colormap(cmapN)
colorbar('Position',[0.8 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Small fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(7))])
%E 
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gS(:,:,9)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Small fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(9))])
%F 
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gS(:,:,11)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Small fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(11))])
print('-dpng',[ppath exper 'map_ps_mean_S.png'])

%% M
figure(8)
%A 
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gM(:,:,1)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Medium fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(1))])
%B 
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gM(:,:,3)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Medium fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(3))])
%C 
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gM(:,:,5)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Medium fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(5))])
%D 
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gM(:,:,7)))
colormap(cmapN)
colorbar('Position',[0.8 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Medium fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(7))])
%E 
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gM(:,:,9)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Medium fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(9))])
%F 
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gM(:,:,11)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Medium fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(11))])
print('-dpng',[ppath exper 'map_ps_mean_M.png'])

%% L
figure(9)
%A 
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gL(:,:,1)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Large fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(1))])
%B 
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gL(:,:,3)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Large fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(3))])
%C 
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gL(:,:,5)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Large fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(5))])
%D 
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gL(:,:,7)))
colormap(cmapN)
colorbar('Position',[0.8 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Large fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(7))])
%E 
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gL(:,:,9)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Large fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(9))])
%F 
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gL(:,:,11)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'Large fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(11))])
print('-dpng',[ppath exper 'map_ps_mean_L.png'])

%% All
figure(10)
%A 
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gall(:,:,1)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'All fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(1))])
%B 
subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gall(:,:,3)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'All fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(3))])
%C 
subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gall(:,:,5)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'All fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(5))])
%D 
subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gall(:,:,7)))
colormap(cmapN)
colorbar('Position',[0.8 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'All fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(7))])
%E 
subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gall(:,:,9)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'All fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(9))])
%F 
subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,squeeze(gall(:,:,11)))
colormap(cmapN)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'All fish','HorizontalAlignment','center')
text(-2.75,1.75,['\alpha = ' num2str(alp(11))])
print('-dpng',[ppath exper 'map_ps_mean_All.png'])

