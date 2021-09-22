% Visualize difference between
% ESM2M Hindcast 50yr mean and
% CESM FOSI mean

clear all
close all

%% COBALT Hindcast grid
hpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load([hpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([hpath 'grid_csv.csv']); %grid
[hi,hj]=size(geolon_t);

%% CESM FOSI grid
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6.mat']);
load([cpath 'Data_grid_POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

%% FEISTY Output
cfile2 = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
mod = 'v13_sMZ090_mMZ045_All_fish03_';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
ppath = [pp cfile2 '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%COBALT
cfile1 = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
gpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile1 '/Historic_ESM2M/'];
load([gpath 'Means_Historic_',harv,'_' cfile1 '.mat'],...
    'sf_mean50','sp_mean50','sd_mean50',...
    'mf_mean50','mp_mean50','md_mean50',...
    'lp_mean50','ld_mean50','b_mean50');

%CESM
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile2 '/'];
load([fpath 'Space_Means_FOSI_' mod cfile2 '.mat']);

%% Put biomass on grid
Csf=NaN*ones(ni,nj);
Csp=NaN*ones(ni,nj);
Csd=NaN*ones(ni,nj);
Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);
Cb =NaN*ones(ni,nj);
Csf(ID)=sf_sbio;
Csp(ID)=sp_sbio;
Csd(ID)=sd_sbio;
Cmf(ID)=mf_sbio;
Cmp(ID)=mp_sbio;
Cmd(ID)=md_sbio;
Clp(ID)=lp_sbio;
Cld(ID)=ld_sbio;
Cb(ID) =b_sbio;

Hsf=NaN*ones(hi,hj);
Hsp=NaN*ones(hi,hj);
Hsd=NaN*ones(hi,hj);
Hmf=NaN*ones(hi,hj);
Hmp=NaN*ones(hi,hj);
Hmd=NaN*ones(hi,hj);
Hlp=NaN*ones(hi,hj);
Hld=NaN*ones(hi,hj);
Hb =NaN*ones(hi,hj);
Hsf(grid(:,1))=sf_mean50;
Hsp(grid(:,1))=sp_mean50;
Hsd(grid(:,1))=sd_mean50;
Hmf(grid(:,1))=mf_mean50;
Hmp(grid(:,1))=mp_mean50;
Hmd(grid(:,1))=md_mean50;
Hlp(grid(:,1))=lp_mean50;
Hld(grid(:,1))=ld_mean50;
Hb(grid(:,1)) =b_mean50;

CF = Csf+Cmf;
CP = Csp+Cmp+Clp;
CD = Csd+Cmd+Cld;
CS = Csp+Csf+Csd;
CM = Cmp+Cmf+Cmd;
CL = Clp+Cld;

HF = Hsf+Hmf;
HP = Hsp+Hmp+Hlp;
HD = Hsd+Hmd+Hld;
HS = Hsp+Hsf+Hsd;
HM = Hmp+Hmf+Hmd;
HL = Hlp+Hld;

%% Interpolate to same grid
%tlat       [-79.2205 89.7064]
%geolat_t   [-81.5 89.4879]
%tlon       [0.0147 359.996]
%geolon_t   [-279.9803 79.9803]

%%
figure
pcolor(TLONG)
shading flat
colorbar

figure
pcolor(geolon_t)
shading flat
colorbar

%% Need to fix both longitudes
test = TLONG;
id=find(test>180);
test(id)=test(id)-360;
lon = test;

clat = TLAT';
clon = lon';

test2=geolon_t;
id=find(test2<-180);
test2(id)=test2(id)+360;
geolon = test2;

geolat = geolat_t';
geolon = geolon';

lats = -89.5:89.5;
lons = -179.5:179.5;
[glon,glat] = meshgrid(lons,lats);

%%
hF = griddata(geolat,geolon,HF',glat,glon);
hP = griddata(geolat,geolon,HP',glat,glon);
hD = griddata(geolat,geolon,HD',glat,glon);
hB = griddata(geolat,geolon,Hb',glat,glon);
hS = griddata(geolat,geolon,HS',glat,glon);
hM = griddata(geolat,geolon,HM',glat,glon);
hL = griddata(geolat,geolon,HL',glat,glon);

cF = griddata(clat,clon,CF',glat,glon);
cP = griddata(clat,clon,CP',glat,glon);
cD = griddata(clat,clon,CD',glat,glon);
cB = griddata(clat,clon,Cb',glat,glon);
cS = griddata(clat,clon,CS',glat,glon);
cM = griddata(clat,clon,CM',glat,glon);
cL = griddata(clat,clon,CL',glat,glon);

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

load coastlines

%%
cAll = cF+cP+cD;
cFracPD = cP ./ (cP+cD);
cFracPF = cP ./ (cP+cF);
cFracLM = cL ./ (cL+cM);

hAll = hF+hP+hD;
hFrahPD = hP ./ (hP+hD);
hFrahPF = hP ./ (hP+hF);
hFrahLM = hL ./ (hL+hM);

%
diffF = (cF./hF);
diffP = (cP./hP);
diffD = (cD./hD);
diffB = (cB./hB);
diffAll = (cAll./hAll);

pdiffF = (cF-hF) ./ hF;
pdiffP = (cP-hP) ./ hP;
pdiffD = (cD-hD) ./ hD;
pdiffB = (cB-hB) ./ hB;
pdiffAll = (cAll-hAll) ./ hAll;

%% Maps
figure(1)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hF)))
colormap('jet')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('COBALT F');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(cF)))
colormap('jet')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('CESM F');
print('-dpng',[ppath 'CESM_COBALT_' mod 'global_F.png'])

%P
figure(2)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hP)))
colormap('jet')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('COBALT P');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(cP)))
colormap('jet')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('CESM P');
print('-dpng',[ppath 'CESM_COBALT_' mod 'global_P.png'])

% D
figure(3)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hD)))
colormap('jet')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('COBALT D');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(cD)))
colormap('jet')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('CESM D');
print('-dpng',[ppath 'CESM_COBALT_' mod 'global_D.png'])

%4
figure(4)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hAll)))
colormap('jet')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
colorbar
set(gcf,'renderer','painters')
title('COBALT All');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(cAll)))
colormap('jet')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
colorbar
set(gcf,'renderer','painters')
title('CESM All');
print('-dpng',[ppath 'CESM_COBALT_' mod 'global_All.png'])

%% B
figure(5)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(hB)))
colormap('jet')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('COBALT Benthos');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) %,'origin',[0 -100 0])
surfm(glat,glon,real(log10(cB)))
colormap('jet')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 0.5]);
colorbar
set(gcf,'renderer','painters')
title('CESM Benthos');
print('-dpng',[ppath 'CESM_COBALT_' mod 'global_Bent.png'])

%% side by side on one fig
cmBP50=cbrewer('seq','BuPu',50,'PCHIP');

f1 = figure('Units','inches','Position',[1 3 7.5 10]);
%f1.Units = 'inches';

%1 - Hist cmcc
subplot('Position',[0.025 0.825 0.4 0.165])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,real(log10(hAll)))
colormap(cmBP50)
caxis([-1 2])
text(0,1.75,'ESM2M-CORE','HorizontalAlignment','center','FontWeight','bold')
text(-1.75,1.75,'All','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%2 - Hist cnrm
subplot('Position',[0.025 0.66 0.4 0.165])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,real(log10(hF)))
colormap(cmBP50)
caxis([-1 2])
text(-1.75,1.75,'Forage','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%3 - Hist gfdl
subplot('Position',[0.025 0.495 0.4 0.165])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,real(log10(hP)))
colormap(cmBP50)
caxis([-1 2])
text(-1.75,1.75,'Lg Pelagic','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%4 - CESM D
subplot('Position',[0.025 0.33 0.4 0.165])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,real(log10(hD)))
colormap(cmBP50)
caxis([-1 2])
text(-1.75,1.75,'Demersal','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%5 - CESM B
subplot('Position',[0.025 0.165 0.4 0.165])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,real(log10(hB)))
colormap(cmBP50)
caxis([-1 2])
text(-1.75,1.75,'Benthos','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%6 - Hist obs
% subplot('Position',[0.025 0.0 0.4 0.165])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(lat_o,lon_o,biomes_o)
% colormap(cmBP50)
% caxis([-1 2])
% text(-1.75,1.75,'obs','HorizontalAlignment','center')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

% 7 - GFDL All
subplot('Position',[0.43 0.825 0.4 0.165])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,real(log10(cAll)))
colormap(cmBP50)
caxis([-1 2])
text(0,1.75,'CESM-FOSI','HorizontalAlignment','center','FontWeight','bold')
text(-1.75,1.75,'All','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%8 - GFDL F
subplot('Position',[0.43 0.66 0.4 0.165])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,real(log10(cF)))
colormap(cmBP50)
caxis([-1 2])
text(-1.75,1.75,'Forage','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%9 - GFDL P
subplot('Position',[0.43 0.495 0.4 0.165])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,real(log10(cP)))
colormap(cmBP50)
caxis([-1 2])
text(-1.75,1.75,'Lg Pelagic','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
cb = colorbar('Position',[0.8 0.45 0.025 0.25]);
xlabel(cb,'biomass (log_1_0 g m^-^2)')

%10 - GFDL D
subplot('Position',[0.43 0.33 0.4 0.165])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,real(log10(cD)))
colormap(cmBP50)
caxis([-1 2])
text(-1.75,1.75,'Demersal','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%11 - GFDL B
subplot('Position',[0.43 0.165 0.4 0.165])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(glat,glon,real(log10(cB)))
colormap(cmBP50)
caxis([-1 2])
text(-1.75,1.75,'Benthos','HorizontalAlignment','center')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[ppath 'CESM_COBALT_' mod 'global_all_types.png']);

%% diffs
figure(7)
% All F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffF)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('CESM / COBALT F');

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffP)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
set(gcf,'renderer','painters')
title('CESM / COBALT P');

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffD)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
set(gcf,'renderer','painters')
title('CESM / COBALT D');

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffAll)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
set(gcf,'renderer','painters')
title('CESM / COBALT All');
stamp('')
print('-dpng',[ppath 'CESM_COBALT_' mod 'global_diff_types.png'])

%% pdiffs relative to COBALT
figure(8)
% All F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffF)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('CESM - COBALT F');

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffP)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('CESM - COBALT P');

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffD)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('CESM - COBALT D');

% All
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffAll)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('CESM - COBALT All');
stamp('')
print('-dpng',[ppath 'CESM_COBALT_' mod 'global_pdiff_types.png'])

%% B
figure(10)
subplot(2,1,1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,diffB)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-10 10]);
colorbar
set(gcf,'renderer','painters')
title('CESM / COBALT B');

subplot(2,1,2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1) 
surfm(glat,glon,pdiffB)
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar
set(gcf,'renderer','painters')
title('CESM - COBALT B');
stamp('')
print('-dpng',[ppath 'CESM_COBALT_' mod 'global_diffs_B.png'])

%% Calc differences in total biomass


