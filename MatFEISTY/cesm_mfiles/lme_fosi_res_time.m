% Calc LME biomass of FEISTY
% CESM FOSI

clear 
close all

%% Map data
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

geolon_t = double(TLONG);
geolat_t = double(TLAT);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%% Fish
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
harv = 'v15_All_fish03_';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%low biomass regions removed
load([fpath 'Residence_time_smeans_FOSI_' harv cfile '.mat'],...
  'sf_mres','sp_mres','sd_mres','mf_mres','mp_mres','md_mres','lp_mres','ld_mres')

%% Put biomass on grid 
Csf=NaN*ones(ni,nj);
Csp=NaN*ones(ni,nj);
Csd=NaN*ones(ni,nj);
Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

Csf(ID)=sf_mres;
Csp(ID)=sp_mres;
Csd(ID)=sd_mres;
Cmf(ID)=mf_mres;
Cmp(ID)=mp_mres;
Cmd(ID)=md_mres;
Clp(ID)=lp_mres;
Cld(ID)=ld_mres;

% g/m2 --> total g
AREA_OCN = TAREA * 1e-4;
Asf_mean = Csf .* AREA_OCN;
Asp_mean = Csp .* AREA_OCN;
Asd_mean = Csd .* AREA_OCN;
Amf_mean = Cmf .* AREA_OCN;
Amp_mean = Cmp .* AREA_OCN;
Amd_mean = Cmd .* AREA_OCN;
Alp_mean = Clp .* AREA_OCN;
Ald_mean = Cld .* AREA_OCN;

%% Calc LMEs
tlme = double(lme_mask);
tlme(tlme<0) = nan;

lme_ares = NaN*ones(66,8);
lme_sres = NaN*ones(66,8);
lme_area = NaN*ones(66,1);

lme_type = NaN*ones(66,4);

for L=1:66
    lid = find(tlme==L);
    %plain means
    lme_ares(L,1) = mean(Csf(lid),'omitnan');
    lme_ares(L,2) = mean(Csp(lid),'omitnan');
    lme_ares(L,3) = mean(Csd(lid),'omitnan');
    lme_ares(L,4) = mean(Cmf(lid),'omitnan');
    lme_ares(L,5) = mean(Cmp(lid),'omitnan');
    lme_ares(L,6) = mean(Cmd(lid),'omitnan');
    lme_ares(L,7) = mean(Clp(lid),'omitnan');
    lme_ares(L,8) = mean(Cld(lid),'omitnan');

    %area-weighted means
    lme_sres(L,1) = sum(Asf_mean(lid),'omitnan');
    lme_sres(L,2) = sum(Asp_mean(lid),'omitnan');
    lme_sres(L,3) = sum(Asd_mean(lid),'omitnan');
    lme_sres(L,4) = sum(Amf_mean(lid),'omitnan');
    lme_sres(L,5) = sum(Amp_mean(lid),'omitnan');
    lme_sres(L,6) = sum(Amd_mean(lid),'omitnan');
    lme_sres(L,7) = sum(Alp_mean(lid),'omitnan');
    lme_sres(L,8) = sum(Ald_mean(lid),'omitnan');

    %LME area
    lme_area(L,1) = sum(AREA_OCN(lid),'omitnan');
end

%% Change to area-weighted means
lme_area_mat = repmat(lme_area,1,8);
lme_awmres = lme_sres ./ lme_area_mat;

%% Change anything over 50 years to 50
%lme_awmres(lme_awmres(:) > (50*365)) = (50*365);

%%
save([fpath 'LME_fosi_fished_',harv,cfile '.mat'],...
    'lme_ares','lme_awmres','-append');

%% Put LME means on grid

Msf=NaN*ones(ni,nj);
Msp=NaN*ones(ni,nj);
Msd=NaN*ones(ni,nj);
Mmf=NaN*ones(ni,nj);
Mmp=NaN*ones(ni,nj);
Mmd=NaN*ones(ni,nj);
Mlp=NaN*ones(ni,nj);
Mld=NaN*ones(ni,nj);

Rsf=NaN*ones(ni,nj);
Rsp=NaN*ones(ni,nj);
Rsd=NaN*ones(ni,nj);
Rmf=NaN*ones(ni,nj);
Rmp=NaN*ones(ni,nj);
Rmd=NaN*ones(ni,nj);
Rlp=NaN*ones(ni,nj);
Rld=NaN*ones(ni,nj);

for i=1:66
    id = find(tlme==i);

    Msf(id) = lme_ares(i,1);
    Msp(id) = lme_ares(i,2);
    Msd(id) = lme_ares(i,3);
    Mmf(id) = lme_ares(i,4);
    Mmp(id) = lme_ares(i,5);
    Mmd(id) = lme_ares(i,6);
    Mlp(id) = lme_ares(i,7);
    Mld(id) = lme_ares(i,8);
   
    Rsf(id) = lme_awmres(i,1);
    Rsp(id) = lme_awmres(i,2);
    Rsd(id) = lme_awmres(i,3);
    Rmf(id) = lme_awmres(i,4);
    Rmp(id) = lme_awmres(i,5);
    Rmd(id) = lme_awmres(i,6);
    Rlp(id) = lme_awmres(i,7);
    Rld(id) = lme_awmres(i,8);
    
end

%% map - 8 plot of area-weighted mean res
f1 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Msf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 5])
set(gcf,'renderer','painters')
text(0,1.75,'SF res','HorizontalAlignment','center')
c1=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c1,'d')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Msp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 5])
set(gcf,'renderer','painters')
text(0,1.75,'SP res','HorizontalAlignment','center')
c2=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c2,'d')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Msd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 5])
set(gcf,'renderer','painters')
text(0,1.75,'SD res','HorizontalAlignment','center')
c3=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c3,'d')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Mmf)./365)
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 2])
set(gcf,'renderer','painters')
text(0,1.75,'MF res','HorizontalAlignment','center')
c4=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c4,'y')

%E
subplot('Position',[0.495 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Mmp)./365)
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'MP res','HorizontalAlignment','center')
c5=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c5,'y')

%F
subplot('Position',[0.495 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Mmd)./365)
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1])
%cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
%xlabel(cb,'y')
set(gcf,'renderer','painters')
text(0,1.75,'MD res','HorizontalAlignment','center')
c6=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c6,'y')

%G
subplot('Position',[0.495 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Mlp)./365)
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 5])
set(gcf,'renderer','painters')
text(0,1.75,'LP res','HorizontalAlignment','center')
c7=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c7,'y')

%H
subplot('Position',[0.495 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Mld)./365)
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 5])
c8=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c8,'y')
set(gcf,'renderer','painters')
text(0,1.75,'LD res','HorizontalAlignment','center')
print('-dpng',[ppath 'FOSI_map_LMEmean_resMean_stages.png'])


%% map - 8 plot of area-weighted mean res
f6 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 5])
set(gcf,'renderer','painters')
text(0,1.75,'SF res','HorizontalAlignment','center')
c1=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c1,'d')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsp))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 5])
set(gcf,'renderer','painters')
text(0,1.75,'SP res','HorizontalAlignment','center')
c2=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c2,'d')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsd))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 5])
set(gcf,'renderer','painters')
text(0,1.75,'SD res','HorizontalAlignment','center')
c3=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c3,'d')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rmf)./365)
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 2])
set(gcf,'renderer','painters')
text(0,1.75,'MF res','HorizontalAlignment','center')
c4=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c4,'y')

%E
subplot('Position',[0.495 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rmp)./365)
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'MP res','HorizontalAlignment','center')
c5=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c5,'y')

%F
subplot('Position',[0.495 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rmd)./365)
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1])
%cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
%xlabel(cb,'y')
set(gcf,'renderer','painters')
text(0,1.75,'MD res','HorizontalAlignment','center')
c6=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c6,'y')

%G
subplot('Position',[0.495 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rlp)./365)
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 5])
set(gcf,'renderer','painters')
text(0,1.75,'LP res','HorizontalAlignment','center')
c7=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c7,'y')

%H
subplot('Position',[0.495 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rld)./365)
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 5])
c8=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c8,'y')
set(gcf,'renderer','painters')
text(0,1.75,'LD res','HorizontalAlignment','center')
print('-dpng',[ppath 'FOSI_map_LMEmean_areaweighted_resMean_stages.png'])



