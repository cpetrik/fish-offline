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

lme_medres = NaN*ones(66,8);
lme_sres = NaN*ones(66,8);
lme_area = NaN*ones(66,1);

lme_type = NaN*ones(66,4);

for L=1:66
    lid = find(tlme==L);
    %plain means
    lme_medres(L,1) = median(Csf(lid),'omitnan');
    lme_medres(L,2) = median(Csp(lid),'omitnan');
    lme_medres(L,3) = median(Csd(lid),'omitnan');
    lme_medres(L,4) = median(Cmf(lid),'omitnan');
    lme_medres(L,5) = median(Cmp(lid),'omitnan');
    lme_medres(L,6) = median(Cmd(lid),'omitnan');
    lme_medres(L,7) = median(Clp(lid),'omitnan');
    lme_medres(L,8) = median(Cld(lid),'omitnan');

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

%%
save([fpath 'LME_fosi_fished_',harv,cfile '.mat'],...
    'lme_medres','lme_awmres','-append');

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

    Msf(id) = lme_medres(i,1);
    Msp(id) = lme_medres(i,2);
    Msd(id) = lme_medres(i,3);
    Mmf(id) = lme_medres(i,4);
    Mmp(id) = lme_medres(i,5);
    Mmd(id) = lme_medres(i,6);
    Mlp(id) = lme_medres(i,7);
    Mld(id) = lme_medres(i,8);
   
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
clim([0 45])
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
clim([0 45])
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
clim([0 45])
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
clim([0 1])
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
clim([0 2])
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
clim([0 4])
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
clim([0 4])
c8=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c8,'y')
set(gcf,'renderer','painters')
text(0,1.75,'LD res','HorizontalAlignment','center')
print('-dpng',[ppath 'FOSI_map_LMEmedian_resMean_stages.png'])


%% map - 8 plot of area-weighted mean res
f6 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,(Rsf))
cmocean('speed')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 45])
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
clim([0 45])
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
clim([0 45])
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
clim([0 1])
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
clim([0 2])
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
clim([0 4])
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
clim([0 4])
c8=colorbar('orientation','vertical','AxisLocation','out');
xlabel(c8,'y')
set(gcf,'renderer','painters')
text(0,1.75,'LD res','HorizontalAlignment','center')
print('-dpng',[ppath 'FOSI_map_LMEmean_areaweighted_resMean_stages.png'])


%% means of each type
lme_type(:,1) = mean(lme_medres(:,1:3),2,'omitnan');
lme_type(:,2) = mean(lme_medres(:,4:6),2,'omitnan');
lme_type(:,3) = mean(lme_medres(:,7:8),2,'omitnan');
lme_type(:,4) = mean(lme_medres(:,[1,4]),2,'omitnan');
lme_type(:,5) = mean(lme_medres(:,[2,5,7]),2,'omitnan');
lme_type(:,6) = mean(lme_medres(:,[3,6,8]),2,'omitnan');

%% Save as table for HGL
lnum = 1:66;
ltex = char(lnum);

Ltab = array2table(lme_medres,'VariableNames',...
    {'SF','SP','SD','MF','MP','MD','LP','LD'});

Ttab = array2table(lme_type,'VariableNames',...
    {'S','M','L','F','P','D'});

writetable(Ltab,[fpath 'Residence_time_smeans_FOSI_' harv 'stages.csv'],'WriteRowNames',true);
writetable(Ttab,[fpath 'Residence_time_smeans_FOSI_' harv 'types.csv'],'WriteRowNames',true);

