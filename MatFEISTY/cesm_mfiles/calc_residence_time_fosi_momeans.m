% Calc residence time
% Residence time = biomass / input
% or             = biomass / output
% Total inputs: rec
% Total outputs: gamma, rep, nmort, die (pred), yield (fishing)
% Exclude cells where biomass < 1 fish per grid cell
% Use monthly means instead of 68-yr mean

clear
close all

%%
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;
AREA_OCN = TAREA * 1e-4;

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
harv = 'v15_All_fish03_';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

load([fpath 'Annual_Means_FOSI_' harv cfile '.mat']);

%% add gains and losses
% rec + ((nu - rep - gamma - nmort - fmort) .* bio_in) - die
% rec, die, prod, yield are g/m2/d
% rep, pred, nu, gamma, nmort, fmort are g/g/m2/d

% sf_in = max(0,sf_anu) + sf_arec;
% sp_in = max(0,sp_anu) + sp_arec;
% sd_in = max(0,sd_anu) + sd_arec;
% mp_in = max(0,mp_anu) + mp_arec;
% md_in = max(0,md_anu) + md_arec;
% mf_in = max(0,mf_anu) + mf_arec;
% lp_in = max(0,lp_anu) + lp_arec;
% ld_in = max(0,ld_anu) + ld_arec;

%All g/m2/d
sf_in = max(0,sf_aprod) + sf_arec;
sp_in = max(0,sp_aprod) + sp_arec;
sd_in = max(0,sd_aprod) + sd_arec;
mp_in = max(0,mp_aprod) + mp_arec;
md_in = max(0,md_aprod) + md_arec;
mf_in = max(0,mf_aprod) + mf_arec;
lp_in = max(0,lp_aprod) + lp_arec;
ld_in = max(0,ld_aprod) + ld_arec;

%Have gamma defined as growth into size class. e.g. mp_agam = sml_p_gamma
% sf_out = mf_agam + sf_amort + sf_adie;
% sp_out = mp_agam + sp_amort + sp_adie;
% sd_out = md_agam + sd_amort + sd_adie;
% mp_out = lp_agam + mp_amort + mp_adie + mp_amy;
% md_out = ld_agam + md_amort + md_adie + md_amy;
% mf_out = mf_arep + mf_amort + mf_adie + mf_amy;
% lp_out = lp_arep + lp_amort + lp_adie + lp_amy;
% ld_out = ld_arep + ld_amort + ld_adie + ld_amy;

%All g/m2/d
sf_out = (mf_agam.*sf_abio) + (sf_amort.*sf_abio) + sf_adie;
sp_out = mp_agam + sp_amort + sp_adie;
sd_out = md_agam + sd_amort + sd_adie;
mp_out = lp_agam + mp_amort + mp_adie + mp_may;
md_out = ld_agam + md_amort + md_adie + md_may;
mf_out = mf_arep + mf_amort + mf_adie + mf_may;
lp_out = lp_arep + lp_amort + lp_adie + lp_may;
ld_out = ld_arep + ld_amort + ld_adie + ld_may;

%% I THINK I NEED TO REMOVE LOW BIOMASS AREAS

nid=length(ID);
varea = AREA_OCN(ID);

% gr_bio_s = repmat([0.001, 0.02, 0.5],nid,1) ./ repmat(varea,1,3);
% gr_bio_m = repmat([0.5, 11.2, 250],nid,1) ./ repmat(varea,1,3);
% gr_bio_l = repmat([250, 5600, 125000],nid,1) ./ repmat(varea,1,3);

gr_bio_s = repmat(0.001,nid,68) ./ repmat(varea,1,68);
gr_bio_m = repmat(0.5,nid,68) ./ repmat(varea,1,68);
gr_bio_l = repmat(250,nid,68) ./ repmat(varea,1,68);

%% test lowest first
sid = (sf_sbio(:)<gr_bio_s(:));
sf_sbio(sid) = nan;

sid = (sp_sbio(:)<gr_bio_s(:));
sp_sbio(sid) = nan;

sid = (sd_sbio(:)<gr_bio_s(:));
sd_sbio(sid) = nan;

sid = (mf_sbio(:)<gr_bio_m(:));
mf_sbio(sid) = nan;

sid = (mp_sbio(:)<gr_bio_m(:));
mp_sbio(sid) = nan;

sid = (md_sbio(:)<gr_bio_m(:));
md_sbio(sid) = nan;

sid = (lp_sbio(:)<gr_bio_l(:));
lp_sbio(sid) = nan;

sid = (ld_sbio(:)<gr_bio_l(:));
ld_sbio(sid) = nan;


%% v1: d =  g/m2   /  g/m2/d
sf_res1 = sf_sbio ./ sf_in;
sp_res1 = sp_sbio ./ sp_in;
sd_res1 = sd_sbio ./ sd_in;
mf_res1 = mf_sbio ./ mf_in;
mp_res1 = mp_sbio ./ mp_in;
md_res1 = md_sbio ./ md_in;
lp_res1 = lp_sbio ./ lp_in;
ld_res1 = ld_sbio ./ ld_in;

sf_res1(isinf(sf_res1(:))) = NaN;
sp_res1(isinf(sp_res1(:))) = NaN;
sd_res1(isinf(sd_res1(:))) = NaN;
mf_res1(isinf(mf_res1(:))) = NaN;
mp_res1(isinf(mp_res1(:))) = NaN;
md_res1(isinf(md_res1(:))) = NaN;
lp_res1(isinf(lp_res1(:))) = NaN;
ld_res1(isinf(ld_res1(:))) = NaN;

%% v2
sf_res2 = sf_sbio ./ sf_out;
sp_res2 = sp_sbio ./ sp_out;
sd_res2 = sd_sbio ./ sd_out;
mf_res2 = mf_sbio ./ mf_out;
mp_res2 = mp_sbio ./ mp_out;
md_res2 = md_sbio ./ md_out;
lp_res2 = lp_sbio ./ lp_out;
ld_res2 = ld_sbio ./ ld_out;

sf_res2(isinf(sf_res2(:))) = NaN;
sp_res2(isinf(sp_res2(:))) = NaN;
sd_res2(isinf(sd_res2(:))) = NaN;
mf_res2(isinf(mf_res2(:))) = NaN;
mp_res2(isinf(mp_res2(:))) = NaN;
md_res2(isinf(md_res2(:))) = NaN;
lp_res2(isinf(lp_res2(:))) = NaN;
ld_res2(isinf(ld_res2(:))) = NaN;

%% means (if using annual values)
sf_mbio = mean(sf_bio,2,'omitnan');
sp_mbio = mean(sp_bio,2,'omitnan');
sd_mbio = mean(sd_bio,2,'omitnan');
mf_mbio = mean(mf_bio,2,'omitnan');
mp_mbio = mean(mp_bio,2,'omitnan');
md_mbio = mean(md_bio,2,'omitnan');
lp_mbio = mean(lp_bio,2,'omitnan');
ld_mbio = mean(ld_bio,2,'omitnan');

% sf_min = mean(sf_in,2);
% sp_min = mean(sp_in,2);
% sd_min = mean(sd_in,2);
% mf_min = mean(mf_in,2);
% mp_min = mean(mp_in,2);
% md_min = mean(md_in,2);
% lp_min = mean(lp_in,2);
% ld_min = mean(ld_in,2);
%
% sf_mout = mean(sf_out,2);
% sp_mout = mean(sp_out,2);
% sd_mout = mean(sd_out,2);
% mf_mout = mean(mf_out,2);
% mp_mout = mean(mp_out,2);
% md_mout = mean(md_out,2);
% lp_mout = mean(lp_out,2);
% ld_mout = mean(ld_out,2);

sf_mres1 = mean(sf_res1,2,'omitnan');
sp_mres1 = mean(sp_res1,2,'omitnan');
sd_mres1 = mean(sd_res1,2,'omitnan');
mf_mres1 = mean(mf_res1,2,'omitnan');
mp_mres1 = mean(mp_res1,2,'omitnan');
md_mres1 = mean(md_res1,2,'omitnan');
lp_mres1 = mean(lp_res1,2,'omitnan');
ld_mres1 = mean(ld_res1,2,'omitnan');

sf_mres2 = mean(sf_res2,2,'omitnan');
sp_mres2 = mean(sp_res2,2,'omitnan');
sd_mres2 = mean(sd_res2,2,'omitnan');
mf_mres2 = mean(mf_res2,2,'omitnan');
mp_mres2 = mean(mp_res2,2,'omitnan');
md_mres2 = mean(md_res2,2,'omitnan');
lp_mres2 = mean(lp_res2,2,'omitnan');
ld_mres2 = mean(ld_res2,2,'omitnan');

%% means of res1 and res2

sf_mres = (sf_mres1 + sf_mres2) ./2;
sp_mres = (sp_mres1 + sp_mres2) ./2;
sd_mres = (sd_mres1 + sd_mres2) ./2;
mf_mres = (mf_mres1 + mf_mres2) ./2;
mp_mres = (mp_mres1 + mp_mres2) ./2;
md_mres = (md_mres1 + md_mres2) ./2;
lp_mres = (lp_mres1 + lp_mres2) ./2;
ld_mres = (ld_mres1 + ld_mres2) ./2;

%% Save
save([fpath 'Residence_time_momeans_FOSI_' harv cfile '.mat'],...
    'sf_mres','sp_mres','sd_mres','mf_mres','mp_mres','md_mres','lp_mres','ld_mres',...
    'sf_mres1','sp_mres1','sd_mres1','mf_mres1','mp_mres1','md_mres1','lp_mres1','ld_mres1',...
    'sf_mres2','sp_mres2','sd_mres2','mf_mres2','mp_mres2','md_mres2','lp_mres2','ld_mres2')

%% Histograms
figure(1)
subplot(3,3,1)
histogram(log10(sf_res1))
title('SF')

subplot(3,3,2)
histogram(log10(sp_res1))
title('SP')

subplot(3,3,3)
histogram(log10(sd_res1))
title('SD')

subplot(3,3,4)
histogram(log10(mf_res1))
title('MF')

subplot(3,3,5)
histogram(log10(mp_res1))
title('MP')

subplot(3,3,6)
histogram(log10(md_res1))
title('MD')

subplot(3,3,8)
histogram(log10(lp_res1))
title('LP')

subplot(3,3,9)
histogram(log10(ld_res1))
title('LD')

figure(2)
subplot(3,3,1)
histogram(log10(sf_res2))
title('SF')

subplot(3,3,2)
histogram(log10(sp_res2))
title('SP')

subplot(3,3,3)
histogram(log10(sd_res2))
title('SD')

subplot(3,3,4)
histogram(log10(mf_res2))
title('MF')

subplot(3,3,5)
histogram(log10(mp_res2))
title('MP')

subplot(3,3,6)
histogram(log10(md_res2))
title('MD')

subplot(3,3,8)
histogram(log10(lp_res2))
title('LP')

subplot(3,3,9)
histogram(log10(ld_res2))
title('LD')

%%
quantile(sf_mres(:),[0.1 0.25 0.5 0.75 0.9])
quantile(sf_mres(:),[0.9 0.95 0.99])

quantile(sp_mres(:),[0.9 0.95 0.99])
quantile(sd_mres(:),[0.9 0.95 0.99])

%% Maps
% plot info
geolon_t = double(TLONG);
geolat_t = double(TLAT);
plotminlat=-90;
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%% Plots in space
%bio
Bsf=NaN*ones(ni,nj);
Bsp=NaN*ones(ni,nj);
Bsd=NaN*ones(ni,nj);
Bmf=NaN*ones(ni,nj);
Bmp=NaN*ones(ni,nj);
Bmd=NaN*ones(ni,nj);
Blp=NaN*ones(ni,nj);
Bld=NaN*ones(ni,nj);
%In = res1
Isf=NaN*ones(ni,nj);
Isp=NaN*ones(ni,nj);
Isd=NaN*ones(ni,nj);
Imf=NaN*ones(ni,nj);
Imp=NaN*ones(ni,nj);
Imd=NaN*ones(ni,nj);
Ilp=NaN*ones(ni,nj);
Ild=NaN*ones(ni,nj);
%Out = res2
Osf=NaN*ones(ni,nj);
Osp=NaN*ones(ni,nj);
Osd=NaN*ones(ni,nj);
Omf=NaN*ones(ni,nj);
Omp=NaN*ones(ni,nj);
Omd=NaN*ones(ni,nj);
Olp=NaN*ones(ni,nj);
Old=NaN*ones(ni,nj);
%Mean res
Rsf=NaN*ones(ni,nj);
Rsp=NaN*ones(ni,nj);
Rsd=NaN*ones(ni,nj);
Rmf=NaN*ones(ni,nj);
Rmp=NaN*ones(ni,nj);
Rmd=NaN*ones(ni,nj);
Rlp=NaN*ones(ni,nj);
Rld=NaN*ones(ni,nj);

%bio
Bsf(ID)=sf_sbio;
Bsp(ID)=sp_sbio;
Bsd(ID)=sd_sbio;
Bmf(ID)=mf_sbio;
Bmp(ID)=mp_sbio;
Bmd(ID)=md_sbio;
Blp(ID)=lp_sbio;
Bld(ID)=ld_sbio;
%in
Isf(ID)=sf_res1;
Isp(ID)=sp_res1;
Isd(ID)=sd_res1;
Imf(ID)=mf_res1;
Imp(ID)=mp_res1;
Imd(ID)=md_res1;
Ilp(ID)=lp_res1;
Ild(ID)=ld_res1;
%out
Osf(ID)=sf_res2;
Osp(ID)=sp_res2;
Osd(ID)=sd_res2;
Omf(ID)=mf_res2;
Omp(ID)=mp_res2;
Omd(ID)=md_res2;
Olp(ID)=lp_res2;
Old(ID)=ld_res2;
%mres
Rsf(ID)=sf_mres;
Rsp(ID)=sp_mres;
Rsd(ID)=sd_mres;
Rmf(ID)=mf_mres;
Rmp(ID)=mp_mres;
Rmd(ID)=md_mres;
Rlp(ID)=lp_mres;
Rld(ID)=ld_mres;

%% 8 plot of bio
f3 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bsf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'SF bio','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bsp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'SP bio','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bsd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'SD bio','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bmf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'MF bio','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bmp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'MP bio','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bmd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'log_1_0 g m^-^2')
set(gcf,'renderer','painters')
text(0,1.75,'MD bio','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Blp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'LP bio','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Bld))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 2])
set(gcf,'renderer','painters')
text(0,1.75,'LD bio','HorizontalAlignment','center')
print('-dpng',[ppath 'FOSI_map_momean_bio_stages.png'])

%% 8 plot of flux in (in)
f4 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Isf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'SF res in','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Isp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'SP res in','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Isd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'SD res in','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Imf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'MF res in','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Imp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'MP res in','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Imd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'log_1_0 d')
set(gcf,'renderer','painters')
text(0,1.75,'MD res in','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Ilp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'LP res in','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Ild))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'LD res in','HorizontalAlignment','center')
print('-dpng',[ppath 'FOSI_map_momean_resIn_stages.png'])

%% 8 plot of flux out
f5 = figure('Units','inches','Position',[1 3 6.5 8]);

%A
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Osf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'SF res out','HorizontalAlignment','center')

%B
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Osp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'SP res out','HorizontalAlignment','center')

%C
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Osd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'SD res out','HorizontalAlignment','center')

%D
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Omf))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'MF res out','HorizontalAlignment','center')

%E
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Omp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'MP res out','HorizontalAlignment','center')

%F
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Omd))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
cb = colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
xlabel(cb,'log_1_0 d')
set(gcf,'renderer','painters')
text(0,1.75,'MD res out','HorizontalAlignment','center')

%G
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Olp))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'LP res out','HorizontalAlignment','center')

%H
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Old))
cmocean('matter')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-2 3])
set(gcf,'renderer','painters')
text(0,1.75,'LD res out','HorizontalAlignment','center')
print('-dpng',[ppath 'FOSI_map_momean_resOut_stages.png'])

%% 8 plot of mean res
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
print('-dpng',[ppath 'FOSI_map_momean_resMean_stages.png'])



