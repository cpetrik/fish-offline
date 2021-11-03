% CESM FEISTY FOSI runs
% calc interann variability by grid cell and lme

clear all
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_noCC_RE00100';
mod = 'v13_All_fish03_';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([fpath 'Annual_Means_FOSI_' mod cfile '.mat']);

% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

[ni,nj]=size(TLONG);
[nid,nt]=size(sf_abio);

%% Biomass on grid
Zsf=NaN*ones(ni*nj,nt);
Zsp=NaN*ones(ni*nj,nt);
Zsd=NaN*ones(ni*nj,nt);
Zmf=NaN*ones(ni*nj,nt);
Zmp=NaN*ones(ni*nj,nt);
Zmd=NaN*ones(ni*nj,nt);
Zlp=NaN*ones(ni*nj,nt);
Zld=NaN*ones(ni*nj,nt);
Zb=NaN*ones(ni*nj,nt);

Zsf(GRD.ID,:)=sf_abio;
Zsp(GRD.ID,:)=sp_abio;
Zsd(GRD.ID,:)=sd_abio;
Zmf(GRD.ID,:)=mf_abio;
Zmp(GRD.ID,:)=mp_abio;
Zmd(GRD.ID,:)=md_abio;
Zlp(GRD.ID,:)=lp_abio;
Zld(GRD.ID,:)=ld_abio;
Zb(GRD.ID,:)=b_abio;

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

%% anomalies
asf = Zsf - nanmean(Zsf,2);
asp = Zsp - nanmean(Zsp,2);
asd = Zsd - nanmean(Zsd,2);
amf = Zmf - nanmean(Zmf,2);
amp = Zmp - nanmean(Zmp,2);
amd = Zmd - nanmean(Zmd,2);
alp = Zlp - nanmean(Zlp,2);
ald = Zld - nanmean(Zld,2);
ab = Zb - nanmean(Zb,2);
aa = All - nanmean(All,2);
as = AllS - nanmean(AllS,2);
am = AllM - nanmean(AllM,2);
al = AllL - nanmean(AllL,2);
af = AllF - nanmean(AllF,2);
ap = AllP - nanmean(AllP,2);
ad = AllD - nanmean(AllD,2);
apf = FracPF - nanmean(FracPF,2);
apd = FracPD - nanmean(FracPD,2);
alm = FracLM - nanmean(FracLM,2);

% %% var of anomalies by grid cell
% vsf = var(asf,0,2,'omitnan');
% vsp = var(asp,0,2,'omitnan');
% vsd = var(asd,0,2,'omitnan');
% vmf = var(amf,0,2,'omitnan');
% vmp = var(amp,0,2,'omitnan');
% vmd = var(amd,0,2,'omitnan');
% vlp = var(alp,0,2,'omitnan');
% vld = var(ald,0,2,'omitnan');
% vb = var(ab,0,2,'omitnan');
% va = var(aa,0,2,'omitnan');
% vs = var(as,0,2,'omitnan');
% vm = var(am,0,2,'omitnan');
% vl = var(al,0,2,'omitnan');
% vf = var(af,0,2,'omitnan');
% vp = var(ap,0,2,'omitnan');
% vd = var(ad,0,2,'omitnan');
% vpf = var(apf,0,2,'omitnan');
% vpd = var(apd,0,2,'omitnan');
% vlm = var(alm,0,2,'omitnan');

%% mean & std by grid cell
msf = nanmean(Zsf,2);
msp = nanmean(Zsp,2);
msd = nanmean(Zsd,2);
mmf = nanmean(Zmf,2);
mmp = nanmean(Zmp,2);
mmd = nanmean(Zmd,2);
mlp = nanmean(Zlp,2);
mld = nanmean(Zld,2);
mb = nanmean(Zb,2);
ma = nanmean(All,2);
ms = nanmean(AllS,2);
mm = nanmean(AllM,2);
ml = nanmean(AllL,2);
mf = nanmean(AllF,2);
mp = nanmean(AllP,2);
md = nanmean(AllD,2);
mpf = nanmean(FracPF,2);
mpd = nanmean(FracPD,2);
mlm = nanmean(FracLM,2);
        
ssf = std(Zsf,0,2,'omitnan');
ssp = std(Zsp,0,2,'omitnan');
ssd = std(Zsd,0,2,'omitnan');
smf = std(Zmf,0,2,'omitnan');
smp = std(Zmp,0,2,'omitnan');
smd = std(Zmd,0,2,'omitnan');
slp = std(Zlp,0,2,'omitnan');
sld = std(Zld,0,2,'omitnan');
sb = std(Zb,0,2,'omitnan');
sa = std(All,0,2,'omitnan');
ss = std(AllS,0,2,'omitnan');
sm = std(AllM,0,2,'omitnan');
sl = std(AllL,0,2,'omitnan');
sf = std(AllF,0,2,'omitnan');
sp = std(AllP,0,2,'omitnan');
sd = std(AllD,0,2,'omitnan');
spf = std(FracPF,0,2,'omitnan');
spd = std(FracPD,0,2,'omitnan');
slm = std(FracLM,0,2,'omitnan');

%% mean & std by lme
cpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'LME-mask-POP_gx1v6.mat']);

tlme = double(lme_mask);
tlme(tlme<0) = nan;

lme_sf_mean = NaN*ones(66,1);
lme_sp_mean = NaN*ones(66,1);
lme_sd_mean = NaN*ones(66,1);
lme_mf_mean = NaN*ones(66,1);
lme_mp_mean = NaN*ones(66,1);
lme_md_mean = NaN*ones(66,1);
lme_lp_mean = NaN*ones(66,1);
lme_ld_mean = NaN*ones(66,1);
lme_b_mean = NaN*ones(66,1);
lme_a_mean = NaN*ones(66,1);
lme_s_mean = NaN*ones(66,1);
lme_m_mean = NaN*ones(66,1);
lme_l_mean = NaN*ones(66,1);
lme_f_mean = NaN*ones(66,1);
lme_p_mean = NaN*ones(66,1);
lme_d_mean = NaN*ones(66,1);
lme_pf_mean = NaN*ones(66,1);
lme_pd_mean = NaN*ones(66,1);
lme_lm_mean = NaN*ones(66,1);
lme_sf_std = NaN*ones(66,1);
lme_sp_std = NaN*ones(66,1);
lme_sd_std = NaN*ones(66,1);
lme_mf_std = NaN*ones(66,1);
lme_mp_std = NaN*ones(66,1);
lme_md_std = NaN*ones(66,1);
lme_lp_std = NaN*ones(66,1);
lme_ld_std = NaN*ones(66,1);
lme_b_std = NaN*ones(66,1);
lme_a_std = NaN*ones(66,1);
lme_s_std = NaN*ones(66,1);
lme_m_std = NaN*ones(66,1);
lme_l_std = NaN*ones(66,1);
lme_f_std = NaN*ones(66,1);
lme_p_std = NaN*ones(66,1);
lme_d_std = NaN*ones(66,1);
lme_pf_std = NaN*ones(66,1);
lme_pd_std = NaN*ones(66,1);
lme_lm_std = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    
    lme_sf_std(L,1) = nanmean(std(asf(lid,:),0,2,'omitnan'));
    lme_sp_std(L,1) = nanmean(std(asp(lid,:),0,2,'omitnan'));
    lme_sd_std(L,1) = nanmean(std(asd(lid,:),0,2,'omitnan'));
    lme_mf_std(L,1) = nanmean(std(amf(lid,:),0,2,'omitnan'));
    lme_mp_std(L,1) = nanmean(std(amp(lid,:),0,2,'omitnan'));
    lme_md_std(L,1) = nanmean(std(amd(lid,:),0,2,'omitnan'));
    lme_lp_std(L,1) = nanmean(std(alp(lid,:),0,2,'omitnan'));
    lme_ld_std(L,1) = nanmean(std(ald(lid,:),0,2,'omitnan'));
    lme_b_std(L,1) = nanmean(std(ab(lid,:),0,2,'omitnan'));
    lme_a_std(L,1) = nanmean(std(aa(lid,:),0,2,'omitnan'));
    lme_s_std(L,1) = nanmean(std(as(lid,:),0,2,'omitnan'));
    lme_m_std(L,1) = nanmean(std(am(lid,:),0,2,'omitnan'));
    lme_l_std(L,1) = nanmean(std(al(lid,:),0,2,'omitnan'));
    lme_f_std(L,1) = nanmean(std(af(lid,:),0,2,'omitnan'));
    lme_p_std(L,1) = nanmean(std(ap(lid,:),0,2,'omitnan'));
    lme_d_std(L,1) = nanmean(std(ad(lid,:),0,2,'omitnan'));
    lme_pf_std(L,1) = nanmean(std(apf(lid,:),0,2,'omitnan'));
    lme_pd_std(L,1) = nanmean(std(apd(lid,:),0,2,'omitnan'));
    lme_lm_std(L,1) = nanmean(std(alm(lid,:),0,2,'omitnan'));
    
    lme_sf_mean(L,1) = nanmean(nanmean(asf(lid,:),2));
    lme_sp_mean(L,1) = nanmean(nanmean(asp(lid,:),2));
    lme_sd_mean(L,1) = nanmean(nanmean(asd(lid,:),2));
    lme_mf_mean(L,1) = nanmean(nanmean(amf(lid,:),2));
    lme_mp_mean(L,1) = nanmean(nanmean(amp(lid,:),2));
    lme_md_mean(L,1) = nanmean(nanmean(amd(lid,:),2));
    lme_lp_mean(L,1) = nanmean(nanmean(alp(lid,:),2));
    lme_ld_mean(L,1) = nanmean(nanmean(ald(lid,:),2));
    lme_b_mean(L,1) = nanmean(nanmean(ab(lid,:),2));
    lme_a_mean(L,1) = nanmean(nanmean(aa(lid,:),2));
    lme_s_mean(L,1) = nanmean(nanmean(as(lid,:),2));
    lme_m_mean(L,1) = nanmean(nanmean(am(lid,:),2));
    lme_l_mean(L,1) = nanmean(nanmean(al(lid,:),2));
    lme_f_mean(L,1) = nanmean(nanmean(af(lid,:),2));
    lme_p_mean(L,1) = nanmean(nanmean(ap(lid,:),2));
    lme_d_mean(L,1) = nanmean(nanmean(ad(lid,:),2));
    lme_pf_mean(L,1) = nanmean(nanmean(apf(lid,:),2));
    lme_pd_mean(L,1) = nanmean(nanmean(apd(lid,:),2));
    lme_lm_mean(L,1) = nanmean(nanmean(alm(lid,:),2));
    
end

%% Coefficient of variance
sf_cv = ssf ./ msf;
sp_cv = ssp ./ msp;
sd_cv = ssd ./ msd;
mf_cv = smf ./ mmf;
mp_cv = smp ./ mmp;
md_cv = smd ./ mmd;
lp_cv = slp ./ mlp;
ld_cv = sld ./ mld;
b_cv = sb ./ mb;
a_cv = sa ./ ma;
s_cv = ss ./ ms;
m_cv = sm ./ mm;
l_cv = sl ./ ml;
f_cv = sf ./ mf;
p_cv = sp ./ mp;
d_cv = sd ./ md;
pf_cv = spf ./ mpf;
pd_cv = spd ./ mpd;
lm_cv = slm ./ mlm;

lme_sf_cv = lme_sf_std ./ lme_sf_mean;
lme_sp_cv = lme_sp_std ./ lme_sp_mean;
lme_sd_cv = lme_sd_std ./ lme_sd_mean;
lme_mf_cv = lme_mf_std ./ lme_mf_mean;
lme_mp_cv = lme_mp_std ./ lme_mp_mean;
lme_md_cv = lme_md_std ./ lme_md_mean;
lme_lp_cv = lme_lp_std ./ lme_lp_mean;
lme_ld_cv = lme_ld_std ./ lme_ld_mean;
lme_b_cv = lme_b_std ./ lme_b_mean;
lme_a_cv = lme_a_std ./ lme_a_mean;
lme_s_cv = lme_s_std ./ lme_s_mean;
lme_m_cv = lme_m_std ./ lme_m_mean;
lme_l_cv = lme_l_std ./ lme_l_mean;
lme_f_cv = lme_f_std ./ lme_f_mean;
lme_p_cv = lme_p_std ./ lme_p_mean;
lme_d_cv = lme_d_std ./ lme_d_mean;
lme_pf_cv = lme_pf_std ./ lme_pf_mean;
lme_pd_cv = lme_pd_std ./ lme_pd_mean;
lme_lm_cv = lme_lm_std ./ lme_lm_mean;

%% map info
latlim=[-90 90];
lonlim=[-280 80];
load coastlines

cmYR=cbrewer('seq','YlOrRd',66,'pchip');
cmYOB=cbrewer('seq','YlOrBr',66,'pchip');
cmOR=cbrewer('seq','OrRd',66,'pchip');

%% reshape
cvsf = reshape(sf_cv,ni,nj);
cvsp = reshape(sp_cv,ni,nj);
cvsd = reshape(sd_cv,ni,nj);
cvmf = reshape(mf_cv,ni,nj);
cvmp = reshape(mp_cv,ni,nj);
cvmd = reshape(md_cv,ni,nj);
cvlp = reshape(lp_cv,ni,nj);
cvld = reshape(ld_cv,ni,nj);
cvb = reshape(b_cv,ni,nj);
cva = reshape(a_cv,ni,nj);
cvs = reshape(s_cv,ni,nj);
cvm = reshape(m_cv,ni,nj);
cvl = reshape(l_cv,ni,nj);
cvf = reshape(f_cv,ni,nj);
cvp = reshape(p_cv,ni,nj);
cvd = reshape(d_cv,ni,nj);
cvpf = reshape(pf_cv,ni,nj);
cvpd = reshape(pd_cv,ni,nj);
cvlm = reshape(lm_cv,ni,nj);

%% map
% 8plot by stage 
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - sf
subplot('Position',[0.025 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvsf)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'SF','HorizontalAlignment','center')

%B - sp
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvsp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'SP','HorizontalAlignment','center')

%C - mp
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvmp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'MP','HorizontalAlignment','center')

%D - lp
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvlp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'LP','HorizontalAlignment','center')

%E - mf
subplot('Position',[0.5 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvmf)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'MF','HorizontalAlignment','center')

%F - sd
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvsd)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
%colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'SD','HorizontalAlignment','center')

%G - md
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvmd)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'MD','HorizontalAlignment','center')

%H - ld
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvld)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'LD','HorizontalAlignment','center')

print('-dpng',[pp 'Map_FEISTY_FOSI_',mod,'interann_coeffvar_stages.png'])

%% 8plot by type 
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - s
subplot('Position',[0.025 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvs)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'Sm','HorizontalAlignment','center')

%B - m
subplot('Position',[0.025 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvm)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'Md','HorizontalAlignment','center')

%C - l
subplot('Position',[0.025 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvl)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'Lg','HorizontalAlignment','center')

%D - F
subplot('Position',[0.025 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvf)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'F','HorizontalAlignment','center')

%E - P
subplot('Position',[0.5 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvp)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'P','HorizontalAlignment','center')

%F - D
subplot('Position',[0.5 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvd)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
%colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'D','HorizontalAlignment','center')

%G - B
subplot('Position',[0.5 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cvb)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'Bent','HorizontalAlignment','center')

%H - all
subplot('Position',[0.5 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT,TLONG,cva)
colormap(cmYR)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'All','HorizontalAlignment','center')

print('-dpng',[pp 'Map_FEISTY_FOSI_',mod,'interann_coeffvar_types.png'])

%% save
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
save([fpath 'FEISTY_FOSI_',mod,'interann_var.mat'],...
    'cvsf','cvsp','cvsd','cvmf','cvmp','cvmd','cvlp','cvld','cvb','cva',...
    'cvs','cvm','cvl','cvf','cvp','cvd','cvpf','cvpd','cvlm',...
    'msf','msp','msd','mmf','mmp','mmd','mlp','mld','mb','ma',...
    'ms','mm','ml','mf','mp','md','mpf','mpd','mlm',...
    'ssf','ssp','ssd','smf','smp','smd','slp','sld','sb','sa',...
    'ss','sm','sl','sf','sp','sd','spf','spd','slm',...
    'lme_sf_std','lme_sp_std','lme_sd_std',...
    'lme_mf_std','lme_mp_std','lme_md_std',...
    'lme_lp_std','lme_ld_std','lme_b_std','lme_a_std',...
    'lme_s_std','lme_m_std','lme_l_std',...
    'lme_f_std','lme_p_std','lme_d_std',...
    'lme_pf_std','lme_pd_std','lme_lm_std',...
    'lme_sf_mean','lme_sp_mean','lme_sd_mean',...
    'lme_mf_mean','lme_mp_mean','lme_md_mean',...
    'lme_lp_mean','lme_ld_mean','lme_b_mean','lme_a_mean',...
    'lme_s_mean','lme_m_mean','lme_l_mean',...
    'lme_f_mean','lme_p_mean','lme_d_mean',...
    'lme_pf_mean','lme_pd_mean','lme_lm_mean',...
    'lme_sf_cv','lme_sp_cv','lme_sd_cv',...
    'lme_mf_cv','lme_mp_cv','lme_md_cv',...
    'lme_lp_cv','lme_ld_cv','lme_b_cv','lme_a_cv',...
    'lme_s_cv','lme_m_cv','lme_l_cv',...
    'lme_f_cv','lme_p_cv','lme_d_cv',...
    'lme_pf_cv','lme_pd_cv','lme_lm_cv');
save([fpath 'FEISTY_FOSI_',mod,'ann_mean_anoms.mat'],...
    'asf','asp','asd','amf','amp','amd','alp','ald','ab','aa','as','am','al',...
    'af','ap','ad','apf','apd','alm');

