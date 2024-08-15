% CESM FEISTY FOSI runs
% calc interann variability by grid cell and lme

clear all
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v15_All_fish03_';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([fpath 'Annual_Means_FOSI_' mod cfile '.mat']);

% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
%cpath='/Volumes/petrik-lab/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

[ni,nj]=size(TLONG);
[nid,nt]=size(sf_abio);

%% Groups
xF = sf_abio + mf_abio;
xP = sp_abio + mp_abio + lp_abio;
xD = sd_abio + md_abio + ld_abio;
xS = sf_abio + sp_abio + sd_abio;
xM = mf_abio + mp_abio + md_abio;
xL = lp_abio + ld_abio;
xB = b_abio;
xall = xF + xP + xD;

%% mean & std by grid cell
msf = nanmean(sf_abio,2);
msp = nanmean(sp_abio,2);
msd = nanmean(sd_abio,2);
mmf = nanmean(mf_abio,2);
mmp = nanmean(mp_abio,2);
mmd = nanmean(md_abio,2);
mlp = nanmean(lp_abio,2);
mld = nanmean(ld_abio,2);
mb = nanmean(xB,2);
ma = nanmean(xall,2);
ms = nanmean(xS,2);
mm = nanmean(xM,2);
ml = nanmean(xL,2);
mf = nanmean(xF,2);
mp = nanmean(xP,2);
md = nanmean(xD,2);

ssf = std(sf_abio,0,2,'omitnan');
ssp = std(sp_abio,0,2,'omitnan');
ssd = std(sd_abio,0,2,'omitnan');
smf = std(mf_abio,0,2,'omitnan');
smp = std(mp_abio,0,2,'omitnan');
smd = std(md_abio,0,2,'omitnan');
slp = std(lp_abio,0,2,'omitnan');
sld = std(ld_abio,0,2,'omitnan');
sb = std(xB,0,2,'omitnan');
sa = std(xall,0,2,'omitnan');
ss = std(xS,0,2,'omitnan');
sm = std(xM,0,2,'omitnan');
sl = std(xL,0,2,'omitnan');
sf = std(xF,0,2,'omitnan');
sp = std(xP,0,2,'omitnan');
sd = std(xD,0,2,'omitnan');

%% mean & std by lme
%cpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'LME-mask-POP_gx1v6.mat']);

tlme = double(lme_mask);
tlme(tlme<0) = nan;
olme = tlme(GRD.ID);

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

for L=1:66
    lid = find(olme==L);

    lme_sf_std(L,1) = nanmean(std(sf_abio(lid,:),0,2,'omitnan'));
    lme_sp_std(L,1) = nanmean(std(sp_abio(lid,:),0,2,'omitnan'));
    lme_sd_std(L,1) = nanmean(std(sd_abio(lid,:),0,2,'omitnan'));
    lme_mf_std(L,1) = nanmean(std(mf_abio(lid,:),0,2,'omitnan'));
    lme_mp_std(L,1) = nanmean(std(mp_abio(lid,:),0,2,'omitnan'));
    lme_md_std(L,1) = nanmean(std(md_abio(lid,:),0,2,'omitnan'));
    lme_lp_std(L,1) = nanmean(std(lp_abio(lid,:),0,2,'omitnan'));
    lme_ld_std(L,1) = nanmean(std(ld_abio(lid,:),0,2,'omitnan'));
    lme_b_std(L,1) = nanmean(std(xB(lid,:),0,2,'omitnan'));
    lme_a_std(L,1) = nanmean(std(xall(lid,:),0,2,'omitnan'));
    lme_s_std(L,1) = nanmean(std(xS(lid,:),0,2,'omitnan'));
    lme_m_std(L,1) = nanmean(std(xM(lid,:),0,2,'omitnan'));
    lme_l_std(L,1) = nanmean(std(xL(lid,:),0,2,'omitnan'));
    lme_f_std(L,1) = nanmean(std(xF(lid,:),0,2,'omitnan'));
    lme_p_std(L,1) = nanmean(std(xP(lid,:),0,2,'omitnan'));
    lme_d_std(L,1) = nanmean(std(xD(lid,:),0,2,'omitnan'));

    lme_sf_mean(L,1) = nanmean(nanmean(sf_abio(lid,:),2));
    lme_sp_mean(L,1) = nanmean(nanmean(sp_abio(lid,:),2));
    lme_sd_mean(L,1) = nanmean(nanmean(sd_abio(lid,:),2));
    lme_mf_mean(L,1) = nanmean(nanmean(mf_abio(lid,:),2));
    lme_mp_mean(L,1) = nanmean(nanmean(mp_abio(lid,:),2));
    lme_md_mean(L,1) = nanmean(nanmean(md_abio(lid,:),2));
    lme_lp_mean(L,1) = nanmean(nanmean(lp_abio(lid,:),2));
    lme_ld_mean(L,1) = nanmean(nanmean(ld_abio(lid,:),2));
    lme_b_mean(L,1) = nanmean(nanmean(xB(lid,:),2));
    lme_a_mean(L,1) = nanmean(nanmean(xall(lid,:),2));
    lme_s_mean(L,1) = nanmean(nanmean(xS(lid,:),2));
    lme_m_mean(L,1) = nanmean(nanmean(xM(lid,:),2));
    lme_l_mean(L,1) = nanmean(nanmean(xL(lid,:),2));
    lme_f_mean(L,1) = nanmean(nanmean(xF(lid,:),2));
    lme_p_mean(L,1) = nanmean(nanmean(xP(lid,:),2));
    lme_d_mean(L,1) = nanmean(nanmean(xD(lid,:),2));

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

%% ANOMALIES -------------------------------------------------

%% remove linear trend
sf = NaN*ones(size(sf_abio));
sp = NaN*ones(size(sp_abio));
sd = NaN*ones(size(sd_abio));
mf = NaN*ones(size(mf_abio));
mp = NaN*ones(size(mp_abio));
md = NaN*ones(size(md_abio));
lp = NaN*ones(size(lp_abio));
ld = NaN*ones(size(ld_abio));

B = NaN*ones(size(xB));
F = NaN*ones(size(xF));
P = NaN*ones(size(xP));
D = NaN*ones(size(xD));
S = NaN*ones(size(xS));
M = NaN*ones(size(xM));
L = NaN*ones(size(xL));
All = NaN*ones(size(xall));

for i = 1:nid
    %STAGES
    xi = sf_abio(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    sf(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = sp_abio(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    sp(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = sd_abio(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    sd(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = mf_abio(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    mf(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = mp_abio(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    mp(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = md_abio(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    md(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = lp_abio(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    lp(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = ld_abio(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    ld(i,:) = dR;
    clear R T t b m tH dR data
    
    %TYPES
    xi = xB(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    B(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = xF(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    F(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = xP(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    P(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = xD(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    D(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = xS(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    S(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = xM(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    M(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = xL(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    L(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = xall(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    All(i,:) = dR;
    clear R T t b m tH dR data
    
end

%% anomalies 
asf = sf - nanmean(sf,2);
asp = sp - nanmean(sp,2);
asd = sd - nanmean(sd,2);
amf = mf - nanmean(mf,2);
amp = mp - nanmean(mp,2);
amd = md - nanmean(md,2);
alp = lp - nanmean(lp,2);
ald = ld - nanmean(ld,2);
aa = All - nanmean(All,2);
as = S - nanmean(S,2);
am = M - nanmean(M,2);
al = L - nanmean(L,2);
af = F - nanmean(F,2);
ap = P - nanmean(P,2);
ad = D - nanmean(D,2);
ab = B - nanmean(B,2);

%% var of anomalies by grid cell
vsf = var(asf,0,2,'omitnan');
vsp = var(asp,0,2,'omitnan');
vsd = var(asd,0,2,'omitnan');
vmf = var(amf,0,2,'omitnan');
vmp = var(amp,0,2,'omitnan');
vmd = var(amd,0,2,'omitnan');
vlp = var(alp,0,2,'omitnan');
vld = var(ald,0,2,'omitnan');
vb = var(ab,0,2,'omitnan');
va = var(aa,0,2,'omitnan');
vs = var(as,0,2,'omitnan');
vm = var(am,0,2,'omitnan');
vl = var(al,0,2,'omitnan');
vf = var(af,0,2,'omitnan');
vp = var(ap,0,2,'omitnan');
vd = var(ad,0,2,'omitnan');

%% Anom on grid
Zsf=NaN*ones(ni*nj,nt);
Zsp=NaN*ones(ni*nj,nt);
Zsd=NaN*ones(ni*nj,nt);
Zmf=NaN*ones(ni*nj,nt);
Zmp=NaN*ones(ni*nj,nt);
Zmd=NaN*ones(ni*nj,nt);
Zlp=NaN*ones(ni*nj,nt);
Zld=NaN*ones(ni*nj,nt);
Zb=NaN*ones(ni*nj,nt);
Zf=NaN*ones(ni*nj,nt);
Zp=NaN*ones(ni*nj,nt);
Zd=NaN*ones(ni*nj,nt);
Zs=NaN*ones(ni*nj,nt);
Zm=NaN*ones(ni*nj,nt);
Zl=NaN*ones(ni*nj,nt);
Za=NaN*ones(ni*nj,nt);

Zsf(GRD.ID,:)=asf;
Zsp(GRD.ID,:)=asp;
Zsd(GRD.ID,:)=asd;
Zmf(GRD.ID,:)=amf;
Zmp(GRD.ID,:)=amp;
Zmd(GRD.ID,:)=amd;
Zlp(GRD.ID,:)=alp;
Zld(GRD.ID,:)=ald;
Zb(GRD.ID,:)=ab;
Zf(GRD.ID,:)=af;
Zp(GRD.ID,:)=ap;
Zd(GRD.ID,:)=ad;
Zs(GRD.ID,:)=as;
Zm(GRD.ID,:)=am;
Zl(GRD.ID,:)=al;
Za(GRD.ID,:)=aa;

Zsf=reshape(Zsf,ni,nj,nt);
Zsp=reshape(Zsp,ni,nj,nt);
Zsd=reshape(Zsd,ni,nj,nt);
Zmf=reshape(Zmf,ni,nj,nt);
Zmp=reshape(Zmp,ni,nj,nt);
Zmd=reshape(Zmd,ni,nj,nt);
Zlp=reshape(Zlp,ni,nj,nt);
Zld=reshape(Zld,ni,nj,nt);
Zb=reshape(Zb,ni,nj,nt);
Zf=reshape(Zf,ni,nj,nt);
Zp=reshape(Zp,ni,nj,nt);
Zd=reshape(Zd,ni,nj,nt);
Zs=reshape(Zs,ni,nj,nt);
Zm=reshape(Zm,ni,nj,nt);
Zl=reshape(Zl,ni,nj,nt);
Za=reshape(Za,ni,nj,nt);

%% map info
latlim=[-90 90];
lonlim=[-280 80];
load coastlines

cmYR=cbrewer('seq','YlOrRd',66,'pchip');
cmYOB=cbrewer('seq','YlOrBr',66,'pchip');
cmOR=cbrewer('seq','OrRd',66,'pchip');

%% reshape CV
cvsf=NaN*ones(ni,nj);
cvsp=NaN*ones(ni,nj);
cvsd=NaN*ones(ni,nj);
cvmf=NaN*ones(ni,nj);
cvmp=NaN*ones(ni,nj);
cvmd=NaN*ones(ni,nj);
cvlp=NaN*ones(ni,nj);
cvld=NaN*ones(ni,nj);
cvb=NaN*ones(ni,nj);
cva = NaN*ones(ni,nj);
cvs = NaN*ones(ni,nj);
cvm = NaN*ones(ni,nj);
cvl = NaN*ones(ni,nj);
cvf = NaN*ones(ni,nj);
cvp = NaN*ones(ni,nj);
cvd = NaN*ones(ni,nj);

cvsf(GRD.ID)=sf_cv;
cvsp(GRD.ID)=sp_cv;
cvsd(GRD.ID)=sd_cv;
cvmf(GRD.ID)=mf_cv;
cvmp(GRD.ID)=mp_cv;
cvmd(GRD.ID)=md_cv;
cvlp(GRD.ID)=lp_cv;
cvld(GRD.ID)=ld_cv;
cvb(GRD.ID)=b_cv;
cva(GRD.ID)=a_cv;
cvs(GRD.ID)=s_cv;
cvm(GRD.ID)=m_cv;
cvl(GRD.ID)=l_cv;
cvf(GRD.ID)=f_cv;
cvp(GRD.ID)=p_cv;
cvd(GRD.ID)=d_cv;

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

print('-dpng',[ppath 'Map_FEISTY_FOSI_',mod,'interann_coeffvar_stages.png'])

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

print('-dpng',[ppath 'Map_FEISTY_FOSI_',mod,'interann_coeffvar_types.png'])

%% save
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
%fpath=['/Volumes/petrik-lab/NC/CESM_MAPP/' cfile '/'];
save([fpath 'FEISTY_FOSI_',mod,'interann_var.mat'],...
    'cvsf','cvsp','cvsd','cvmf','cvmp','cvmd','cvlp','cvld','cvb','cva',...
    'cvs','cvm','cvl','cvf','cvp','cvd',...
    'msf','msp','msd','mmf','mmp','mmd','mlp','mld','mb','ma',...
    'ms','mm','ml','mf','mp','md',...
    'ssf','ssp','ssd','smf','smp','smd','slp','sld','sb','sa',...
    'ss','sm','sl','sf','sp','sd',...
    'lme_sf_std','lme_sp_std','lme_sd_std',...
    'lme_mf_std','lme_mp_std','lme_md_std',...
    'lme_lp_std','lme_ld_std','lme_b_std','lme_a_std',...
    'lme_s_std','lme_m_std','lme_l_std',...
    'lme_f_std','lme_p_std','lme_d_std',...
    'lme_sf_mean','lme_sp_mean','lme_sd_mean',...
    'lme_mf_mean','lme_mp_mean','lme_md_mean',...
    'lme_lp_mean','lme_ld_mean','lme_b_mean','lme_a_mean',...
    'lme_s_mean','lme_m_mean','lme_l_mean',...
    'lme_f_mean','lme_p_mean','lme_d_mean',...
    'lme_sf_cv','lme_sp_cv','lme_sd_cv',...
    'lme_mf_cv','lme_mp_cv','lme_md_cv',...
    'lme_lp_cv','lme_ld_cv','lme_b_cv','lme_a_cv',...
    'lme_s_cv','lme_m_cv','lme_l_cv',...
    'lme_f_cv','lme_p_cv','lme_d_cv');

%%
save([fpath 'FEISTY_FOSI_',mod,'ann_mean_anoms.mat'],...
    'asf','asp','asd','amf','amp','amd','alp','ald','ab','aa','as','am','al',...
    'af','ap','ad',...
    'vsf','vsp','vsd','vmf','vmp','vmd','vlp','vld','vb','va','vs','vm','vl',...
    'vf','vp','vd',...
    'Zsf','Zsp','Zsd','Zmf','Zmp','Zmd','Zlp','Zld','Zb','Za','Zs','Zm','Zl',...
    'Zf','Zp','Zd');
