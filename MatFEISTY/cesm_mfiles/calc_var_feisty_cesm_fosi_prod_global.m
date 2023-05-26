% CESM FEISTY FOSI runs
% calc interann variability by grid cell of prod

clear 
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
sims = {'v15_All_fish03_';'v15_climatol_';'v15_varTemp_';'v15_varFood_'};
mod = sims{1};

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
ppath = [pp cfile '/FOSI/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([fpath 'Annual_Means_FOSI_' mod cfile '.mat'],...
    'sf_aprod','sp_aprod','sd_aprod',...
    'mf_aprod','mp_aprod','md_aprod',...
    'lp_aprod','ld_aprod');

% Map data
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

[ni,nj]=size(TLONG);
[nid,nt]=size(sf_aprod);

%% Groups
xF = mf_aprod;
xP = lp_aprod;
xD = ld_aprod;
xS = sf_aprod + sp_aprod + sd_aprod;
xM = mf_aprod + mp_aprod + md_aprod;
xL = lp_aprod + ld_aprod;
xall = xF + xP + xD;

%% mean & std by grid cell
msf = mean(sf_aprod,2,'omitnan');
msp = mean(sp_aprod,2,'omitnan');
msd = mean(sd_aprod,2,'omitnan');
mmf = mean(mf_aprod,2,'omitnan');
mmp = mean(mp_aprod,2,'omitnan');
mmd = mean(md_aprod,2,'omitnan');
mlp = mean(lp_aprod,2,'omitnan');
mld = mean(ld_aprod,2,'omitnan');
ma = mean(xall,2,'omitnan');
ms = mean(xS,2,'omitnan');
mm = mean(xM,2,'omitnan');
ml = mean(xL,2,'omitnan');
mf = mean(xF,2,'omitnan');
mp = mean(xP,2,'omitnan');
md = mean(xD,2,'omitnan');

ssf = std(sf_aprod,0,2,'omitnan');
ssp = std(sp_aprod,0,2,'omitnan');
ssd = std(sd_aprod,0,2,'omitnan');
smf = std(mf_aprod,0,2,'omitnan');
smp = std(mp_aprod,0,2,'omitnan');
smd = std(md_aprod,0,2,'omitnan');
slp = std(lp_aprod,0,2,'omitnan');
sld = std(ld_aprod,0,2,'omitnan');
sa = std(xall,0,2,'omitnan');
ss = std(xS,0,2,'omitnan');
sm = std(xM,0,2,'omitnan');
sl = std(xL,0,2,'omitnan');
sf = std(xF,0,2,'omitnan');
sp = std(xP,0,2,'omitnan');
sd = std(xD,0,2,'omitnan');

%% Coefficient of variance
sf_cv = ssf ./ msf;
sp_cv = ssp ./ msp;
sd_cv = ssd ./ msd;
mf_cv = smf ./ mmf;
mp_cv = smp ./ mmp;
md_cv = smd ./ mmd;
lp_cv = slp ./ mlp;
ld_cv = sld ./ mld;
a_cv = sa ./ ma;
s_cv = ss ./ ms;
m_cv = sm ./ mm;
l_cv = sl ./ ml;
f_cv = sf ./ mf;
p_cv = sp ./ mp;
d_cv = sd ./ md;

%% ANOMALIES -------------------------------------------------

%% remove linear trend
sf = NaN*ones(size(sf_aprod));
sp = NaN*ones(size(sp_aprod));
sd = NaN*ones(size(sd_aprod));
mf = NaN*ones(size(mf_aprod));
mp = NaN*ones(size(mp_aprod));
md = NaN*ones(size(md_aprod));
lp = NaN*ones(size(lp_aprod));
ld = NaN*ones(size(ld_aprod));

S = NaN*ones(size(xS));
M = NaN*ones(size(xM));
L = NaN*ones(size(xL));
All = NaN*ones(size(xall));

for i = 1:nid
    %STAGES
    xi = sf_aprod(i,:);
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
    
    xi = sp_aprod(i,:);
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
    
    xi = sd_aprod(i,:);
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
    
    xi = mf_aprod(i,:);
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
    
    xi = mp_aprod(i,:);
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
    
    xi = md_aprod(i,:);
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
    
    xi = lp_aprod(i,:);
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
    
    xi = ld_aprod(i,:);
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
asf = sf - mean(sf,2,'omitnan');
asp = sp - mean(sp,2,'omitnan');
asd = sd - mean(sd,2,'omitnan');
amf = mf - mean(mf,2,'omitnan');
amp = mp - mean(mp,2,'omitnan');
amd = md - mean(md,2,'omitnan');
alp = lp - mean(lp,2,'omitnan');
ald = ld - mean(ld,2,'omitnan');
aa = All - mean(All,2,'omitnan');
as = S - mean(S,2,'omitnan');
am = M - mean(M,2,'omitnan');
al = L - mean(L,2,'omitnan');

%% var of anomalies by grid cell
vsf = var(asf,0,2,'omitnan');
vsp = var(asp,0,2,'omitnan');
vsd = var(asd,0,2,'omitnan');
vmf = var(amf,0,2,'omitnan');
vmp = var(amp,0,2,'omitnan');
vmd = var(amd,0,2,'omitnan');
vlp = var(alp,0,2,'omitnan');
vld = var(ald,0,2,'omitnan');
va = var(aa,0,2,'omitnan');
vs = var(as,0,2,'omitnan');
vm = var(am,0,2,'omitnan');
vl = var(al,0,2,'omitnan');

%% Anom on grid
Zsf=NaN*ones(ni*nj,nt);
Zsp=NaN*ones(ni*nj,nt);
Zsd=NaN*ones(ni*nj,nt);
Zmf=NaN*ones(ni*nj,nt);
Zmp=NaN*ones(ni*nj,nt);
Zmd=NaN*ones(ni*nj,nt);
Zlp=NaN*ones(ni*nj,nt);
Zld=NaN*ones(ni*nj,nt);
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
cva(GRD.ID)=a_cv;
cvs(GRD.ID)=s_cv;
cvm(GRD.ID)=m_cv;
cvl(GRD.ID)=l_cv;
cvf(GRD.ID)=f_cv;
cvp(GRD.ID)=p_cv;
cvd(GRD.ID)=d_cv;

%% save
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
save([fpath 'FEISTY_FOSI_',mod,'prod_interann_var.mat'],...
    'cvsf','cvsp','cvsd','cvmf','cvmp','cvmd','cvlp','cvld','cva',...
    'cvs','cvm','cvl','cvf','cvp','cvd',...
    'msf','msp','msd','mmf','mmp','mmd','mlp','mld','ma',...
    'ms','mm','ml','mf','mp','md',...
    'ssf','ssp','ssd','smf','smp','smd','slp','sld','sa',...
    'ss','sm','sl','sf','sp','sd');

%%
save([fpath 'FEISTY_FOSI_',mod,'prod_ann_mean_anoms.mat'],...
    'asf','asp','asd','amf','amp','amd','alp','ald','aa','as','am','al',...
    'vsf','vsp','vsd','vmf','vmp','vmd','vlp','vld','va','vs','vm','vl',...
    'Zsf','Zsp','Zsd','Zmf','Zmp','Zmd','Zlp','Zld','Za','Zs','Zm','Zl');
