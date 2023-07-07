% CESM FEISTY FOSI runs
% calc interann variability by grid cell of nu
% convert from per day to per year

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
    'mf_anu','lp_anu','ld_anu');

% Map data
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

[ni,nj]=size(TLONG);
[nid,nt]=size(mf_anu);

%% Per day to per year
mf_anu = 365 * mf_anu;
lp_anu = 365 * lp_anu;
ld_anu = 365 * ld_anu;

%% Groups
xF = mf_anu;
xP = lp_anu;
xD = ld_anu;

xA = (mf_anu+lp_anu+ld_anu)./3;

%% mean & std by grid cell
mmf = mean(xF,2,'omitnan');
mlp = mean(xP,2,'omitnan');
mld = mean(xD,2,'omitnan');
mall = mean(xA,2,'omitnan');

sf = std(xF,0,2,'omitnan');
sp = std(xP,0,2,'omitnan');
sd = std(xD,0,2,'omitnan');
sall = std(xA,0,2,'omitnan');

%% Coefficient of variance
f_cv = sf ./ mmf;
p_cv = sp ./ mlp;
d_cv = sd ./ mld;
a_cv = sall ./ mall;

%% ANOMALIES -------------------------------------------------

%% remove linear trend
mf = NaN*ones(size(mf_anu));
lp = NaN*ones(size(lp_anu));
ld = NaN*ones(size(ld_anu));
all = NaN*ones(size(xA));

for i = 1:nid
    %STAGES 
    xi = mf_anu(i,:);
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
    
    xi = lp_anu(i,:);
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
    
    xi = ld_anu(i,:);
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

    xi = xA(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    all(i,:) = dR;
    clear R T t b m tH dR data
    
end

%% anomalies 
amf = mf - mean(mf,2,'omitnan');
alp = lp - mean(lp,2,'omitnan');
ald = ld - mean(ld,2,'omitnan');
aa = all - mean(all,2,'omitnan');

%% var of anomalies by grid cell
vmf = var(amf,0,2,'omitnan');
vlp = var(alp,0,2,'omitnan');
vld = var(ald,0,2,'omitnan');
va  = var(aa,0,2,'omitnan');

%% Anom on grid
Zmf=NaN*ones(ni*nj,nt);
Zlp=NaN*ones(ni*nj,nt);
Zld=NaN*ones(ni*nj,nt);
ZA=NaN*ones(ni*nj,nt);

Zmf(GRD.ID,:)=amf;
Zlp(GRD.ID,:)=alp;
Zld(GRD.ID,:)=ald;
ZA(GRD.ID,:)=aa;

Zmf=reshape(Zmf,ni,nj,nt);
Zlp=reshape(Zlp,ni,nj,nt);
Zld=reshape(Zld,ni,nj,nt);
ZA=reshape(ZA,ni,nj,nt);

%% map info
latlim=[-90 90];
lonlim=[-280 80];
load coastlines

cmYR=cbrewer('seq','YlOrRd',66,'pchip');
cmYOB=cbrewer('seq','YlOrBr',66,'pchip');
cmOR=cbrewer('seq','OrRd',66,'pchip');

%% reshape CV
cvf = NaN*ones(ni,nj);
cvp = NaN*ones(ni,nj);
cvd = NaN*ones(ni,nj);
cva = NaN*ones(ni,nj);

cvf(GRD.ID)=f_cv;
cvp(GRD.ID)=p_cv;
cvd(GRD.ID)=d_cv;
cva(GRD.ID)=a_cv;

%% save
units = 'per year';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
save([fpath 'FEISTY_FOSI_',mod,'nu_interann_var.mat'],...
    'cvf','cvp','cvd','cva',...
    'mmf','mlp','mld','mall',...
    'sf','sp','sd','sall','units');

%%
save([fpath 'FEISTY_FOSI_',mod,'nu_ann_mean_anoms.mat'],'units',...
    'amf','alp','ald','aa',...
    'vmf','vlp','vld','va',...
    'Zmf','Zlp','Zld','ZA');
