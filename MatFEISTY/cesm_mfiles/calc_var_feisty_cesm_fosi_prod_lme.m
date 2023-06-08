% CESM FEISTY FOSI runs
% calc interann variability of prod by lme
% units from per day to per year

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
load([cpath 'LME-mask-POP_gx1v6.mat']);

AREA_OCN = TAREA * 1e-4; %grid cell area in m2
area_vec = AREA_OCN(GRD.ID);

[ni,nj]=size(TLONG);
[nid,nyr]=size(sf_aprod);

%% Per day to per year
sf_aprod = 365 * sf_aprod;
sp_aprod = 365 * sp_aprod;
sd_aprod = 365 * sd_aprod;
mf_aprod = 365 * mf_aprod;
mp_aprod = 365 * mp_aprod;
md_aprod = 365 * md_aprod;
lp_aprod = 365 * lp_aprod;
ld_aprod = 365 * ld_aprod;

%% Groups
xF = mf_aprod;
xP = lp_aprod;
xD = ld_aprod;
xS = sf_aprod + sp_aprod + sd_aprod;
xM = mf_aprod + mp_aprod + md_aprod;
xL = lp_aprod + ld_aprod;
xall = xF + xP + xD;

%% mean & std by lme
tlme = double(lme_mask);
tlme(tlme<0) = nan;
olme = tlme(GRD.ID);

% First create area-weighted mean time series for each LME
lme_sf_mean_ts = NaN*ones(66,nyr);
lme_sp_mean_ts = NaN*ones(66,nyr);
lme_sd_mean_ts = NaN*ones(66,nyr);
lme_mf_mean_ts = NaN*ones(66,nyr);
lme_mp_mean_ts = NaN*ones(66,nyr);
lme_md_mean_ts = NaN*ones(66,nyr);
lme_lp_mean_ts = NaN*ones(66,nyr);
lme_ld_mean_ts = NaN*ones(66,nyr);
lme_a_mean_ts = NaN*ones(66,nyr);
lme_s_mean_ts = NaN*ones(66,nyr);
lme_m_mean_ts = NaN*ones(66,nyr);
lme_l_mean_ts = NaN*ones(66,nyr);
lme_f_mean_ts = NaN*ones(66,nyr);
lme_p_mean_ts = NaN*ones(66,nyr);
lme_d_mean_ts = NaN*ones(66,nyr);

for L=1:66
    lid = find(olme==L);
    
    lme_sf_mean_ts(L,:) = (sum(sf_aprod(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_sp_mean_ts(L,:) = (sum(sp_aprod(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_sd_mean_ts(L,:) = (sum(sd_aprod(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_mf_mean_ts(L,:) = (sum(mf_aprod(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_mp_mean_ts(L,:) = (sum(mp_aprod(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_md_mean_ts(L,:) = (sum(md_aprod(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_lp_mean_ts(L,:) = (sum(lp_aprod(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_ld_mean_ts(L,:) = (sum(ld_aprod(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_a_mean_ts(L,:)  = (sum(xall(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_s_mean_ts(L,:)  = (sum(xS(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_m_mean_ts(L,:)  = (sum(xM(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_l_mean_ts(L,:)  = (sum(xL(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_f_mean_ts(L,:)  = (sum(xF(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_p_mean_ts(L,:)  = (sum(xP(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_d_mean_ts(L,:)  = (sum(xD(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');


end

%% Then calc mean & std dev
% Mean
lme_sf_mean = mean(lme_sf_mean_ts,2,'omitnan');
lme_sp_mean = mean(lme_sp_mean_ts,2,'omitnan');
lme_sd_mean = mean(lme_sd_mean_ts,2,'omitnan');
lme_mf_mean = mean(lme_mf_mean_ts,2,'omitnan');
lme_mp_mean = mean(lme_mp_mean_ts,2,'omitnan');
lme_md_mean = mean(lme_md_mean_ts,2,'omitnan');
lme_lp_mean = mean(lme_lp_mean_ts,2,'omitnan');
lme_ld_mean = mean(lme_ld_mean_ts,2,'omitnan');
lme_a_mean  = mean(lme_a_mean_ts,2,'omitnan');
lme_s_mean  = mean(lme_s_mean_ts,2,'omitnan');
lme_m_mean  = mean(lme_m_mean_ts,2,'omitnan');
lme_l_mean  = mean(lme_l_mean_ts,2,'omitnan');
lme_f_mean  = mean(lme_f_mean_ts,2,'omitnan');
lme_p_mean  = mean(lme_p_mean_ts,2,'omitnan');
lme_d_mean  = mean(lme_d_mean_ts,2,'omitnan');

% Std dev
lme_sf_std = std(lme_sf_mean_ts,0,2,'omitnan');
lme_sp_std = std(lme_sp_mean_ts,0,2,'omitnan');
lme_sd_std = std(lme_sd_mean_ts,0,2,'omitnan');
lme_mf_std = std(lme_mf_mean_ts,0,2,'omitnan');
lme_mp_std = std(lme_mp_mean_ts,0,2,'omitnan');
lme_md_std = std(lme_md_mean_ts,0,2,'omitnan');
lme_lp_std = std(lme_lp_mean_ts,0,2,'omitnan');
lme_ld_std = std(lme_ld_mean_ts,0,2,'omitnan');
lme_a_std  = std(lme_a_mean_ts,0,2,'omitnan');
lme_s_std  = std(lme_s_mean_ts,0,2,'omitnan');
lme_m_std  = std(lme_m_mean_ts,0,2,'omitnan');
lme_l_std  = std(lme_l_mean_ts,0,2,'omitnan');
lme_f_std  = std(lme_f_mean_ts,0,2,'omitnan');
lme_p_std  = std(lme_p_mean_ts,0,2,'omitnan');
lme_d_std  = std(lme_d_mean_ts,0,2,'omitnan');

% Coefficient of variance
lme_sf_cv = lme_sf_std ./ lme_sf_mean;
lme_sp_cv = lme_sp_std ./ lme_sp_mean;
lme_sd_cv = lme_sd_std ./ lme_sd_mean;
lme_mf_cv = lme_mf_std ./ lme_mf_mean;
lme_mp_cv = lme_mp_std ./ lme_mp_mean;
lme_md_cv = lme_md_std ./ lme_md_mean;
lme_lp_cv = lme_lp_std ./ lme_lp_mean;
lme_ld_cv = lme_ld_std ./ lme_ld_mean;
lme_a_cv = lme_a_std ./ lme_a_mean;
lme_s_cv = lme_s_std ./ lme_s_mean;
lme_m_cv = lme_m_std ./ lme_m_mean;
lme_l_cv = lme_l_std ./ lme_l_mean;
lme_f_cv = lme_f_std ./ lme_f_mean;
lme_p_cv = lme_p_std ./ lme_p_mean;
lme_d_cv = lme_d_std ./ lme_d_mean;

%% ANOMALIES -------------------------------------------------

%% remove linear trend
sf = NaN*ones(66,nyr);
sp = NaN*ones(66,nyr);
sd = NaN*ones(66,nyr);
mf = NaN*ones(66,nyr);
mp = NaN*ones(66,nyr);
md = NaN*ones(66,nyr);
lp = NaN*ones(66,nyr);
ld = NaN*ones(66,nyr);

B = NaN*ones(66,nyr);
F = NaN*ones(66,nyr);
P = NaN*ones(66,nyr);
D = NaN*ones(66,nyr);
S = NaN*ones(66,nyr);
M = NaN*ones(66,nyr);
L = NaN*ones(66,nyr);
All = NaN*ones(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62
    %STAGES
    xi = lme_sf_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        sf(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_sp_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        sp(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_sd_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        sd(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_mf_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        mf(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_mp_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        mp(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_md_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        md(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_lp_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        lp(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_ld_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        ld(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    %TYPES
    xi = lme_f_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        F(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_p_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        P(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_d_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        D(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_s_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        S(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_m_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        M(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_l_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        L(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_a_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        All(i,:) = dR;
    end
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
af = F - mean(F,2,'omitnan');
ap = P - mean(P,2,'omitnan');
ad = D - mean(D,2,'omitnan');

%% var of anomalies by lme
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
vf = var(af,0,2,'omitnan');
vp = var(ap,0,2,'omitnan');
vd = var(ad,0,2,'omitnan');

%% save
units = 'per year';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
save([fpath 'FEISTY_FOSI_',mod,'lme_prod_interann_var.mat'],...
    'lme_sf_std','lme_sp_std','lme_sd_std',...
    'lme_mf_std','lme_mp_std','lme_md_std',...
    'lme_lp_std','lme_ld_std','lme_a_std',...
    'lme_s_std','lme_m_std','lme_l_std',...
    'lme_f_std','lme_p_std','lme_d_std',...
    'lme_sf_mean','lme_sp_mean','lme_sd_mean',...
    'lme_mf_mean','lme_mp_mean','lme_md_mean',...
    'lme_lp_mean','lme_ld_mean','lme_a_mean',...
    'lme_s_mean','lme_m_mean','lme_l_mean',...
    'lme_f_mean','lme_p_mean','lme_d_mean',...
    'lme_sf_mean_ts','lme_sp_mean_ts','lme_sd_mean_ts',...
    'lme_mf_mean_ts','lme_mp_mean_ts','lme_md_mean_ts',...
    'lme_lp_mean_ts','lme_ld_mean_ts','lme_a_mean_ts',...
    'lme_s_mean_ts','lme_m_mean_ts','lme_l_mean_ts',...
    'lme_f_mean_ts','lme_p_mean_ts','lme_d_mean_ts',...
    'lme_sf_cv','lme_sp_cv','lme_sd_cv',...
    'lme_mf_cv','lme_mp_cv','lme_md_cv',...
    'lme_lp_cv','lme_ld_cv','lme_a_cv',...
    'lme_s_cv','lme_m_cv','lme_l_cv',...
    'lme_f_cv','lme_p_cv','lme_d_cv','units');

%%
save([fpath 'FEISTY_FOSI_',mod,'lme_prod_ann_mean_anoms.mat'],...
    'asf','asp','asd','amf','amp','amd','alp','ald','aa','as','am','al',...
    'af','ap','ad',...
    'vsf','vsp','vsd','vmf','vmp','vmd','vlp','vld','va','vs','vm','vl',...
    'vf','vp','vd','units');
