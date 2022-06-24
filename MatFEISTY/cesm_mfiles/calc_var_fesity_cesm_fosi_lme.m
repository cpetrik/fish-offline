% CESM FEISTY FOSI runs
% calc interann variability by lme

clear all
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v15_All_fish03_';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
%fpath=['/Volumes/petrik-lab/NC/CESM_MAPP/' cfile '/'];
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
[nid,nyr]=size(sf_abio);

%% Groups
xF = sf_abio + mf_abio;
xP = sp_abio + mp_abio + lp_abio;
xD = sd_abio + md_abio + ld_abio;
xS = sf_abio + sp_abio + sd_abio;
xM = mf_abio + mp_abio + md_abio;
xL = lp_abio + ld_abio;
xB = b_abio;
xall = xF + xP + xD;

%% mean & std by lme
%cpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'LME-mask-POP_gx1v6.mat']);

tlme = double(lme_mask);
tlme(tlme<0) = nan;
olme = tlme(GRD.ID);

lme_sf_mean1 = NaN*ones(66,1);
lme_sp_mean1 = NaN*ones(66,1);
lme_sd_mean1 = NaN*ones(66,1);
lme_mf_mean1 = NaN*ones(66,1);
lme_mp_mean1 = NaN*ones(66,1);
lme_md_mean1 = NaN*ones(66,1);
lme_lp_mean1 = NaN*ones(66,1);
lme_ld_mean1 = NaN*ones(66,1);
lme_b_mean1 = NaN*ones(66,1);
lme_a_mean1 = NaN*ones(66,1);
lme_s_mean1 = NaN*ones(66,1);
lme_m_mean1 = NaN*ones(66,1);
lme_l_mean1 = NaN*ones(66,1);
lme_f_mean1 = NaN*ones(66,1);
lme_p_mean1 = NaN*ones(66,1);
lme_d_mean1 = NaN*ones(66,1);

lme_sf_mean2 = NaN*ones(66,nyr);
lme_sp_mean2 = NaN*ones(66,nyr);
lme_sd_mean2 = NaN*ones(66,nyr);
lme_mf_mean2 = NaN*ones(66,nyr);
lme_mp_mean2 = NaN*ones(66,nyr);
lme_md_mean2 = NaN*ones(66,nyr);
lme_lp_mean2 = NaN*ones(66,nyr);
lme_ld_mean2 = NaN*ones(66,nyr);
lme_b_mean2 = NaN*ones(66,nyr);
lme_a_mean2 = NaN*ones(66,nyr);
lme_s_mean2 = NaN*ones(66,nyr);
lme_m_mean2 = NaN*ones(66,nyr);
lme_l_mean2 = NaN*ones(66,nyr);
lme_f_mean2 = NaN*ones(66,nyr);
lme_p_mean2 = NaN*ones(66,nyr);
lme_d_mean2 = NaN*ones(66,nyr);

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
    
    lme_sf_mean1(L,1) = nanmean(nanmean(sf_abio(lid,:),2));
    lme_sp_mean1(L,1) = nanmean(nanmean(sp_abio(lid,:),2));
    lme_sd_mean1(L,1) = nanmean(nanmean(sd_abio(lid,:),2));
    lme_mf_mean1(L,1) = nanmean(nanmean(mf_abio(lid,:),2));
    lme_mp_mean1(L,1) = nanmean(nanmean(mp_abio(lid,:),2));
    lme_md_mean1(L,1) = nanmean(nanmean(md_abio(lid,:),2));
    lme_lp_mean1(L,1) = nanmean(nanmean(lp_abio(lid,:),2));
    lme_ld_mean1(L,1) = nanmean(nanmean(ld_abio(lid,:),2));
    lme_b_mean1(L,1) = nanmean(nanmean(xB(lid,:),2));
    lme_a_mean1(L,1) = nanmean(nanmean(xall(lid,:),2));
    lme_s_mean1(L,1) = nanmean(nanmean(xS(lid,:),2));
    lme_m_mean1(L,1) = nanmean(nanmean(xM(lid,:),2));
    lme_l_mean1(L,1) = nanmean(nanmean(xL(lid,:),2));
    lme_f_mean1(L,1) = nanmean(nanmean(xF(lid,:),2));
    lme_p_mean1(L,1) = nanmean(nanmean(xP(lid,:),2));
    lme_d_mean1(L,1) = nanmean(nanmean(xD(lid,:),2));
    
    lme_sf_mean2(L,:) = (nanmean(sf_abio(lid,:),1));
    lme_sp_mean2(L,:) = (nanmean(sp_abio(lid,:),1));
    lme_sd_mean2(L,:) = (nanmean(sd_abio(lid,:),1));
    lme_mf_mean2(L,:) = (nanmean(mf_abio(lid,:),1));
    lme_mp_mean2(L,:) = (nanmean(mp_abio(lid,:),1));
    lme_md_mean2(L,:) = (nanmean(md_abio(lid,:),1));
    lme_lp_mean2(L,:) = (nanmean(lp_abio(lid,:),1));
    lme_ld_mean2(L,:) = (nanmean(ld_abio(lid,:),1));
    lme_b_mean2(L,:) = (nanmean(xB(lid,:),1));
    lme_a_mean2(L,:) = (nanmean(xall(lid,:),1));
    lme_s_mean2(L,:) = (nanmean(xS(lid,:),1));
    lme_m_mean2(L,:) = (nanmean(xM(lid,:),1));
    lme_l_mean2(L,:) = (nanmean(xL(lid,:),1));
    lme_f_mean2(L,:) = (nanmean(xF(lid,:),1));
    lme_p_mean2(L,:) = (nanmean(xP(lid,:),1));
    lme_d_mean2(L,:) = (nanmean(xD(lid,:),1));
    
end

%% Coefficient of variance
lme_sf_cv = lme_sf_std ./ lme_sf_mean1;
lme_sp_cv = lme_sp_std ./ lme_sp_mean1;
lme_sd_cv = lme_sd_std ./ lme_sd_mean1;
lme_mf_cv = lme_mf_std ./ lme_mf_mean1;
lme_mp_cv = lme_mp_std ./ lme_mp_mean1;
lme_md_cv = lme_md_std ./ lme_md_mean1;
lme_lp_cv = lme_lp_std ./ lme_lp_mean1;
lme_ld_cv = lme_ld_std ./ lme_ld_mean1;
lme_b_cv = lme_b_std ./ lme_b_mean1;
lme_a_cv = lme_a_std ./ lme_a_mean1;
lme_s_cv = lme_s_std ./ lme_s_mean1;
lme_m_cv = lme_m_std ./ lme_m_mean1;
lme_l_cv = lme_l_std ./ lme_l_mean1;
lme_f_cv = lme_f_std ./ lme_f_mean1;
lme_p_cv = lme_p_std ./ lme_p_mean1;
lme_d_cv = lme_d_std ./ lme_d_mean1;

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
    xi = lme_sf_mean2(i,:);
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
    
    xi = lme_sp_mean2(i,:);
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
    
    xi = lme_sd_mean2(i,:);
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
    
    xi = lme_mf_mean2(i,:);
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
    
    xi = lme_mp_mean2(i,:);
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
    
    xi = lme_md_mean2(i,:);
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
    
    xi = lme_lp_mean2(i,:);
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
    
    xi = lme_ld_mean2(i,:);
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
    xi = lme_b_mean2(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        B(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_f_mean2(i,:);
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
    
    xi = lme_p_mean2(i,:);
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
    
    xi = lme_d_mean2(i,:);
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
    
    xi = lme_s_mean2(i,:);
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
    
    xi = lme_m_mean2(i,:);
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
    
    xi = lme_l_mean2(i,:);
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
    
    xi = lme_a_mean2(i,:);
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

%% save
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
%fpath=['/Volumes/petrik-lab/NC/CESM_MAPP/' cfile '/'];
save([fpath 'FEISTY_FOSI_',mod,'lme_interann_var.mat'],...
    'lme_sf_std','lme_sp_std','lme_sd_std',...
    'lme_mf_std','lme_mp_std','lme_md_std',...
    'lme_lp_std','lme_ld_std','lme_b_std','lme_a_std',...
    'lme_s_std','lme_m_std','lme_l_std',...
    'lme_f_std','lme_p_std','lme_d_std',...
    'lme_sf_mean1','lme_sp_mean1','lme_sd_mean1',...
    'lme_mf_mean1','lme_mp_mean1','lme_md_mean1',...
    'lme_lp_mean1','lme_ld_mean1','lme_b_mean1','lme_a_mean1',...
    'lme_s_mean1','lme_m_mean1','lme_l_mean1',...
    'lme_f_mean1','lme_p_mean1','lme_d_mean1',...
    'lme_sf_mean2','lme_sp_mean2','lme_sd_mean2',...
    'lme_mf_mean2','lme_mp_mean2','lme_md_mean2',...
    'lme_lp_mean2','lme_ld_mean2','lme_b_mean2','lme_a_mean2',...
    'lme_s_mean2','lme_m_mean2','lme_l_mean2',...
    'lme_f_mean2','lme_p_mean2','lme_d_mean2',...
    'lme_sf_cv','lme_sp_cv','lme_sd_cv',...
    'lme_mf_cv','lme_mp_cv','lme_md_cv',...
    'lme_lp_cv','lme_ld_cv','lme_b_cv','lme_a_cv',...
    'lme_s_cv','lme_m_cv','lme_l_cv',...
    'lme_f_cv','lme_p_cv','lme_d_cv');

%%
save([fpath 'FEISTY_FOSI_',mod,'lme_ann_mean_anoms.mat'],...
    'asf','asp','asd','amf','amp','amd','alp','ald','ab','aa','as','am','al',...
    'af','ap','ad',...
    'vsf','vsp','vsd','vmf','vmp','vmd','vlp','vld','vb','va','vs','vm','vl',...
    'vf','vp','vd');
