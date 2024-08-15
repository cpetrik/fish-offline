% CESM FEISTY FOSI runs
% calc interann variability by lme

clear
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v15_obsfish_';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
ppath = [pp cfile '/FOSI/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([fpath 'Annual_Means_FOSI_' mod cfile '.mat']);

% Map data
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

AREA_OCN = TAREA * 1e-4; %grid cell area in m2
area_vec = AREA_OCN(GRD.ID);

[ni,nj]=size(TLONG);
[nid,nyr]=size(mf_tac);

%% Groups
xF = mf_tac;
xP = mp_tac + lp_tac;
xD = md_tac + ld_tac;
xM = mf_tac + mp_tac + md_tac;
xL = lp_tac + ld_tac;
xall = xF + xP + xD;

%% mean & std by lme
tlme = double(lme_mask);
tlme(tlme<0) = nan;
olme = tlme(GRD.ID);

lme_mf_mean_ts = NaN*ones(66,nyr);
lme_mp_mean_ts = NaN*ones(66,nyr);
lme_md_mean_ts = NaN*ones(66,nyr);
lme_lp_mean_ts = NaN*ones(66,nyr);
lme_ld_mean_ts = NaN*ones(66,nyr);
lme_a_mean_ts  = NaN*ones(66,nyr);
lme_m_mean_ts  = NaN*ones(66,nyr);
lme_l_mean_ts  = NaN*ones(66,nyr);
lme_f_mean_ts  = NaN*ones(66,nyr);
lme_p_mean_ts  = NaN*ones(66,nyr);
lme_d_mean_ts  = NaN*ones(66,nyr);

% First create area-weighted mean time series for each LME
for L=1:66
    lid = find(olme==L);

    lme_mf_mean_ts(L,:) = (sum(mf_tac(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_mp_mean_ts(L,:) = (sum(mp_tac(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_md_mean_ts(L,:) = (sum(md_tac(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_lp_mean_ts(L,:) = (sum(lp_tac(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_ld_mean_ts(L,:) = (sum(ld_tac(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_a_mean_ts(L,:)  = (sum(xall(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_m_mean_ts(L,:)  = (sum(xM(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_l_mean_ts(L,:)  = (sum(xL(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_f_mean_ts(L,:)  = (sum(xF(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_p_mean_ts(L,:)  = (sum(xP(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_d_mean_ts(L,:)  = (sum(xD(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');

end

%% Then calc mean & std dev
% Mean
lme_mf_mean = mean(lme_mf_mean_ts,2,'omitnan');
lme_mp_mean = mean(lme_mp_mean_ts,2,'omitnan');
lme_md_mean = mean(lme_md_mean_ts,2,'omitnan');
lme_lp_mean = mean(lme_lp_mean_ts,2,'omitnan');
lme_ld_mean = mean(lme_ld_mean_ts,2,'omitnan');
lme_a_mean  = mean(lme_a_mean_ts,2,'omitnan');
lme_m_mean  = mean(lme_m_mean_ts,2,'omitnan');
lme_l_mean  = mean(lme_l_mean_ts,2,'omitnan');
lme_f_mean  = mean(lme_f_mean_ts,2,'omitnan');
lme_p_mean  = mean(lme_p_mean_ts,2,'omitnan');
lme_d_mean  = mean(lme_d_mean_ts,2,'omitnan');

% Std dev
lme_mf_std = std(lme_mf_mean_ts,0,2,'omitnan');
lme_mp_std = std(lme_mp_mean_ts,0,2,'omitnan');
lme_md_std = std(lme_md_mean_ts,0,2,'omitnan');
lme_lp_std = std(lme_lp_mean_ts,0,2,'omitnan');
lme_ld_std = std(lme_ld_mean_ts,0,2,'omitnan');
lme_a_std  = std(lme_a_mean_ts,0,2,'omitnan');
lme_m_std  = std(lme_m_mean_ts,0,2,'omitnan');
lme_l_std  = std(lme_l_mean_ts,0,2,'omitnan');
lme_f_std  = std(lme_f_mean_ts,0,2,'omitnan');
lme_p_std  = std(lme_p_mean_ts,0,2,'omitnan');
lme_d_std  = std(lme_d_mean_ts,0,2,'omitnan');

% Coefficient of variance
lme_mf_cv = lme_mf_std ./ lme_mf_mean;
lme_mp_cv = lme_mp_std ./ lme_mp_mean;
lme_md_cv = lme_md_std ./ lme_md_mean;
lme_lp_cv = lme_lp_std ./ lme_lp_mean;
lme_ld_cv = lme_ld_std ./ lme_ld_mean;
lme_a_cv = lme_a_std ./ lme_a_mean;
lme_m_cv = lme_m_std ./ lme_m_mean;
lme_l_cv = lme_l_std ./ lme_l_mean;
lme_f_cv = lme_f_std ./ lme_f_mean;
lme_p_cv = lme_p_std ./ lme_p_mean;
lme_d_cv = lme_d_std ./ lme_d_mean;

%% ANOMALIES -------------------------------------------------

%% remove linear trend
mf = NaN*ones(66,nyr);
mp = NaN*ones(66,nyr);
md = NaN*ones(66,nyr);
lp = NaN*ones(66,nyr);
ld = NaN*ones(66,nyr);

F = NaN*ones(66,nyr);
P = NaN*ones(66,nyr);
D = NaN*ones(66,nyr);
M = NaN*ones(66,nyr);
L = NaN*ones(66,nyr);
All = NaN*ones(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62
    %STAGES
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
amf = mf - mean(mf,2,'omitnan');
amp = mp - mean(mp,2,'omitnan');
amd = md - mean(md,2,'omitnan');
alp = lp - mean(lp,2,'omitnan');
ald = ld - mean(ld,2,'omitnan');
aa = All - mean(All,2,'omitnan');
am = M - mean(M,2,'omitnan');
al = L - mean(L,2,'omitnan');
af = F - mean(F,2,'omitnan');
ap = P - mean(P,2,'omitnan');
ad = D - mean(D,2,'omitnan');

%% var of anomalies by grid cell
vmf = var(amf,0,2,'omitnan');
vmp = var(amp,0,2,'omitnan');
vmd = var(amd,0,2,'omitnan');
vlp = var(alp,0,2,'omitnan');
vld = var(ald,0,2,'omitnan');
va = var(aa,0,2,'omitnan');
vm = var(am,0,2,'omitnan');
vl = var(al,0,2,'omitnan');
vf = var(af,0,2,'omitnan');
vp = var(ap,0,2,'omitnan');
vd = var(ad,0,2,'omitnan');

%% save
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
save([fpath 'FEISTY_FOSI_',mod,'lme_catch_interann_var.mat'],...
    'lme_mf_std','lme_mp_std','lme_md_std',...
    'lme_lp_std','lme_ld_std','lme_a_std',...
    'lme_m_std','lme_l_std',...
    'lme_f_std','lme_p_std','lme_d_std',...
    'lme_mf_mean','lme_mp_mean','lme_md_mean',...
    'lme_lp_mean','lme_ld_mean','lme_a_mean',...
    'lme_m_mean','lme_l_mean',...
    'lme_f_mean','lme_p_mean','lme_d_mean',...
    'lme_mf_mean_ts','lme_mp_mean_ts','lme_md_mean_ts',...
    'lme_lp_mean_ts','lme_ld_mean_ts','lme_a_mean_ts',...
    'lme_m_mean_ts','lme_l_mean_ts',...
    'lme_f_mean_ts','lme_p_mean_ts','lme_d_mean_ts',...
    'lme_mf_cv','lme_mp_cv','lme_md_cv',...
    'lme_lp_cv','lme_ld_cv','lme_a_cv',...
    'lme_m_cv','lme_l_cv',...
    'lme_f_cv','lme_p_cv','lme_d_cv');

%%
save([fpath 'FEISTY_FOSI_',mod,'lme_catch_ann_mean_anoms.mat'],...
    'amf','amp','amd','alp','ald','aa','am','al',...
    'af','ap','ad',...
    'vmf','vmp','vmd','vlp','vld','va','vm','vl',...
    'vf','vp','vd');
