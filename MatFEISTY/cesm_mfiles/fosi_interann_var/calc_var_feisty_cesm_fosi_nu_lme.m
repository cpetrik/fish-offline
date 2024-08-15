% CESM FEISTY FOSI runs
% calc interann variability of nu by lme
% units from per day to per year

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
load([fpath 'Annual_Means_FOSI_' mod cfile '.mat'],...
    'mf_anu','lp_anu','ld_anu');

% Map data
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

AREA_OCN = TAREA * 1e-4; %grid cell area in m2
area_vec = AREA_OCN(GRD.ID);

[ni,nj]=size(TLONG);
[nid,nyr]=size(mf_anu);

%% Per day to per year
mf_anu = 365 * mf_anu;
lp_anu = 365 * lp_anu;
ld_anu = 365 * ld_anu;

%% Groups
xF = mf_anu;
xP = lp_anu;
xD = ld_anu;

xA = (mf_anu+lp_anu+ld_anu);

%% mean & std by lme
tlme = double(lme_mask);
tlme(tlme<0) = nan;
olme = tlme(GRD.ID);

% First create area-weighted mean time series for each LME
lme_f_mean_ts = NaN*ones(66,nyr);
lme_p_mean_ts = NaN*ones(66,nyr);
lme_d_mean_ts = NaN*ones(66,nyr);
lme_a_mean_ts = NaN*ones(66,nyr);

for L=1:66
    lid = find(olme==L);
    
    lme_f_mean_ts(L,:)  = (sum(xF(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_p_mean_ts(L,:)  = (sum(xP(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_d_mean_ts(L,:)  = (sum(xD(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_a_mean_ts(L,:)  = (sum(xA(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');

end

%% Then calc mean & std dev
% Mean
lme_f_mean  = mean(lme_f_mean_ts,2,'omitnan');
lme_p_mean  = mean(lme_p_mean_ts,2,'omitnan');
lme_d_mean  = mean(lme_d_mean_ts,2,'omitnan');
lme_a_mean  = mean(lme_a_mean_ts,2,'omitnan');

% Std dev
lme_f_std  = std(lme_f_mean_ts,0,2,'omitnan');
lme_p_std  = std(lme_p_mean_ts,0,2,'omitnan');
lme_d_std  = std(lme_d_mean_ts,0,2,'omitnan');
lme_a_std  = std(lme_a_mean_ts,0,2,'omitnan');

% Coefficient of variance
lme_f_cv = lme_f_std ./ lme_f_mean;
lme_p_cv = lme_p_std ./ lme_p_mean;
lme_d_cv = lme_d_std ./ lme_d_mean;
lme_a_cv = lme_a_std ./ lme_a_mean;

%% ANOMALIES -------------------------------------------------

%% remove linear trend
F = NaN*ones(66,nyr);
P = NaN*ones(66,nyr);
D = NaN*ones(66,nyr);
A = NaN*ones(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62
    
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
        A(i,:) = dR;
    end
    clear R T t b m tH dR data
   
end

%% anomalies
af = F - mean(F,2,'omitnan');
ap = P - mean(P,2,'omitnan');
ad = D - mean(D,2,'omitnan');
aa = A - mean(A,2,'omitnan');

%% var of anomalies by lme
vf = var(af,0,2,'omitnan');
vp = var(ap,0,2,'omitnan');
vd = var(ad,0,2,'omitnan');
va = var(aa,0,2,'omitnan');

%% save
units = 'per year';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
save([fpath 'FEISTY_FOSI_',mod,'lme_nu_interann_var.mat'],...
    'lme_f_std','lme_p_std','lme_d_std','lme_a_std',...
    'lme_f_mean','lme_p_mean','lme_d_mean','lme_a_mean',...
    'lme_f_mean_ts','lme_p_mean_ts','lme_d_mean_ts','lme_a_mean_ts',...
    'lme_f_cv','lme_p_cv','lme_d_cv','lme_a_cv','units');

%%
save([fpath 'FEISTY_FOSI_',mod,'lme_nu_ann_mean_anoms.mat'],...
    'af','ap','ad','aa',...
    'vf','vp','vd','va','units');
