% CESM FEISTY FOSI runs
% calc interann variability by lme
% Calculate anomaly time series for different ranges
% 1997-2015 

clear
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v15_obsfish2015_'; % v15_All_fish03; 'v15_obsfish_'

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
[nid,nyr]=size(sf_abio);

%% Groups
xF = sf_abio + mf_abio;
xP = sp_abio + mp_abio + lp_abio;
xD = sd_abio + md_abio + ld_abio;
xS = sf_abio + sp_abio + sd_abio;
xall = xF + xP + xD;

%% area-weighted LME means 
tlme = double(lme_mask);
tlme(tlme<0) = nan;
olme = tlme(GRD.ID);

lme_a_awmean_ts  = NaN*ones(66,nyr);
lme_f_awmean_ts  = NaN*ones(66,nyr);
lme_p_awmean_ts  = NaN*ones(66,nyr);
lme_d_awmean_ts  = NaN*ones(66,nyr);

% First create area-weighted mean time series for each LME
for L=1:66
    lid = find(olme==L);

    lme_a_awmean_ts(L,:)  = (sum(xall(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_f_awmean_ts(L,:)  = (sum(xF(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_p_awmean_ts(L,:)  = (sum(xP(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_d_awmean_ts(L,:)  = (sum(xD(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');

end

%% Next select specific yrs
fyr = 1948:2015;
eyr = 1997:2015; % subset effort years

lid=1:66;

[~,fide] = intersect(fyr,eyr);

%% Then calc mean & std dev
% Mean
lme_a15_mean  = mean(lme_a_awmean_ts(:,fide),2,'omitnan');
lme_f15_mean  = mean(lme_f_awmean_ts(:,fide),2,'omitnan');
lme_p15_mean  = mean(lme_p_awmean_ts(:,fide),2,'omitnan');
lme_d15_mean  = mean(lme_d_awmean_ts(:,fide),2,'omitnan');

% Std dev
lme_a15_std  = std(lme_a_awmean_ts(:,fide),0,2,'omitnan');
lme_f15_std  = std(lme_f_awmean_ts(:,fide),0,2,'omitnan');
lme_p15_std  = std(lme_p_awmean_ts(:,fide),0,2,'omitnan');
lme_d15_std  = std(lme_d_awmean_ts(:,fide),0,2,'omitnan');

% Coefficient of variance
lme_a15_cv = lme_a15_std ./ lme_a15_mean;
lme_f15_cv = lme_f15_std ./ lme_f15_mean;
lme_p15_cv = lme_p15_std ./ lme_p15_mean;
lme_d15_cv = lme_d15_std ./ lme_d15_mean;

%% ANOMALIES -------------------------------------------------

nte = length(eyr); %2015

%% remove linear trend
F15 = NaN*ones(66,nte);
P15 = NaN*ones(66,nte);
D15 = NaN*ones(66,nte);
A15 = NaN*ones(66,nte);

for i = 1:66 %in future should remove interior seas 23, 33, 62
    
    %TYPES - CPUE 
    xi = lme_f_awmean_ts(i,fide);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        F15(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = lme_p_awmean_ts(i,fide);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        P15(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = lme_d_awmean_ts(i,fide);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        D15(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = lme_a_awmean_ts(i,fide);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        A15(i,:) = dR;
    end
    clear R T t b m tH dR data

end

%% anomalies
aba15 = A15 - mean(A15,2,'omitnan');
abf15 = F15 - mean(F15,2,'omitnan');
abp15 = P15 - mean(P15,2,'omitnan');
abd15 = D15 - mean(D15,2,'omitnan');

%% var of anomalies by grid cell
vba15 = var(aba15,0,2,'omitnan');
vbf15 = var(abf15,0,2,'omitnan');
vbp15 = var(abp15,0,2,'omitnan');
vbd15 = var(abd15,0,2,'omitnan');

%% save
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
save([fpath 'FEISTY_FOSI_',mod,'lme_biom_interann_var_1997_2015.mat'],...
    'lme_a_awmean_ts','lme_f_awmean_ts','lme_p_awmean_ts','lme_d_awmean_ts',...
    'lme_a15_mean','lme_f15_mean','lme_p15_mean','lme_d15_mean',...
    'lme_a15_std','lme_f15_std','lme_p15_std','lme_d15_std',...
    'lme_a15_cv','lme_f15_cv','lme_p15_cv','lme_d15_cv','eyr');

%%
save([fpath 'FEISTY_FOSI_',mod,'lme_biom_ann_mean_anoms_1997_2015.mat'],...
    'aba15','abf15','abp15','abd15',...
    'vba15','vbf15','vbp15','vbd15','eyr');
