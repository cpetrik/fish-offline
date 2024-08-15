% CESM FEISTY FOSI runs
% calc interann variability by lme
% Calculate anomaly time series for different ranges
% 1982-2010 &
% 1982-2015

clear
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v15_obsfish_'; % v15_All_fish03; 'v15_obsfish_'

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
xall = xF + xP + xD;

%% area-weighted means per LME
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
fyr = 1948:2010;
eyr = 1982:2010; % subset effort years

lid=1:66;

[~,fide] = intersect(fyr,eyr);

%% Then calc mean & std dev
% Mean
lme_a10_mean  = mean(lme_a_awmean_ts(:,fide),2,'omitnan');
lme_f10_mean  = mean(lme_f_awmean_ts(:,fide),2,'omitnan');
lme_p10_mean  = mean(lme_p_awmean_ts(:,fide),2,'omitnan');
lme_d10_mean  = mean(lme_d_awmean_ts(:,fide),2,'omitnan');

% Std dev
lme_a10_std  = std(lme_a_awmean_ts(:,fide),0,2,'omitnan');
lme_f10_std  = std(lme_f_awmean_ts(:,fide),0,2,'omitnan');
lme_p10_std  = std(lme_p_awmean_ts(:,fide),0,2,'omitnan');
lme_d10_std  = std(lme_d_awmean_ts(:,fide),0,2,'omitnan');

% Coefficient of variance
lme_a10_cv = lme_a10_std ./ lme_a10_mean;
lme_f10_cv = lme_f10_std ./ lme_f10_mean;
lme_p10_cv = lme_p10_std ./ lme_p10_mean;
lme_d10_cv = lme_d10_std ./ lme_d10_mean;

%% ANOMALIES -------------------------------------------------

nte = length(eyr); %2010

%% remove linear trend
F10 = NaN*ones(66,nte);
P10 = NaN*ones(66,nte);
D10 = NaN*ones(66,nte);
A10 = NaN*ones(66,nte);

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
        F10(i,:) = dR;
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
        P10(i,:) = dR;
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
        D10(i,:) = dR;
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
        A10(i,:) = dR;
    end
    clear R T t b m tH dR data

end

%% anomalies
aya10 = A10 - mean(A10,2,'omitnan');
ayf10 = F10 - mean(F10,2,'omitnan');
ayp10 = P10 - mean(P10,2,'omitnan');
ayd10 = D10 - mean(D10,2,'omitnan');

%% var of anomalies by grid cell
vya10 = var(aya10,0,2,'omitnan');
vyf10 = var(ayf10,0,2,'omitnan');
vyp10 = var(ayp10,0,2,'omitnan');
vyd10 = var(ayd10,0,2,'omitnan');

%% save
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
save([fpath 'FEISTY_FOSI_',mod,'lme_catch_interann_var_1982_2010.mat'],...
    'lme_f_awmean_ts','lme_p_awmean_ts','lme_d_awmean_ts','lme_a_awmean_ts',...
    'lme_f10_std','lme_p10_std','lme_d10_std','lme_a10_std',...
    'lme_f10_mean','lme_p10_mean','lme_d10_mean','lme_a10_mean',...
    'lme_f10_cv','lme_p10_cv','lme_d10_cv','lme_a10_cv','eyr');

%%
save([fpath 'FEISTY_FOSI_',mod,'lme_catch_ann_mean_anoms_1982_2010.mat'],...
    'aya10','ayf10','ayp10','ayd10',...
    'vya10','vyf10','vyp10','vyd10','eyr');
