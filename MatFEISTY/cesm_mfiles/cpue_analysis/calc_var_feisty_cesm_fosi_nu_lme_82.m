% CESM FEISTY FOSI runs
% calc interann variability of nu by lme
% units from per day to per year
% Calculate anomaly time series for different ranges
% 1982-2010 &
% 1982-2015

clear 
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v15_All_fish03_'; % v15_All_fish03; 'v15_obsfish_'

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
lme_f_awmean_ts = NaN*ones(66,nyr);
lme_p_awmean_ts = NaN*ones(66,nyr);
lme_d_awmean_ts = NaN*ones(66,nyr);
lme_a_awmean_ts = NaN*ones(66,nyr);

for L=1:66
    lid = find(olme==L);
    
    lme_f_awmean_ts(L,:)  = (sum(xF(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_p_awmean_ts(L,:)  = (sum(xP(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_d_awmean_ts(L,:)  = (sum(xD(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_a_awmean_ts(L,:)  = (sum(xA(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');

end

%% Next select specific yrs
fyr = 1948:2015;
eyr = 1982:2010; % subset effort years
yyr = 1982:2015; % subset catch years

lid=1:66;

[~,fide] = intersect(fyr,eyr);
[~,fidc] = intersect(fyr,yyr);

%% Then calc mean & std dev
% Mean
lme_f10_mean  = mean(lme_f_awmean_ts(:,fide),2,'omitnan');
lme_p10_mean  = mean(lme_p_awmean_ts(:,fide),2,'omitnan');
lme_d10_mean  = mean(lme_d_awmean_ts(:,fide),2,'omitnan');
lme_a10_mean  = mean(lme_a_awmean_ts(:,fide),2,'omitnan');

% Std dev
lme_f10_std  = std(lme_f_awmean_ts(:,fide),0,2,'omitnan');
lme_p10_std  = std(lme_p_awmean_ts(:,fide),0,2,'omitnan');
lme_d10_std  = std(lme_d_awmean_ts(:,fide),0,2,'omitnan');
lme_a10_std  = std(lme_a_awmean_ts(:,fide),0,2,'omitnan');

% Coefficient of variance
lme_f10_cv = lme_f10_std ./ lme_f10_mean;
lme_p10_cv = lme_p10_std ./ lme_p10_mean;
lme_d10_cv = lme_d10_std ./ lme_d10_mean;
lme_a10_cv = lme_a10_std ./ lme_a10_mean;


% Mean
lme_f15_mean  = mean(lme_f_awmean_ts(:,fidc),2,'omitnan');
lme_p15_mean  = mean(lme_p_awmean_ts(:,fidc),2,'omitnan');
lme_d15_mean  = mean(lme_d_awmean_ts(:,fidc),2,'omitnan');
lme_a15_mean  = mean(lme_a_awmean_ts(:,fidc),2,'omitnan');

% Std dev
lme_f15_std  = std(lme_f_awmean_ts(:,fidc),0,2,'omitnan');
lme_p15_std  = std(lme_p_awmean_ts(:,fidc),0,2,'omitnan');
lme_d15_std  = std(lme_d_awmean_ts(:,fidc),0,2,'omitnan');
lme_a15_std  = std(lme_a_awmean_ts(:,fidc),0,2,'omitnan');

% Coefficient of variance
lme_f15_cv = lme_f15_std ./ lme_f15_mean;
lme_p15_cv = lme_p15_std ./ lme_p15_mean;
lme_d15_cv = lme_d15_std ./ lme_d15_mean;
lme_a15_cv = lme_a15_std ./ lme_a15_mean;

%% ANOMALIES -------------------------------------------------

nte = length(eyr); %2010
nty = length(yyr); %2015

%% remove linear trend
F10 = NaN*ones(66,nte);
P10 = NaN*ones(66,nte);
D10 = NaN*ones(66,nte);
A10 = NaN*ones(66,nte);
F15 = NaN*ones(66,nty);
P15 = NaN*ones(66,nty);
D15 = NaN*ones(66,nty);
A15 = NaN*ones(66,nty);

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

    %TYPES - Catch
    xi = lme_f_awmean_ts(i,fidc);
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
    
    xi = lme_p_awmean_ts(i,fidc);
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
    
    xi = lme_d_awmean_ts(i,fidc);
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

    xi = lme_a_awmean_ts(i,fidc);
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
anf10 = F10 - mean(F10,2,'omitnan');
anp10 = P10 - mean(P10,2,'omitnan');
and10 = D10 - mean(D10,2,'omitnan');
ana10 = A10 - mean(A10,2,'omitnan');

anf15 = F15 - mean(F15,2,'omitnan');
anp15 = P15 - mean(P15,2,'omitnan');
and15 = D15 - mean(D15,2,'omitnan');
ana15 = A15 - mean(A15,2,'omitnan');

%% var of anomalies by lme
vnf10 = var(anf10,0,2,'omitnan');
vnp10 = var(anp10,0,2,'omitnan');
vnd10 = var(and10,0,2,'omitnan');
vna10 = var(ana10,0,2,'omitnan');

vnf15 = var(anf15,0,2,'omitnan');
vnp15 = var(anp15,0,2,'omitnan');
vnd15 = var(and15,0,2,'omitnan');
vna15 = var(ana15,0,2,'omitnan');

%% save
units = 'per year';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
save([fpath 'FEISTY_FOSI_',mod,'lme_nu_interann_var_1982_2010_2015.mat'],...
    'lme_f_awmean_ts','lme_p_awmean_ts','lme_d_awmean_ts','lme_a_awmean_ts',...
    'lme_f10_std','lme_p10_std','lme_d10_std','lme_a10_std',...
    'lme_f10_mean','lme_p10_mean','lme_d10_mean','lme_a10_mean',...
    'lme_f10_cv','lme_p10_cv','lme_d10_cv','lme_a10_cv','eyr','yyr',...
    'lme_f15_std','lme_p15_std','lme_d15_std','lme_a15_std',...
    'lme_f15_mean','lme_p15_mean','lme_d15_mean','lme_a15_mean',...
    'lme_f15_cv','lme_p15_cv','lme_d15_cv','lme_a15_cv','units');

%%
save([fpath 'FEISTY_FOSI_',mod,'lme_nu_ann_mean_anoms_1982_2010_2015.mat'],...
    'anf10','anp10','and10','ana10',...
    'vnf10','vnp10','vnd10','vna10',...
    'anf15','anp15','and15','ana15',...
    'vnf15','vnp15','vnd15','vna15','units','eyr','yyr');
