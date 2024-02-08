% CESM FEISTY FOSI runs
% calc interann variability of gam & rec by lme
% units from per day to per year

clear 
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
sims = {'v15_All_fish03_';'v15_climatol_';'v15_varTemp_';'v15_varFood_'};
mod = sims{1};

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
ppath = [pp cfile '/FOSI/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([fpath 'Annual_Means_FOSI_' mod cfile '.mat'],...
    'mf_agam','lp_agam','ld_agam','mf_arec','lp_arec','ld_arec');

% Map data
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

AREA_OCN = TAREA * 1e-4; %grid cell area in m2
area_vec = AREA_OCN(GRD.ID);

[ni,nj]=size(TLONG);
[nid,nyr]=size(mf_agam);

%% Per day to per year
mf_agam = 365 * mf_agam;
lp_agam = 365 * lp_agam;
ld_agam = 365 * ld_agam;

mf_arec = 365 * mf_arec;
lp_arec = 365 * lp_arec;
ld_arec = 365 * ld_arec;

%% Groups
xF = mf_agam;
xP = lp_agam;
xD = ld_agam;

rF = mf_arec;
rP = lp_arec;
rD = ld_arec;

xA = (mf_agam+lp_agam+ld_agam)./3;
rA = (mf_arec+lp_arec+ld_arec)./3;

%% mean & std by lme
tlme = double(lme_mask);
tlme(tlme<0) = nan;
olme = tlme(GRD.ID);

% First create area-weighted mean time series for each LME
lme_gf_mean_ts = NaN*ones(66,nyr);
lme_gp_mean_ts = NaN*ones(66,nyr);
lme_gd_mean_ts = NaN*ones(66,nyr);
lme_ga_mean_ts = NaN*ones(66,nyr);

lme_rf_mean_ts = NaN*ones(66,nyr);
lme_rp_mean_ts = NaN*ones(66,nyr);
lme_rd_mean_ts = NaN*ones(66,nyr);
lme_ra_mean_ts = NaN*ones(66,nyr);

for L=1:66
    lid = find(olme==L);
    
    lme_gf_mean_ts(L,:)  = (sum(xF(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_gp_mean_ts(L,:)  = (sum(xP(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_gd_mean_ts(L,:)  = (sum(xD(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_ga_mean_ts(L,:)  = (sum(xA(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');

    lme_rf_mean_ts(L,:)  = (sum(rF(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_rp_mean_ts(L,:)  = (sum(rP(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_rd_mean_ts(L,:)  = (sum(rD(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_ra_mean_ts(L,:)  = (sum(rA(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');

end

%% Then calc mean & std dev
% Mean
lme_gf_mean  = mean(lme_gf_mean_ts,2,'omitnan');
lme_gp_mean  = mean(lme_gp_mean_ts,2,'omitnan');
lme_gd_mean  = mean(lme_gd_mean_ts,2,'omitnan');
lme_ga_mean  = mean(lme_ga_mean_ts,2,'omitnan');

% Std dev
lme_gf_std  = std(lme_gf_mean_ts,0,2,'omitnan');
lme_gp_std  = std(lme_gp_mean_ts,0,2,'omitnan');
lme_gd_std  = std(lme_gd_mean_ts,0,2,'omitnan');
lme_ga_std  = std(lme_ga_mean_ts,0,2,'omitnan');

% Coefficient of variance
lme_gf_cv = lme_gf_std ./ lme_gf_mean;
lme_gp_cv = lme_gp_std ./ lme_gp_mean;
lme_gd_cv = lme_gd_std ./ lme_gd_mean;
lme_ga_cv = lme_ga_std ./ lme_ga_mean;

% Mean
lme_rf_mean  = mean(lme_rf_mean_ts,2,'omitnan');
lme_rp_mean  = mean(lme_rp_mean_ts,2,'omitnan');
lme_rd_mean  = mean(lme_rd_mean_ts,2,'omitnan');
lme_ra_mean  = mean(lme_ra_mean_ts,2,'omitnan');

% Std dev
lme_rf_std  = std(lme_rf_mean_ts,0,2,'omitnan');
lme_rp_std  = std(lme_rp_mean_ts,0,2,'omitnan');
lme_rd_std  = std(lme_rd_mean_ts,0,2,'omitnan');
lme_ra_std  = std(lme_ra_mean_ts,0,2,'omitnan');

% Coefficient of variance
lme_rf_cv = lme_rf_std ./ lme_rf_mean;
lme_rp_cv = lme_rp_std ./ lme_rp_mean;
lme_rd_cv = lme_rd_std ./ lme_rd_mean;
lme_ra_cv = lme_ra_std ./ lme_ra_mean;

%% ANOMALIES -------------------------------------------------

%% remove linear trend
gF = NaN*ones(66,nyr);
gP = NaN*ones(66,nyr);
gD = NaN*ones(66,nyr);
gA = NaN*ones(66,nyr);
rF = NaN*ones(66,nyr);
rP = NaN*ones(66,nyr);
rD = NaN*ones(66,nyr);
rA = NaN*ones(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62
    
    %TYPES - gamma
    xi = lme_gf_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        gF(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_gp_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        gP(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_gd_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        gD(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = lme_ga_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        gA(i,:) = dR;
    end
    clear R T t b m tH dR data
   
    %TYPES - rec
    xi = lme_rf_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        rF(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_rp_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        rP(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_rd_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        rD(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = lme_ra_mean_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        rA(i,:) = dR;
    end
    clear R T t b m tH dR data
end

%% anomalies
agf = gF - mean(gF,2,'omitnan');
agp = gP - mean(gP,2,'omitnan');
agd = gD - mean(gD,2,'omitnan');
aga = gA - mean(gA,2,'omitnan');

arf = rF - mean(rF,2,'omitnan');
arp = rP - mean(rP,2,'omitnan');
ard = rD - mean(rD,2,'omitnan');
ara = rA - mean(rA,2,'omitnan');

%% var of anomalies by lme
vgf = var(agf,0,2,'omitnan');
vgp = var(agp,0,2,'omitnan');
vgd = var(agd,0,2,'omitnan');
vga = var(aga,0,2,'omitnan');

vrf = var(arf,0,2,'omitnan');
vrp = var(arp,0,2,'omitnan');
vrd = var(ard,0,2,'omitnan');
vra = var(ara,0,2,'omitnan');

%% save
units = 'per year';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
save([fpath 'FEISTY_FOSI_',mod,'lme_gam_rec_interann_var.mat'],...
    'lme_gf_std','lme_gp_std','lme_gd_std','lme_ga_std',...
    'lme_gf_mean','lme_gp_mean','lme_gd_mean','lme_ga_mean',...
    'lme_gf_mean_ts','lme_gp_mean_ts','lme_gd_mean_ts','lme_ga_mean_ts',...
    'lme_gf_cv','lme_gp_cv','lme_gd_cv','lme_ga_cv',...
    'lme_rf_std','lme_rp_std','lme_rd_std','lme_ra_std',...
    'lme_rf_mean','lme_rp_mean','lme_rd_mean','lme_ra_mean',...
    'lme_rf_mean_ts','lme_rp_mean_ts','lme_rd_mean_ts','lme_ra_mean_ts',...
    'lme_rf_cv','lme_rp_cv','lme_rd_cv','lme_ra_cv','units');

%%
save([fpath 'FEISTY_FOSI_',mod,'lme_gam_rec_ann_mean_anoms.mat'],...
    'agf','agp','agd','aga',...
    'vgf','vgp','vgd','vga',...
    'arf','arp','ard','ara',...
    'vrf','vrp','vrd','vra','units');
