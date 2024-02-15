% Satellite chl data

clear 
close all

cpath = '/Volumes/petrik-lab/Feisty/Obs_data/Chl/';
tpath = '/Volumes/petrik-lab/Feisty/Obs_data/OISST/';

%%
load([cpath 'SEAWIFS_1997_2001_MODIS_2002_2022.YR.CHL.chlor_a_1deg.mat']);
cyr = yrs;
clear yrs

load([tpath 'oisst.annual.mean.1deg_1982_2020.mat']);
tyr = yrs;
clear yrs

%% LME data for regular 1 deg grid
opath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';

load([opath 'lme_gfdl-mom6-cobalt2_onedeg_temporary.mat'])
load([opath 'cellarea_onedeg.mat']) %area in m2
tlme = fliplr(tlme);
cell_area = fliplr(cell_area);

%% LME only cells
[ni,nj,nts] = size(tos);
[~,~,ntc] = size(chl);
sst2 = reshape(tos,ni*nj,nts);
chl2 = reshape(chl,ni*nj,ntc);

oid = find(~isnan(tlme(:)));
olme = tlme(oid);
sst2 = sst2(oid,:);
chl2 = chl2(oid,:);
area_vec = cell_area(oid);

%% First create area-weighted mean time series for each LME
lme_sst_ts = NaN*ones(66,length(tyr));
lme_chl_ts = NaN*ones(66,length(cyr));

for L=1:66
    lid = find(olme==L);
    
    lme_sst_ts(L,:)  = (sum(sst2(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    lme_chl_ts(L,:)  = (sum(chl2(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');

end

%% Then calc mean & std dev
% Mean
lme_sst_mean  = mean(lme_sst_ts,2,'omitnan');
lme_chl_mean  = mean(lme_chl_ts,2,'omitnan');

% Std dev
lme_sst_std  = std(lme_sst_ts,0,2,'omitnan');
lme_chl_std  = std(lme_chl_ts,0,2,'omitnan');

% Coefficient of variance
lme_sst_cv = lme_sst_std ./ lme_sst_mean;
lme_chl_cv = lme_chl_std ./ lme_chl_mean;

%% ANOMALIES -------------------------------------------------

% remove linear trend
tSST = NaN*ones(66,length(tyr));
tCHL = NaN*ones(66,length(cyr));

for i = 1:66 %in future should remove interior seas 23, 33, 62
    
    %SST
    xi = lme_sst_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        tSST(i,:) = dR;
    end
    clear R T t b m tH dR data
   
    %Chl
    xi = lme_chl_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        tCHL(i,:) = dR;
    end
    clear R T t b m tH dR data
   
end

%% anomalies
asst = tSST - mean(tSST,2,'omitnan');
achl = tCHL - mean(tCHL,2,'omitnan');

%% var of anomalies by lme
vsst = var(asst,0,2,'omitnan');
vchl = var(achl,0,2,'omitnan');

%% save
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

save([fpath 'lme_satellite_sst_chl_interann_var.mat'],...
    'lme_sst_std','lme_chl_std',...
    'lme_sst_mean','lme_chl_mean',...
    'lme_sst_ts','lme_chl_ts',...
    'lme_sst_cv','lme_chl_cv','cyr','tyr');

save([fpath 'lme_satellite_sst_chl_ann_mean_anoms.mat'],...
    'asst','achl','vsst','vchl','cyr','tyr');

save([cpath 'lme_satellite_chl_interann_var.mat'],...
    'lme_chl_std','lme_chl_mean','lme_chl_ts','lme_chl_cv','cyr');

save([cpath 'lme_satellite_chl_ann_mean_anoms.mat'],...
    'achl','vchl','cyr');

save([tpath 'lme_satellite_sst_interann_var.mat'],...
    'lme_sst_std','lme_sst_mean','lme_sst_ts','lme_sst_cv','tyr');

save([tpath 'lme_satellite_sst_ann_mean_anoms.mat'],...
    'asst','vsst','tyr');

