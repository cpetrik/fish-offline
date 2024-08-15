% Satellite chl & SST data
% Calculate anomaly time series for different ranges
% 1997-2010 &
% 1997-2015

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

%% Next select specific yrs
eyr = 1997:2010; % subset effort years
yyr = 1997:2015; % subset catch years

lid=1:66;

[~,cide] = intersect(cyr,eyr);
[~,side] = intersect(tyr,eyr);
[~,cidc] = intersect(cyr,yyr);
[~,sidc] = intersect(tyr,yyr);

lme_sst_ts_10 = lme_sst_ts(:,side);
lme_sst_ts_15 = lme_sst_ts(:,sidc);
lme_chl_ts_10 = lme_chl_ts(:,cide);
lme_chl_ts_15 = lme_chl_ts(:,cidc);

%% Then calc mean & std dev
% Mean
lme_sst10_mean  = mean(lme_sst_ts_10,2,'omitnan');
lme_sst15_mean  = mean(lme_sst_ts_15,2,'omitnan');
lme_chl10_mean  = mean(lme_chl_ts_10,2,'omitnan');
lme_chl15_mean  = mean(lme_chl_ts_15,2,'omitnan');

% Std dev
lme_sst10_std  = std(lme_sst_ts_10,0,2,'omitnan');
lme_sst15_std  = std(lme_sst_ts_15,0,2,'omitnan');
lme_chl10_std  = std(lme_chl_ts_10,0,2,'omitnan');
lme_chl15_std  = std(lme_chl_ts_15,0,2,'omitnan');

% Coefficient of variance
lme_sst10_cv = lme_sst10_std ./ lme_sst10_mean;
lme_sst15_cv = lme_sst15_std ./ lme_sst15_mean;
lme_chl10_cv = lme_chl10_std ./ lme_chl10_mean;
lme_chl15_cv = lme_chl15_std ./ lme_chl15_mean;

%% ANOMALIES -------------------------------------------------

% remove linear trend
tSST10 = NaN*ones(66,length(eyr));
tSST15 = NaN*ones(66,length(yyr));
tCHL10 = NaN*ones(66,length(eyr));
tCHL15 = NaN*ones(66,length(yyr));

for i = 1:66 %in future should remove interior seas 23, 33, 62
    
    %SST
    xi = lme_sst_ts_10(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        tSST10(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = lme_sst_ts_15(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        tSST15(i,:) = dR;
    end
    clear R T t b m tH dR data
   
    %Chl
    xi = lme_chl_ts_10(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        tCHL10(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = lme_chl_ts_15(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        tCHL15(i,:) = dR;
    end
    clear R T t b m tH dR data
   
end

%% anomalies
asst10 = tSST10 - mean(tSST10,2,'omitnan');
asst15 = tSST15 - mean(tSST15,2,'omitnan');
achl10 = tCHL10 - mean(tCHL10,2,'omitnan');
achl15 = tCHL15 - mean(tCHL15,2,'omitnan');

%% var of anomalies by lme
vsst10 = var(asst10,0,2,'omitnan');
vsst15 = var(asst15,0,2,'omitnan');
vchl10 = var(achl10,0,2,'omitnan');
vchl15 = var(achl15,0,2,'omitnan');

%% save
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

save([fpath 'lme_satellite_sst_chl_interann_var_1997_2010_2015.mat'],...
    'lme_sst10_std','lme_sst15_std','lme_chl10_std','lme_chl15_std',...
    'lme_sst10_mean','lme_sst15_mean','lme_chl10_mean','lme_chl15_mean',...
    'lme_sst_ts_10','lme_sst_ts_15','lme_chl_ts_10','lme_chl_ts_15',...
    'lme_sst10_cv','lme_sst15_cv','lme_chl10_cv','lme_chl15_cv','eyr','yyr');

save([fpath 'lme_satellite_sst_chl_ann_mean_anoms_1997_2010_2015.mat'],...
    'asst10','asst15','achl10','achl15','vsst10','vchl10','vsst15','vchl15','eyr','yyr');

save([cpath 'lme_satellite_chl_interann_var_1997_2010_2015.mat'],...
    'lme_chl10_std','lme_chl10_mean','lme_chl_ts_10','lme_chl10_cv','eyr','yyr',...
    'lme_chl15_std','lme_chl15_mean','lme_chl_ts_15','lme_chl15_cv');

save([cpath 'lme_satellite_chl_ann_mean_anoms_1997_2010_2015.mat'],...
    'achl10','achl15','vchl10','vchl15','eyr','yyr');

save([tpath 'lme_satellite_sst_interann_var_1997_2010_2015.mat'],...
    'lme_sst10_std','lme_sst10_mean','lme_sst_ts_10','lme_sst10_cv','eyr','yyr',...
    'lme_sst15_std','lme_sst15_mean','lme_sst_ts_15','lme_sst15_cv');

save([tpath 'lme_satellite_sst_ann_mean_anoms_1997_2010_2015.mat'],...
    'asst10','asst15','vsst10','vsst15','eyr','yyr');

