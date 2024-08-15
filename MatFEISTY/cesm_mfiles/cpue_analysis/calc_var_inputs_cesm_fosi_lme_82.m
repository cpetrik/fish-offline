% CESM FOSI output
% calc interann variability by lme
% Calculate anomaly time series for different ranges
% 1982-2010 &
% 1982-2015

clear 
close all

%% Paths
fpath='/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6_noSeas.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');
load([fpath 'LME-mask-POP_gx1v6.mat']);

AREA_OCN = TAREA * 1e-4; %grid cell area in m2

%% FEISTY Inputs
fpath='/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([fpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat'],...
    'tp','tb','dety','zlosy') 
%also annual means to be like fish

[ni,nj,nyr] = size(dety);

%% lme area-weighted means
tlme = double(lme_mask);
tlme(tlme<0) = nan;
olme = tlme(GRD.ID);

lme_tp_mean = NaN*ones(66,nyr);
lme_tb_mean = NaN*ones(66,nyr);
lme_dety_mean = NaN*ones(66,nyr);
lme_mzly_mean = NaN*ones(66,nyr);

%vectorize
tp2 = reshape(tp,ni*nj,nyr);
tb2 = reshape(tb,ni*nj,nyr);
dety = reshape(dety,ni*nj,nyr);
mzly = reshape(zlosy,ni*nj,nyr);
area_vec = reshape(AREA_OCN,ni*nj,1);
area_vec = repmat(area_vec,1,nyr);

for L=1:66
    lid = find(tlme==L);
    if ~isempty(lid)
        
        lme_tp_mean(L,:)  = (sum(tp2(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
        lme_tb_mean(L,:)  = (sum(tb2(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
        lme_dety_mean(L,:) = (sum(dety(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
        lme_mzly_mean(L,:) = (sum(mzly(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
        
    end
end

%% Next select specific yrs
fyr = 1948:2015;
eyr = 1982:2010; % subset effort years
yyr = 1982:2015; % subset catch years

lid=1:66;

[~,fide] = intersect(fyr,eyr);
[~,fidc] = intersect(fyr,yyr);

lme_tp_ts_10 = lme_tp_mean(:,fide);
lme_tp_ts_15 = lme_tp_mean(:,fidc);
lme_tb_ts_10 = lme_tb_mean(:,fide);
lme_tb_ts_15 = lme_tb_mean(:,fidc);
lme_det_ts_10 = lme_dety_mean(:,fide);
lme_det_ts_15 = lme_dety_mean(:,fidc);
lme_mzl_ts_10 = lme_mzly_mean(:,fide);
lme_mzl_ts_15 = lme_mzly_mean(:,fidc);

%% remove linear trend
nte = length(eyr);
nty = length(yyr);
xtp10 = NaN*ones(66,nte);
xtp15 = NaN*ones(66,nty);
xtb10 = NaN*ones(66,nte);
xtb15 = NaN*ones(66,nty);
xdet10 = NaN*ones(66,nte);
xdet15 = NaN*ones(66,nty);
xzl10 = NaN*ones(66,nte);
xzl15 = NaN*ones(66,nty);

for i = 1:66 %in future should remove interior seas 23, 33, 62
    % Drivers
    xi = lme_tp_ts_10(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); %0.651784 seconds.
    tH = m*t + b;
    dR = R - tH;
    xtp10(i,:) = dR;
    clear R T t b m tH dR data

    xi = lme_tp_ts_15(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); %0.651784 seconds.
    tH = m*t + b;
    dR = R - tH;
    xtp15(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = lme_tb_ts_10(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data);
    tH = m*t + b;
    dR = R - tH;
    xtb10(i,:) = dR;
    clear R T t b m tH dR data

    xi = lme_tb_ts_15(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data);
    tH = m*t + b;
    dR = R - tH;
    xtb15(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = lme_det_ts_10(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data);
    tH = m*t + b;
    dR = R - tH;
    xdet10(i,:) = dR;
    clear R T t b m tH dR data

    xi = lme_det_ts_15(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data);
    tH = m*t + b;
    dR = R - tH;
    xdet15(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = lme_mzl_ts_10(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        xzl10(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = lme_mzl_ts_15(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        xzl15(i,:) = dR;
    end
    clear R T t b m tH dR data
    
end

%% anomalies
atp10 = xtp10 - mean(xtp10,2,'omitnan');
atp15 = xtp15 - mean(xtp15,2,'omitnan');
atb10 = xtb10 - mean(xtb10,2,'omitnan');
atb15 = xtb15 - mean(xtb15,2,'omitnan');
adet10 = xdet10 - mean(xdet10,2,'omitnan');
adet15 = xdet15 - mean(xdet15,2,'omitnan');
azlos10 = xzl10 - mean(xzl10,2,'omitnan');
azlos15 = xzl15 - mean(xzl15,2,'omitnan');

%% std by lme before putting back on grid
lme_tp10_stda = std(atp10,0,2,'omitnan');
lme_tp15_stda = std(atp15,0,2,'omitnan');
lme_tb10_stda = std(atb10,0,2,'omitnan');
lme_tb15_stda = std(atb15,0,2,'omitnan');
lme_det10_stda = std(adet10,0,2,'omitnan');
lme_det15_stda = std(adet15,0,2,'omitnan');
lme_mzl10_stda = std(azlos10,0,2,'omitnan');
lme_mzl15_stda = std(azlos15,0,2,'omitnan');

%% save anoms
save([fpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom_1982_2010_2015.mat'],...
    'atp10','atp15','atb10','atb15','adet10','adet15','azlos10','azlos15',...
    'lme_tp10_stda','lme_tp15_stda','lme_tb10_stda','lme_tb15_stda','lme_det10_stda','lme_det15_stda',...
    'lme_mzl10_stda','lme_mzl15_stda');

