% CESM FOSI output
% calc interann variability by lme

clear all
close all

%% Paths
fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
%fpath='/Volumes/petrik-lab/GCM_Data/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6_noSeas.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');
load([fpath 'LME-mask-POP_gx1v6.mat']);

%% FEISTY Inputs
fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
%fpath='/Volumes/petrik-lab/GCM_Data/CESM/FOSI/';
load([fpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat'],...
    'tp','tb','det','zoo','zlos')

[ni,nj,nyr] = size(det);

%% lme means
tlme = double(lme_mask);
tlme(tlme<0) = nan;
olme = tlme(GRD.ID);

lme_tp_mean = NaN*ones(66,nyr);
lme_tb_mean = NaN*ones(66,nyr);
lme_det_mean = NaN*ones(66,nyr);
lme_mz_mean = NaN*ones(66,nyr);
lme_mzl_mean = NaN*ones(66,nyr);

%vectorize
tp2 = reshape(tp,ni*nj,nyr);
tb2 = reshape(tb,ni*nj,nyr);
det2 = reshape(det,ni*nj,nyr);
mz = reshape(zoo,ni*nj,nyr);
mzl = reshape(zlos,ni*nj,nyr);

for L=1:66
    lid = find(tlme==L);
    if ~isempty(lid)
        
        ltp = tp2(lid,:);
        ltb = tb2(lid,:);
        ldet = det2(lid,:);
        lmz = mz(lid,:);
        lmzl = mzl(lid,:);
        
        lme_tp_mean(L,:) = nanmean(ltp);
        lme_tb_mean(L,:) = nanmean(ltb);
        lme_det_mean(L,:) = nanmean(ldet);
        lme_mz_mean(L,:) = nanmean(lmz);
        lme_mzl_mean(L,:) = nanmean(lmzl);
        
    end
end

%% remove linear trend
xtp = NaN*ones(66,nyr);
xtb = NaN*ones(66,nyr);
xdet = NaN*ones(66,nyr);
xz = NaN*ones(66,nyr);
xzl = NaN*ones(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62
    % Drivers
    xi = lme_tp_mean(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); %0.651784 seconds.
    tH = m*t + b;
    dR = R - tH;
    xtp(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = lme_tb_mean(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data);
    tH = m*t + b;
    dR = R - tH;
    xtb(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = lme_det_mean(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data);
    tH = m*t + b;
    dR = R - tH;
    xdet(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = lme_mz_mean(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        xz(i,:) = dR;
    end
    clear R T t b m tH dR data
    
    xi = lme_mzl_mean(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        xzl(i,:) = dR;
    end
    clear R T t b m tH dR data
    
end

%% anomalies
atp = xtp - nanmean(xtp,2);
atb = xtb - nanmean(xtb,2);
adet = xdet - nanmean(xdet,2);
azoo = xz - nanmean(xz,2);
azlos = xzl - nanmean(xzl,2);

%% std by lme before putting back on grid
lme_tp_stda = std(atp,0,2,'omitnan');
lme_tb_stda = std(atb,0,2,'omitnan');
lme_det_stda = std(adet,0,2,'omitnan');
lme_mz_stda = std(amz,0,2,'omitnan');
lme_mzl_stda = std(amzl,0,2,'omitnan');


%% save anoms
fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
save([fpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'atp','atb','adet','azoo','azlos',...
    'lme_tp_stda','lme_tb_stda','lme_det_stda','lme_mz_stda','lme_mzl_stda');
