% CESM FOSI output
% calc interann mean globally and by lme

clear all
close all

%% Paths
fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6.mat'],'GRD');
load([fpath 'LME-mask-POP_gx1v6.mat']);
load([fpath 'CESM_FOSI_interann_mean_forcings_anom.mat']);

tlme = double(lme_mask);
tlme(tlme<0) = nan;

%% time means 
[ni,nj,nyr] = size(tp);
tp2 = reshape(tp,ni*nj,nyr);
tb2 = reshape(tb,ni*nj,nyr);
det2 = reshape(det,ni*nj,nyr);
mz = reshape(zoo,ni*nj,nyr);
mzl = reshape(zlos,ni*nj,nyr);

%% time mean anomaly
mtp = nanmean(tp2);
mtb = nanmean(tb2);
mdet = nanmean(det2);
mzoo = nanmean(mz);
mzlos = nanmean(mzl);

amtp = mtp - nanmean(mtp);
amtb = mtb - nanmean(mtb);
amdet = mdet - nanmean(mdet);
amzoo = mzoo - nanmean(mzoo);
amzlos = mzlos - nanmean(mzlos);

%% time means in LMEs
ltp  = NaN*ones(66,nyr);
ltb  = NaN*ones(66,nyr);
ldet = NaN*ones(66,nyr);
lmz  = NaN*ones(66,nyr);
lmzl = NaN*ones(66,nyr);

for L=1:66
    lid = find(tlme==L);
    
    ltp(L,:) = nanmean(tp2(lid,:));
    ltb(L,:) = nanmean(tb2(lid,:));
    ldet(L,:) = nanmean(det2(lid,:));
    lmz(L,:) = nanmean(mz(lid,:));
    lmzl(L,:) = mean(mzl(lid,:));
    
end

%% LME anomalies
altp = ltp - nanmean(ltp,2);
altb = ltb - nanmean(ltb,2);
aldet = ldet - nanmean(ldet,2);
alzm = lmz - nanmean(lmz,2);
alzml = lmzl - nanmean(lmzl,2);

%%
save([fpath 'CESM_FOSI_interann_mean_forcings_anom.mat'],...
    'mtp','mtb','mdet','mzoo','mzlos',...
    'amtp','amtb','amdet','amzoo','amzlos',...
    'ltp','ltb','ldet','lmz','lmzl',...
    'altp','altb','aldet','alzm','alzml','-append')


