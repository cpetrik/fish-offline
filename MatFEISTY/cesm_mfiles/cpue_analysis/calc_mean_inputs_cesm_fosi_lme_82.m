% CESM FOSI output
% calc interann means by lme

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

fyr = 1948:2015;


%% save means
save([fpath 'CESM_FOSI_v15_lme_interann_mean_drivers_1948_2015.mat'],...
    'lme_tp_mean','lme_tb_mean','lme_dety_mean','lme_mzly_mean','fyr');

