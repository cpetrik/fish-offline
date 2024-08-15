% CESM FOSI output
% calc interann variability by lme

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
    'dety','zlosy') 
%also annual means to be like fish

[ni,nj,nyr] = size(dety);

%% lme area-weighted means
tlme = double(lme_mask);
tlme(tlme<0) = nan;
olme = tlme(GRD.ID);

lme_dety_mean = NaN*ones(66,nyr);
lme_mzly_mean = NaN*ones(66,nyr);

%vectorize
dety = reshape(dety,ni*nj,nyr);
mzly = reshape(zlosy,ni*nj,nyr);
area_vec = reshape(AREA_OCN,ni*nj,1);
area_vec = repmat(area_vec,1,nyr);

for L=1:66
    lid = find(tlme==L);
    if ~isempty(lid)
        
        lme_dety_mean(L,:) = (sum(dety(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
        lme_mzly_mean(L,:) = (sum(mzly(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
        
    end
end

lme_mzl_det = (lme_mzly_mean+eps) ./ (lme_dety_mean+eps);

%% remove linear trend
xzd = NaN*ones(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62
    % Drivers
    xi = lme_mzl_det(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        xzd(i,:) = dR;
    end
    clear R T t b m tH dR data
end

%% anomalies
azld = xzd - mean(xzd,2,'omitnan');

% std by lme before putting back on grid
lme_zld_stda = std(azld,0,2,'omitnan');

%% save anoms
save([fpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'azld','lme_zld_stda','-append');

