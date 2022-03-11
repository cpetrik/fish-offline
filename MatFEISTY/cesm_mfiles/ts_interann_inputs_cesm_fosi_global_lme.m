% CESM FOSI output
% means by LME over time

clear all
close all

%% Paths

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([fpath 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');
load([fpath 'LME-mask-POP_gx1v6.mat']);

ID = GRD.ID;

% AREA_OCN = TAREA * 1e-4;
% tlme = double(lme_mask);
% tlme(tlme<0) = nan;

%%
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'time','yr');
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo_totloss_allphytoC.mat'],...
    'LzooC_150m','Lzoo_loss_150m');

%% nans
TEMP_bottom(TEMP_bottom >= 9.9e+36) = nan;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;

% zeros
LzooC_150m(LzooC_150m<0) = 0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0;

% doubles
tp = double(TEMP_150m);
tb = double(TEMP_bottom);
det_btm = double(POC_FLUX_IN_bottom);
mz = double(LzooC_150m);
loss_mz = double(Lzoo_loss_150m);

%% reshape and vectorize to ocean cells
[ni,nj,nt] = size(tp);

tp = reshape(tp,ni*nj,nt);
tb = reshape(tb,ni*nj,nt);
det_btm = reshape(det_btm,ni*nj,nt);
mz = reshape(mz,ni*nj,nt);
loss_mz = reshape(loss_mz,ni*nj,nt);

% % g/m2 --> total g
% AREA = TAREA(ID) * 1e-4;
% AREA_OCN = repmat(AREA,1,nt);
% tp_mean = tp(ID,:) .* AREA_OCN;
% tb_mean = tb(ID,:) .* AREA_OCN;
% det_mean = det_btm(ID,:) .* AREA_OCN;
% mz_mean = mz(ID,:) .* AREA_OCN;
% loss_mean = loss_mz(ID,:) .* AREA_OCN;

%% 
tp_mean = tp(ID,:);
tb_mean = tb(ID,:);
det_mean = det_btm(ID,:);
mz_mean = mz(ID,:);
loss_mean = loss_mz(ID,:);

%% Calc LMEs
glme = double(lme_mask);
glme(glme<0) = nan;
tlme = glme(ID);

lme_mtp = NaN*ones(66,nt);
lme_mtb = lme_mtp;
lme_mdet = lme_mtp;
lme_mmz = lme_mtp;
lme_mloss = lme_mtp;

for L=1:66
    lid = find(tlme==L);
    if (~isempty(lid))
        %mean biomass
        lme_mtp(L,:) = nanmean(tp_mean(lid,:));
        lme_mtb(L,:) = nanmean(tb_mean(lid,:));
        lme_mdet(L,:) = nanmean(det_mean(lid,:));
        lme_mmz(L,:) = nanmean(mz_mean(lid,:));
        lme_mloss(L,:) = nanmean(loss_mean(lid,:));
    end
end

%%
save([fpath 'FOSI_lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'],...
    'lme_mtp','lme_mtb','lme_mdet','lme_mmz','lme_mloss','tlme','ID');

spath='/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/FOSI/';
if (~isfolder(spath))
    mkdir(spath)
end
save([spath 'FOSI_lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'],...
    'lme_mtp','lme_mtb','lme_mdet','lme_mmz','lme_mloss','tlme','ID');
