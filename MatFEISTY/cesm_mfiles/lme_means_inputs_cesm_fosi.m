% CESM FOSI output
% means by LME

clear all
close all

%% Paths

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([fpath 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');
load([fpath 'LME-mask-POP_gx1v6.mat']);

AREA_OCN = TAREA * 1e-4;
tlme = double(lme_mask);
tlme(tlme<0) = nan;

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

%zeros
LzooC_150m(LzooC_150m<0) = 0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0;

%%
tp = double(TEMP_150m);
tb = double(TEMP_bottom);
det_btm = double(POC_FLUX_IN_bottom) * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
mz = double(LzooC_150m) * 1e-9 * 1e4 * 12.01 * 9.0;
hploss_mz = double(Lzoo_loss_150m) * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

% means in space
lme_tp_fosi = NaN*ones(66,1);
lme_tb_fosi = NaN*ones(66,1);
lme_det_fosi = NaN*ones(66,1);
lme_mz_fosi = NaN*ones(66,1);
lme_mzloss_fosi = NaN*ones(66,1);
lme_area = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);

    lme_area(L,1) = nansum(AREA_OCN(lid));
    lme_tp_fosi(L,1) = nanmean(tp(lid));
    lme_tb_fosi(L,1) = nanmean(tb(lid));
    lme_det_fosi(L,1) = nanmean(det_btm(lid));
    lme_mz_fosi(L,1) = nanmean(mz(lid));
    lme_mzloss_fosi(L,1) = nanmean(hploss_mz(lid));

end

%% save
save([fpath 'lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'],'lme_area',...
    'lme_tp_fosi','lme_tb_fosi','lme_det_fosi','lme_mz_fosi','lme_mzloss_fosi')
