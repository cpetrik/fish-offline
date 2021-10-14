% FEISTY forced by FOSI
% A param using mzpref S=1, M=0.5
% initialized with spinup pref=1
% mean of all 68 yrs saved

clear all
close all

cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6.mat']);
load([cpath 'Data_grid_POP_gx1v6.mat']);

mod = 'v13_All_fish03';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
spath=['/Volumes/MIP/NC/CESM_MAPP/'];

%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
ppath = [pp 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_sMZ100_mMZ050_nmort1_BE08_noCC_RE00100/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
%%
mz = 0.55:0.05:0.7;
nid = GRD.N;

SF = NaN*ones(nid,length(mz));
SP = SF;
SD = SF;
MF = SF;
MP = SF;
MD = SF;
LP = SF;
LD = SF;
BI = SF;
MFc = SF;
MPc = SF;
MDc = SF;
LPc = SF;
LDc = SF;
lme_Fmcatch = NaN*ones(66,length(mz));
lme_Pmcatch = NaN*ones(66,length(mz));
lme_Dmcatch = NaN*ones(66,length(mz));
lme_AllF = NaN*ones(66,length(mz));
lme_AllP = NaN*ones(66,length(mz));
r_all = NaN*ones(5,length(mz));
ss_all = NaN*ones(5,length(mz));
rmse_all = NaN*ones(5,length(mz));
mis_all = NaN*ones(length(mz),45,5);

%%
for M=1:length(mz)
    Apar = mz(M);
    
    % Equal
    tap = num2str(1000+int64(100*Apar));
    sname = ['Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A',tap(2:end),...
        '_sMZ100_mMZ050_nmort1_BE08_noCC_RE00100'];
    load([spath sname '/FOSI_' mod '_means.mat']);
    
    %% Means
    sf_mean = FOSI_Sml_f_bio;
    sp_mean = FOSI_Sml_p_bio;
    sd_mean = FOSI_Sml_d_bio;
    mf_mean = FOSI_Med_f_bio;
    mp_mean = FOSI_Med_p_bio;
    md_mean = FOSI_Med_d_bio;
    lp_mean = FOSI_Lrg_p_bio;
    ld_mean = FOSI_Lrg_d_bio;
    b_mean = FOSI_Bent_bio;
    
    mf_my = FOSI_Med_f_yield;
    mp_my = FOSI_Med_p_yield;
    md_my = FOSI_Med_d_yield;
    lp_my = FOSI_Lrg_p_yield;
    ld_my = FOSI_Lrg_d_yield;
    
    SF(:,M) = FOSI_Sml_f_bio;
    SP(:,M) = FOSI_Sml_p_bio;
    SD(:,M) = FOSI_Sml_d_bio;
    MF(:,M) = FOSI_Med_f_bio;
    MP(:,M) = FOSI_Med_p_bio;
    MD(:,M) = FOSI_Med_d_bio;
    LP(:,M) = FOSI_Lrg_p_bio;
    LD(:,M) = FOSI_Lrg_d_bio;
    BI(:,M) = FOSI_Bent_bio;
    
    MFc(:,M) = FOSI_Med_f_yield;
    MPc(:,M) = FOSI_Med_p_yield;
    MDc(:,M) = FOSI_Med_d_yield;
    LPc(:,M) = FOSI_Lrg_p_yield;
    LDc(:,M) = FOSI_Lrg_d_yield;
    
    %% Maps
    vis_fosi_ensem(ppath,sname,sf_mean,sp_mean,sd_mean,...
    mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean);
    
    %% LME
    [lme_mcatch,lme_mbio,lme_area] = lme_fosi_ensem(sf_mean,sp_mean,sd_mean,...
    mf_mean,mp_mean,md_mean,b_mean,lp_mean,ld_mean,mf_my,mp_my,md_my,...
    lp_my,ld_my);
    
    lme_Fmcatch(:,M) = lme_mcatch(:,1);
    lme_Pmcatch(:,M) = (lme_mcatch(:,2)+lme_mcatch(:,4));
    lme_Dmcatch(:,M) = (lme_mcatch(:,3)+lme_mcatch(:,5));
    
    lme_AllF(:,M) = lme_mbio(:,1)+lme_mbio(:,4);
    lme_AllP(:,M) = lme_mbio(:,2)+lme_mbio(:,5)+lme_mbio(:,7);
    
    %% SAU comparison
    % ADD CALC OF RESIDUALS
    [r,rmse,ss,mis] = lme_saup_corr_stock_ensem(lme_mcatch,lme_area);
    r_all(:,M) = r;
    rmse_all(:,M) = rmse;
    ss_all(:,M) = ss;
    mis_all(M,:,:) = mis;
    
    %% Save
    save([spath sname '/FOSI_' mod '_means.mat'],...
        'lme_mcatch','lme_mbio','lme_area','-append');
    
end

%%
% [spath 'FOSI_v13_ensemble_mzpref_equal.mat']
% [spath 'FOSI_v13_ensemble_mzpref_MhalfS.mat']
% [spath 'FOSI_v13_ensemble_mzpref_ShalfM.mat']
% [spath 'FOSI_v13_ensemble_mzpref_S1_Mmzpref.mat']
% [spath 'FOSI_v13_ensemble_mzpref_M1_Smzpref.mat']
% [spath 'FOSI_v13_ensemble_mzpref_SMsum1.mat']

save([spath 'FOSI_v13_ensemble_Apref_MhalfS.mat'],'SF','SP','SD',...
    'MF','MP','MD','BI','LP','LD','MFc','MPc','MDc','LPc','LDc',...
    'lme_Fmcatch','lme_Pmcatch','lme_Dmcatch','lme_AllF','lme_AllP',...
    'r_all','rmse_all','ss_all','mis_all','mz')

