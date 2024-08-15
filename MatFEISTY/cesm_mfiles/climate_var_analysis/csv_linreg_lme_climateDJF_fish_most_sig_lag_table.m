clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];

%%
% climate
apath = '/Users/cpetrik/Dropbox/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'Climate_anomalies_DJF_means.mat'],'tanom');

% FOSI input forcing
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

% LMEs
lid = [54,1:2,65,10,3,5:7]; %ADD 65 = Aleutian Islands
lname = {'CHK','EBS','GAK','AI','HI','CCE','GMX','SE','NE'};

%Lags
yr = 0:5;
ltex = {'lag0','lag1','lag2','lag3','lag4','lag5'};

% Climate
jid = [1:2,4,6,10];
tanom=tanom(jid)';

%% Loop over climate indices and responses
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = lname{L};
    %%
    load([spath,ilme,'_regressDJF_climate5_div2SD_0_5_lag.mat'])

    %% make into tables
    tanom=tanom(jid)';
    FtC = array2table(FtabC,"RowNames",tanom,'VariableNames',ltex);
    FtP = array2table(FtabP,"RowNames",tanom,'VariableNames',ltex);
    PtC = array2table(PtabC,"RowNames",tanom,'VariableNames',ltex);
    PtP = array2table(PtabP,"RowNames",tanom,'VariableNames',ltex);
    DtC = array2table(DtabC,"RowNames",tanom,'VariableNames',ltex);
    DtP = array2table(DtabP,"RowNames",tanom,'VariableNames',ltex);
    AtC = array2table(AtabC,"RowNames",tanom,'VariableNames',ltex);
    AtP = array2table(AtabP,"RowNames",tanom,'VariableNames',ltex);
    BtC = array2table(BtabC,"RowNames",tanom,'VariableNames',ltex);
    BtP = array2table(BtabP,"RowNames",tanom,'VariableNames',ltex);
    TPtC = array2table(TPtabC,"RowNames",tanom,'VariableNames',ltex);
    TPtP = array2table(TPtabP,"RowNames",tanom,'VariableNames',ltex);
    TBtC = array2table(TBtabC,"RowNames",tanom,'VariableNames',ltex);
    TBtP = array2table(TBtabP,"RowNames",tanom,'VariableNames',ltex);
    DEtC = array2table(POCtabC,"RowNames",tanom,'VariableNames',ltex);
    DEtP = array2table(POCtabP,"RowNames",tanom,'VariableNames',ltex);
    ZMtC = array2table(ZMtabC,"RowNames",tanom,'VariableNames',ltex);
    ZMtP = array2table(ZMtabP,"RowNames",tanom,'VariableNames',ltex);
    ZLtC = array2table(ZLtabC,"RowNames",tanom,'VariableNames',ltex);
    ZLtP = array2table(ZLtabP,"RowNames",tanom,'VariableNames',ltex);

    %%
    
    writetable(FtC,[spath,ilme,'_regressDJF_climate_div2SD_coefs_F.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(PtC,[spath,ilme,'_regressDJF_climate_div2SD_coefs_P.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(DtC,[spath,ilme,'_regressDJF_climate_div2SD_coefs_D.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(AtC,[spath,ilme,'_regressDJF_climate_div2SD_coefs_A.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(BtC,[spath,ilme,'_regressDJF_climate_div2SD_coefs_B.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(TPtC,[spath,ilme,'_regressDJF_climate_div2SD_coefs_TP.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(TBtC,[spath,ilme,'_regressDJF_climate_div2SD_coefs_TB.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(DEtC,[spath,ilme,'_regressDJF_climate_div2SD_coefs_POC.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(ZMtC,[spath,ilme,'_regressDJF_climate_div2SD_coefs_ZM.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(ZLtC,[spath,ilme,'_regressDJF_climate_div2SD_coefs_ZL.csv'],...
        'Delimiter',',','WriteRowNames',true);

    writetable(FtP,[spath,ilme,'_regressDJF_climate_div2SD_pvals_F.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(PtP,[spath,ilme,'_regressDJF_climate_div2SD_pvals_P.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(DtP,[spath,ilme,'_regressDJF_climate_div2SD_pvals_D.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(AtP,[spath,ilme,'_regressDJF_climate_div2SD_pvals_A.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(BtP,[spath,ilme,'_regressDJF_climate_div2SD_pvals_B.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(TPtP,[spath,ilme,'_regressDJF_climate_div2SD_pvals_TP.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(TBtP,[spath,ilme,'_regressDJF_climate_div2SD_pvals_TB.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(DEtP,[spath,ilme,'_regressDJF_climate_div2SD_pvals_POC.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(ZMtP,[spath,ilme,'_regressDJF_climate_div2SD_pvals_ZM.csv'],...
        'Delimiter',',','WriteRowNames',true);
    writetable(ZLtP,[spath,ilme,'_regressDJF_climate_div2SD_pvals_ZL.csv'],...
        'Delimiter',',','WriteRowNames',true);

end