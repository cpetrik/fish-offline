% Find patterns in forcing-fish correlations

clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";

%% % ------------------------------------------------------------
% load data
apath = '/Users/cpetrik/Dropbox/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'Climate_anomalies_seasonal_means.mat']);

%% FOSI input forcing
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adet','atb','atp','azlos','azoo');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

% Anoms with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'_lme_ann_mean_anoms.mat'],...
    'aa','ab','ad','af','ap','as','am','al')

%% LM of forcing ---------------------------------------------------------
cnam = {'coef','p','lag','iclimate','climate'};

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = aa(:,1);
lid = find(~isnan(AA));

%Lags
yr = 0:5;

% Climate
jid = 1:5;
tanom2=tanom(jid)';
tanom2(:,2)=tanom2(:,1);
tanom2(:,3)=tanom2(:,1);
tanom2(:,4)=tanom2(:,1);
tanom2(:,5)=tanom2(:,1);
tanom2(:,6)=tanom2(:,1);

[Ymat,Jmat] = meshgrid(yr,1:length(jid));

LFtab = nan*ones(length(lid),4);
LPtab = nan*ones(length(lid),4);
LDtab = nan*ones(length(lid),4);
LAtab = nan*ones(length(lid),4);
LBtab = nan*ones(length(lid),4);
LTPtab = nan*ones(length(lid),4);
LTBtab = nan*ones(length(lid),4);
LPOCtab = nan*ones(length(lid),4);
LZMtab = nan*ones(length(lid),4);
LZLtab = nan*ones(length(lid),4);

LFt = cell(length(lid),1);
LPt = cell(length(lid),1);
LDt = cell(length(lid),1);
LAt = cell(length(lid),1);
LBt = cell(length(lid),1);
LTPt = cell(length(lid),1);
LTBt = cell(length(lid),1);
LPOCt = cell(length(lid),1);
LZMt = cell(length(lid),1);
LZLt = cell(length(lid),1);

FtabC = nan*ones(length(jid),length(yr));
FtabP = nan*ones(length(jid),length(yr));
PtabC = FtabC;
PtabP = FtabC;
DtabC = FtabC;
DtabP = FtabC;
AtabC = FtabC;
AtabP = FtabC;
BtabC = FtabC;
BtabP = FtabC;
TPtabC = FtabC;
TPtabP = FtabC;
TBtabC = FtabC;
TBtabP = FtabC;
POCtabC = FtabC;
POCtabP = FtabC;
ZMtabC = FtabC;
ZMtabP = FtabC;
ZLtabC = FtabC;
ZLtabP = FtabC;

%% Loop over climate indices and responses
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = num2str(i);

    %%
    for J = 1:length(jid)
        j = jid(J);

        %climate
        driver = tanom{j};

        for k=1:length(yr) %Linear regression at diff lags

            t = yr(k);

            sclim = ((manom(j,yst(j):yen(j)-t))') ./ (2*std((manom(j,yst(j):yen(j)-t))));

            %Fish
            mdl0 = fitlm(sclim , (af(i,yst(j)+t:yen(j)))');
            FtabC(J,k) = mdl0.Coefficients.Estimate(2);
            FtabP(J,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (ap(i,yst(j)+t:yen(j)))');
            PtabC(J,k) = mdl0.Coefficients.Estimate(2);
            PtabP(J,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (ad(i,yst(j)+t:yen(j)))');
            DtabC(J,k) = mdl0.Coefficients.Estimate(2);
            DtabP(J,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (aa(i,yst(j)+t:yen(j)))');
            AtabC(J,k) = mdl0.Coefficients.Estimate(2);
            AtabP(J,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (ab(i,yst(j)+t:yen(j)))');
            BtabC(J,k) = mdl0.Coefficients.Estimate(2);
            BtabP(J,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            % Inputs / Forcing
            mdl0 = fitlm(sclim , (atp(i,yst(j)+t:yen(j)))');
            TPtabC(J,k) = mdl0.Coefficients.Estimate(2);
            TPtabP(J,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (atb(i,yst(j)+t:yen(j)))');
            TBtabC(J,k) = mdl0.Coefficients.Estimate(2);
            TBtabP(J,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (adet(i,yst(j)+t:yen(j)))');
            POCtabC(J,k) = mdl0.Coefficients.Estimate(2);
            POCtabP(J,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (azoo(i,yst(j)+t:yen(j)))');
            ZMtabC(J,k) = mdl0.Coefficients.Estimate(2);
            ZMtabP(J,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (azlos(i,yst(j)+t:yen(j)))');
            ZLtabC(J,k) = mdl0.Coefficients.Estimate(2);
            ZLtabP(J,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

        end % Lag

    end % Climate

    %%
    %save([spath,ilme,'_regress_seasonal_climate5_div2SD_0_5_lag.mat'])

    %%
    maxC = max(abs(FtabC(:)));
    cid = find(abs(FtabC(:))==maxC);
    LFtab(L,1) = FtabC(cid);
    LFtab(L,2) = FtabP(cid);
    LFtab(L,3) = Ymat(cid);
    LFtab(L,4) = Jmat(cid);
    LFt(L) = tanom2(cid);

    maxC = max(abs(PtabC(:)));
    cid = find(abs(PtabC(:))==maxC);
    LPtab(L,1) = PtabC(cid);
    LPtab(L,2) = PtabP(cid);
    LPtab(L,3) = Ymat(cid);
    LPtab(L,4) = Jmat(cid);
    LPt(L) = tanom2(cid);

    maxC = max(abs(DtabC(:)));
    cid = find(abs(DtabC(:))==maxC);
    LDtab(L,1) = DtabC(cid);
    LDtab(L,2) = DtabP(cid);
    LDtab(L,3) = Ymat(cid);
    LDtab(L,4) = Jmat(cid);
    LDt(L) = tanom2(cid);

    maxC = max(abs(AtabC(:)));
    cid = find(abs(AtabC(:))==maxC);
    LAtab(L,1) = AtabC(cid);
    LAtab(L,2) = AtabP(cid);
    LAtab(L,3) = Ymat(cid);
    LAtab(L,4) = Jmat(cid);
    LAt(L) = tanom2(cid);

    maxC = max(abs(BtabC(:)));
    cid = find(abs(BtabC(:))==maxC);
    LBtab(L,1) = BtabC(cid);
    LBtab(L,2) = BtabP(cid);
    LBtab(L,3) = Ymat(cid);
    LBtab(L,4) = Jmat(cid);
    LBt(L) = tanom2(cid);

    maxC = max(abs(TPtabC(:)));
    cid = find(abs(TPtabC(:))==maxC);
    LTPtab(L,1) = TPtabC(cid);
    LTPtab(L,2) = TPtabP(cid);
    LTPtab(L,3) = Ymat(cid);
    LTPtab(L,4) = Jmat(cid);
    LTPt(L) = tanom2(cid);

    maxC = max(abs(TBtabC(:)));
    cid = find(abs(TBtabC(:))==maxC);
    LTBtab(L,1) = TBtabC(cid);
    LTBtab(L,2) = TBtabP(cid);
    LTBtab(L,3) = Ymat(cid);
    LTBtab(L,4) = Jmat(cid);
    LTBt(L) = tanom2(cid);

    maxC = max(abs(POCtabC(:)));
    cid = find(abs(POCtabC(:))==maxC);
    LPOCtab(L,1) = POCtabC(cid);
    LPOCtab(L,2) = POCtabP(cid);
    LPOCtab(L,3) = Ymat(cid);
    LPOCtab(L,4) = Jmat(cid);
    LPOCt(L) = tanom2(cid);

    maxC = max(abs(ZMtabC(:)));
    cid = find(abs(ZMtabC(:))==maxC);
    LZMtab(L,1) = ZMtabC(cid);
    LZMtab(L,2) = ZMtabP(cid);
    LZMtab(L,3) = Ymat(cid);
    LZMtab(L,4) = Jmat(cid);
    LZMt(L) = tanom2(cid);

    maxC = max(abs(ZLtabC(:)));
    cid = find(abs(ZLtabC(:))==maxC);
    LZLtab(L,1) = ZLtabC(cid);
    LZLtab(L,2) = ZLtabP(cid);
    LZLtab(L,3) = Ymat(cid);
    LZLtab(L,4) = Jmat(cid);
    LZLt(L) = tanom2(cid);

end %LME

%%
lname = cellstr(num2str(lid));
% cnam, lname, tanom2
Ftab1 = array2table(LFtab,"RowNames",lname);
Ftab1(:,5) = LFt;
Ftab1.Properties.VariableNames = cnam;

Ptab1 = array2table(LPtab,"RowNames",lname);
Ptab1(:,5) = LPt;
Ptab1.Properties.VariableNames = cnam;

Dtab1 = array2table(LDtab,"RowNames",lname);
Dtab1(:,5) = LDt;
Dtab1.Properties.VariableNames = cnam;

Atab1 = array2table(LAtab,"RowNames",lname);
Atab1(:,5) = LAt;
Atab1.Properties.VariableNames = cnam;

Btab1 = array2table(LBtab,"RowNames",lname);
Btab1(:,5) = LBt;
Btab1.Properties.VariableNames = cnam;

TPtab1 = array2table(LTPtab,"RowNames",lname);
TPtab1(:,5) = LTPt;
TPtab1.Properties.VariableNames = cnam;

TBtab1 = array2table(LTBtab,"RowNames",lname);
TBtab1(:,5) = LTBt;
TBtab1.Properties.VariableNames = cnam;

POCtab1 = array2table(LPOCtab,"RowNames",lname);
POCtab1(:,5) = LPOCt;
POCtab1.Properties.VariableNames = cnam;

ZMtab1 = array2table(LZMtab,"RowNames",lname);
ZMtab1(:,5) = LZMt;
ZMtab1.Properties.VariableNames = cnam;

ZLtab1 = array2table(LZLtab,"RowNames",lname);
ZLtab1(:,5) = LZLt;
ZLtab1.Properties.VariableNames = cnam;

%%
writetable(Ftab1,[spath,'LMEs_regress_seasonal_climate_div2SD_maxcoef_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'LMEs_regress_seasonal_climate_div2SD_maxcoef_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'LMEs_regress_seasonal_climate_div2SD_maxcoef_D.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Atab1,[spath,'LMEs_regress_seasonal_climate_div2SD_maxcoef_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Btab1,[spath,'LMEs_regress_seasonal_climate_div2SD_maxcoef_B.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(TPtab1,[spath,'LMEs_regress_seasonal_climate_div2SD_maxcoef_TP.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(TBtab1,[spath,'LMEs_regress_seasonal_climate_div2SD_maxcoef_TB.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(POCtab1,[spath,'LMEs_regress_seasonal_climate_div2SD_maxcoef_POC.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(ZMtab1,[spath,'LMEs_regress_seasonal_climate_div2SD_maxcoef_ZM.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(ZLtab1,[spath,'LMEs_regress_seasonal_climate_div2SD_maxcoef_ZL.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'LMEs_regress_seasonal_climate_div2SD_maxcoefs.mat'],...
    'LFtab','LPtab','LDtab','LAtab','LBtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','Btab1',...
    'TPtab1','TBtab1','POCtab1','ZMtab1','ZLtab1',...
    'LTPtab','LTBtab','LPOCtab','LZMtab','LZLtab','lid');


%% Plots














