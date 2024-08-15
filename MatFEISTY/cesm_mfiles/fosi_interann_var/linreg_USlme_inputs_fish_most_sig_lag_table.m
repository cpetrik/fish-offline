% Find patterns in forcing-fish correlations

clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";

% LMEs
lid = [54,1:2,65,10,3,5:7]; %ADD 65 = Aleutian Islands
lname = {'CHK','EBS','GAK','AI','HI','CCE','GMX','SE','NE'};

%% FOSI input forcing

%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adet','adety','atb','atp','azlos','azlosy','azoo');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

% put into a matrix & use annual production
manom(:,:,1) = atp;
manom(:,:,2) = atb;
manom(:,:,3) = adety;
manom(:,:,4) = azoo;
manom(:,:,5) = azlosy;

tanom = {'TP','TB','Det','Zmeso','ZmLoss'};

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

%% % LM of forcing ---------------------------------------------------------
yst = 1;
yen = 68;

cnam = {'coef','p','lag','idriver','driver'};

% LMEs
lid = [54,1:2,65,10,3,5:7]; %ADD 65 = Aleutian Islands
lname = {'CHK','EBS','GAK','AI','HI','CCE','GMX','SE','NE'};

%Lags
yr = 0:5;

% Drivers
tanom2=tanom';
tanom2(:,2)=tanom2(:,1);
tanom2(:,3)=tanom2(:,1);
tanom2(:,4)=tanom2(:,1);
tanom2(:,5)=tanom2(:,1);
tanom2(:,6)=tanom2(:,1);

[Ymat,Jmat] = meshgrid(yr,1:length(tanom));

LFtab = nan*ones(length(lid),4);
LPtab = nan*ones(length(lid),4);
LDtab = nan*ones(length(lid),4);
LAtab = nan*ones(length(lid),4);
LBtab = nan*ones(length(lid),4);

LFt = cell(length(lid),1);
LPt = cell(length(lid),1);
LDt = cell(length(lid),1);
LAt = cell(length(lid),1);
LBt = cell(length(lid),1);

FtabC = nan*ones(length(tanom),length(yr));
FtabP = nan*ones(length(tanom),length(yr));
PtabC = FtabC;
PtabP = FtabC;
DtabC = FtabC;
DtabP = FtabC;
AtabC = FtabC;
AtabP = FtabC;
BtabC = FtabC;
BtabP = FtabC;

%%
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = lname{L};

    for j = 1:length(tanom)

        %input forcing
        driver = tanom{j};

        for k=1:length(yr) %Linear regression at diff lags
            t = yr(k);

            %               LME     time      driver                LME     time      driver
            sclim = ((manom(i,yst:yen-t,j))') ./ (2*std((manom(i,yst:yen-t,j))));

            %Fish
            mdl0 = fitlm(sclim , (af(i,yst+t:yen))');
            FtabC(j,k) = mdl0.Coefficients.Estimate(2);
            FtabP(j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (ap(i,yst+t:yen))');
            PtabC(j,k) = mdl0.Coefficients.Estimate(2);
            PtabP(j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (ad(i,yst+t:yen))');
            DtabC(j,k) = mdl0.Coefficients.Estimate(2);
            DtabP(j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (aa(i,yst+t:yen))');
            AtabC(j,k) = mdl0.Coefficients.Estimate(2);
            AtabP(j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (ab(i,yst+t:yen))');
            BtabC(j,k) = mdl0.Coefficients.Estimate(2);
            BtabP(j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0
        end % time lag

    end % driver
    %%
    save([spath,ilme,'_regress_drivers_div2SD_0_5_lag.mat'])

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

end %LME

%%
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

%%
writetable(Ftab1,[spath,'USLMEs_regress_driver_div2SD_maxcoef_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'USLMEs_regress_driver_div2SD_maxcoef_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'USLMEs_regress_driver_div2SD_maxcoef_D.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Atab1,[spath,'USLMEs_regress_driver_div2SD_maxcoef_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Btab1,[spath,'USLMEs_regress_driver_div2SD_maxcoef_B.csv'],...
    'Delimiter',',','WriteRowNames',true);

