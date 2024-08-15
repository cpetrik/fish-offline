% Calc single linear regression of forcing-fish for LMEs

clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";

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
manom(:,:,4) = azlosy;

tanom = {'TP','TB','Det','ZmLoss'};

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

% Anoms with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'_lme_ann_mean_anoms.mat'],...
    'aa','ab','ad','af','ap')

%% Divide anoms by 2SD
sdm = std(manom,0,2);
sdm2 = repmat(sdm,1,68,1);
manom = manom ./ (2*sdm2);

%Fish
af = af ./ (2*std(af,0,2));
ap = ap ./ (2*std(ap,0,2));
ad = ad ./ (2*std(ad,0,2));
aa = aa ./ (2*std(aa,0,2));
ab = ab ./ (2*std(ab,0,2));

%% % Corr of forcing ---------------------------------------------------------
yst = 1;
yen = 68;

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = aa(:,1);
lid = find(~isnan(AA));

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

FtabC = nan*ones(length(tanom),length(yr),length(lid));
FtabP = nan*ones(length(tanom),length(yr),length(lid));
PtabC = FtabC;
PtabP = FtabC;
DtabC = FtabC;
DtabP = FtabC;
AtabC = FtabC;
AtabP = FtabC;
BtabC = FtabC;
BtabP = FtabC;
FtabR2 = FtabC;
PtabR2 = FtabC;
DtabR2 = FtabC;
AtabR2 = FtabC;
BtabR2 = FtabC;

%%
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = num2str(i);

    for j = 1:length(tanom)

        %input forcing
        driver = tanom{j};

        for k=1:length(yr) %Correlation at diff lags
            t = yr(k);

            %               LME time    driver       
            sclim = ((manom(i,yst:yen-t,j))') ;

            %Fish
            mdl0 = fitlm(sclim , (af(i,yst+t:yen))');
            FtabC(j,k,L) = mdl0.Coefficients.Estimate(2);
            FtabP(j,k,L) = mdl0.Coefficients.pValue(2);
            FtabR2(j,k,L) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(sclim , (ap(i,yst+t:yen))');
            PtabC(j,k,L) = mdl0.Coefficients.Estimate(2);
            PtabP(j,k,L) = mdl0.Coefficients.pValue(2);
            PtabR2(j,k,L) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(sclim , (ad(i,yst+t:yen))');
            DtabC(j,k,L) = mdl0.Coefficients.Estimate(2);
            DtabP(j,k,L) = mdl0.Coefficients.pValue(2);
            DtabR2(j,k,L) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(sclim , (aa(i,yst+t:yen))');
            AtabC(j,k,L) = mdl0.Coefficients.Estimate(2);
            AtabP(j,k,L) = mdl0.Coefficients.pValue(2);
            AtabR2(j,k,L) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(sclim , (ab(i,yst+t:yen))');
            BtabC(j,k,L) = mdl0.Coefficients.Estimate(2);
            BtabP(j,k,L) = mdl0.Coefficients.pValue(2);
            BtabR2(j,k,L) = mdl0.Rsquared.Ordinary;
            clear mdl0

        end % time lag

    end % driver

end %LME

%%
save([spath,'LMEs_SLR_biom_drivers_0_5_lag.mat'],'tanom','lid','yr',...
    'FtabC','FtabP','PtabC','PtabP','DtabC','DtabP','BtabC','BtabP',...
    'AtabC','AtabP','FtabR2','PtabR2','DtabR2','AtabR2','BtabR2');

