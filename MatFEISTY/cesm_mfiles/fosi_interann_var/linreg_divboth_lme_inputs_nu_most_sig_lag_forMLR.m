% Find most signif lag in single linear regression
% For each driver and fn type - nu
% To use in mult linear regression
% divide both ts by 2SD


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
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms.mat'],...
    'aa','ad','af','ap')

%% Divide anoms by 2SD
sdm = std(manom,0,2);
sdm2 = repmat(sdm,1,68,1);
manom = manom ./ (2*sdm2);

%Fish
af = af ./ (2*std(af,0,2));
ap = ap ./ (2*std(ap,0,2));
ad = ad ./ (2*std(ad,0,2));
aa = aa ./ (2*std(aa,0,2));

%% % LM of forcing ---------------------------------------------------------
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

LFtab = nan*ones(length(lid),5);
LPtab = nan*ones(length(lid),5);
LDtab = nan*ones(length(lid),5);
LAtab = nan*ones(length(lid),5);

FtabC = nan*ones(length(tanom),length(yr));
FtabP = nan*ones(length(tanom),length(yr));
PtabC = FtabC;
PtabP = FtabC;
DtabC = FtabC;
DtabP = FtabC;
AtabC = FtabC;
AtabP = FtabC;

[Ymat,Jmat] = meshgrid(yr,1:length(tanom));
Ymat = Ymat';
Jmat = Jmat';

%%
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = num2str(i);

    for j = 1:length(tanom)

        %input forcing
        driver = tanom{j};

        for k=1:length(yr) %Linear regression at diff lags
            t = yr(k);

            %              LME  time    driver                
            sclim = ((manom(i,yst:yen-t,j))') ;

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

        end % time lag

    end % driver
    %%
    %save([spath,ilme,'_regress_drivers_nu_ALLdiv2SD_0_5_lag.mat'])

    %%
    maxC = max(abs(FtabC'));
    cid = find(abs(FtabC')==maxC);
    LFtab(L,:) = Ymat(cid);

    maxC = max(abs(PtabC'));
    cid = find(abs(PtabC')==maxC);
    LPtab(L,:) = Ymat(cid);
    
    maxC = max(abs(DtabC'));
    cid = find(abs(DtabC')==maxC);
    LDtab(L,:) = Ymat(cid);
    
    maxC = max(abs(AtabC'));
    cid = find(abs(AtabC')==maxC);
    LAtab(L,:) = Ymat(cid);
    
end %LME

%%
lname = cellstr(num2str(lid));
% cnam, lname, tanom2
Ftab1 = array2table(LFtab,"RowNames",lname);
Ftab1.Properties.VariableNames = tanom;

Ptab1 = array2table(LPtab,"RowNames",lname);
Ptab1.Properties.VariableNames = tanom;

Dtab1 = array2table(LDtab,"RowNames",lname);
Dtab1.Properties.VariableNames = tanom;

Atab1 = array2table(LAtab,"RowNames",lname);
Atab1.Properties.VariableNames = tanom;

%%
writetable(Ftab1,[spath,'LMEs_regress_drivers_ALLdiv2SD_siglag_Fnu.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'LMEs_regress_drivers_ALLdiv2SD_siglag_Pnu.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'LMEs_regress_drivers_ALLdiv2SD_siglag_Dnu.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Atab1,[spath,'LMEs_regress_drivers_ALLdiv2SD_siglag_Anu.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'LMEs_nu_regress_drivers_ALLdiv2SD_siglags.mat'],...
    'LFtab','LPtab','LDtab','LAtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','lid');

