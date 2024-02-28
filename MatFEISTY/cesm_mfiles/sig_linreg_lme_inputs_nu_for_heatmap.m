% Find patterns in forcing-fish linear regressions
% Biomass
% Divide by 2SD

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

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
rpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];

mod = 'v15_All_fish03';

% Anoms with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms.mat'],...
    'af','ap','ad','aa');

%% put into a matrix 
manom(:,:,1) = atp;
manom(:,:,2) = atb;
manom(:,:,3) = adety;
manom(:,:,4) = azlosy;

tanom = {'TP','TB','Det','ZmLoss'};

%% Divide anoms by 2SD
sdm = std(manom,0,2);
sdm2 = repmat(sdm,1,68,1);
manom = manom ./ (2*sdm2);

%Fish
af = af ./ (2*std(af,0,2));
ap = ap ./ (2*std(ap,0,2));
ad = ad ./ (2*std(ad,0,2));
aa = aa ./ (2*std(aa,0,2));

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = aa(:,1);
lid = find(~isnan(AA));

%% % LM of forcing ---------------------------------------------------------
%Loop over drivers and responses
yr = 0:5;
yst = 1;
yen = 68;

LFtab = nan*ones(length(lid),4);
LPtab = nan*ones(length(lid),4);
LDtab = nan*ones(length(lid),4);
LAtab = nan*ones(length(lid),4);

LFlag = nan*ones(length(lid),4);
LPlag = nan*ones(length(lid),4);
LDlag = nan*ones(length(lid),4);
LAlag = nan*ones(length(lid),4);

FtabC = nan*ones(length(yr),1);
FtabP = nan*ones(length(yr),1);
PtabC = FtabC;
PtabP = FtabC;
DtabC = FtabC;
DtabP = FtabC;
AtabC = FtabC;
AtabP = FtabC;
FtabR2 = FtabC;
PtabR2 = FtabC;
DtabR2 = FtabC;
AtabR2 = FtabC;

for L = 1:length(lid) %LME

    %LME
    i = lid(L);
    ilme = num2str(i);

    for j = 1:length(tanom) %drivers
        driver = tanom{j};

        n = 0;

        for k=1:length(yr) %Linear regression at diff lags
            t = yr(k);

            %               LME  time   driver
            sclim = ((manom(i,yst:yen-t,j))');

            %Fish
            mdl0 = fitlm(sclim, (af(i,yst+t:yen))');
            FtabC(k) = mdl0.Coefficients.Estimate(2);
            FtabP(k) = mdl0.Coefficients.pValue(2);
            FtabR2(k) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(sclim, (ap(i,yst+t:yen))');
            PtabC(k) = mdl0.Coefficients.Estimate(2);
            PtabP(k) = mdl0.Coefficients.pValue(2);
            PtabR2(k) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(sclim, (ad(i,yst+t:yen))');
            DtabC(k) = mdl0.Coefficients.Estimate(2);
            DtabP(k) = mdl0.Coefficients.pValue(2);
            DtabR2(k) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(sclim, (aa(i,yst+t:yen))');
            AtabC(k) = mdl0.Coefficients.Estimate(2);
            AtabP(k) = mdl0.Coefficients.pValue(2);
            AtabR2(k) = mdl0.Rsquared.Ordinary;
            clear mdl0

        end % time lag

        %% Find lag with max corr for that driver
        maxC = max(FtabR2);
        if(~isnan(maxC))
            fid = find(abs(FtabR2)==maxC);
            if (FtabP(fid)<=0.05)
                LFtab(L,j) = FtabC(fid);
                LFlag(L,j) = yr(fid);
            end
        end
        clear maxC

        maxC = max(abs(PtabR2));
        if(~isnan(maxC))
            fid = find(abs(PtabR2)==maxC);
            if (PtabP(fid)<=0.05)
                LPtab(L,j) = PtabC(fid);
                LPlag(L,j) = yr(fid);
            end
        end
        clear maxC

        maxC = max(abs(DtabR2));
        if(~isnan(maxC))
            fid = find(abs(DtabR2)==maxC);
            if (DtabP(fid)<=0.05)
                LDtab(L,j) = DtabC(fid);
                LDlag(L,j) = yr(fid);
            end
        end
        clear maxC

        maxC = max(abs(AtabR2));
        if(~isnan(maxC))
            fid = find(abs(AtabR2)==maxC);
            if (AtabP(fid)<=0.05)
                LAtab(L,j) = AtabC(fid);
                LAlag(L,j) = yr(fid);
            end
        end
        clear maxC

    end % driver

end % LME

LFtab(:,5) = lid;
LPtab(:,5) = lid;
LDtab(:,5) = lid;
LAtab(:,5) = lid;

%%
cname = {'TP','TB','Det','ZmLoss','LME'};
Ftab = array2table(LFtab,"VariableNames",cname);
Ptab = array2table(LPtab,"VariableNames",cname);
Dtab = array2table(LDtab,"VariableNames",cname);
Atab = array2table(LAtab,"VariableNames",cname);

writetable(Ftab,[rpath,'LME_sig_SLR_maxR2lag_driver_nu_F.csv'],...
    'Delimiter',',');
writetable(Ptab,[rpath,'LME_sig_SLR_maxR2lag_driver_nu_P.csv'],...
    'Delimiter',',');
writetable(Dtab,[rpath,'LME_sig_SLR_maxR2lag_driver_nu_D.csv'],...
    'Delimiter',',');
writetable(Atab,[rpath,'LME_sig_SLR_maxR2lag_driver_nu_A.csv'],...
    'Delimiter',',');

save([rpath,'LME_sig_SLR_maxR2lag_driver_nu.mat'],'lid',...
    'LFtab','LFlag','LPtab','LPlag','LDtab','LDlag','LAtab','LAlag',...
    'tanom');


