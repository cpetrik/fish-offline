% Forcing, fish, CPUE linear regression
% Divide by 2SD

clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";

%% FOSI input forcing
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adety','atb','atp','azlosy');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

% put into a matrix & use annual production
manom(:,:,1) = atp;
manom(:,:,2) = atb;
manom(:,:,3) = adety;
manom(:,:,4) = azlosy;

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
rpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];

mod = 'v15_All_fish03';

% Anoms with linear trend removed
%Biomass
load([fpath 'FEISTY_FOSI_',mod,'_lme_ann_mean_anoms.mat'],...
    'aa','ad','af','ap');

aba = aa;
abd = ad;
abf = af;
abp = ap;

clear aa ad af ap

%% Nu
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms.mat'],...
    'aa','ad','af','ap');

ana = aa;
and = ad;
anf = af;
anp = ap;

clear aa ad af ap

%% Fishing data
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_CPUE_1961-2010_ann_mean_anoms.mat'])

%% subset effort years
fyr = 1948:2015;
cyr = 1961:2010;
[yr,fid] = intersect(fyr,cyr);

manom = manom(:,fid,:);

aba    = aba(:,fid);
ana    = ana(:,fid);
abf    = abf(:,fid);
anf    = anf(:,fid);
abp    = abp(:,fid);
anp    = anp(:,fid);
abd    = abd(:,fid);
and    = and(:,fid);

% Divide anoms by 2SD
sdm = std(manom,0,2);
sdm2 = repmat(sdm,1,50,1);
manom = manom ./ (2*sdm2);

%Fish
abf = abf ./ (2*std(abf,0,2));
abp = abp ./ (2*std(abp,0,2));
abd = abd ./ (2*std(abd,0,2));
aba = aba ./ (2*std(aba,0,2));

anf = anf ./ (2*std(anf,0,2));
anp = anp ./ (2*std(anp,0,2));
and = and ./ (2*std(and,0,2));
ana = ana ./ (2*std(ana,0,2));

af = af ./ (2*std(af,0,2));
ap = ap ./ (2*std(ap,0,2));
ad = ad ./ (2*std(ad,0,2));
aall = aall ./ (2*std(aall,0,2));

%% put into a matrix & use annual production
fanom = manom;
fanom(:,:,5) = abf;
fanom(:,:,6) = anf;

panom = manom;
panom(:,:,5) = abp;
panom(:,:,6) = anp;

danom = manom;
danom(:,:,5) = abd;
danom(:,:,6) = and;

manom(:,:,5) = aba;
manom(:,:,6) = ana;

tanom = {'TP','TB','Det','ZmLoss','Biom','Prod'};

%% % LM of forcing ---------------------------------------------------------

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = aba(:,1);
lid = find(~isnan(AA));

%Loop over drivers and responses
yr = 0:5;
yst = 1;
yen = length(cyr);

LFtab = nan*ones(length(lid),length(tanom));
LPtab = nan*ones(length(lid),length(tanom));
LDtab = nan*ones(length(lid),length(tanom));
LAtab = nan*ones(length(lid),length(tanom));

LFsig = nan*ones(length(lid),length(tanom));
LPsig = nan*ones(length(lid),length(tanom));
LDsig = nan*ones(length(lid),length(tanom));
LAsig = nan*ones(length(lid),length(tanom));

LFlag = nan*ones(length(lid),length(tanom));
LPlag = nan*ones(length(lid),length(tanom));
LDlag = nan*ones(length(lid),length(tanom));
LAlag = nan*ones(length(lid),length(tanom));

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

for i = 1:length(lid) %LME

    L = lid(i);
    
    for j = 1:length(tanom) %drivers
        driver = tanom{j};

        for k=1:length(yr) %Linear regression at diff lags
            t = yr(k);

            %               LME  time   driver
            aclim = ((manom(L,yst:yen-t,j))');
            fclim = ((fanom(L,yst:yen-t,j))');
            pclim = ((panom(L,yst:yen-t,j))');
            dclim = ((danom(L,yst:yen-t,j))');

            %Fish
            mdl0 = fitlm(fclim , (af(L,yst+t:yen))');
            FtabC(k) = mdl0.Coefficients.Estimate(2);
            FtabP(k) = mdl0.Coefficients.pValue(2);
            FtabR2(k) = mdl0.Rsquared.Ordinary;
            clear mdl0

            if (L~=64)
                mdl0 = fitlm(pclim , (ap(L,yst+t:yen))');
                PtabC(k) = mdl0.Coefficients.Estimate(2);
                PtabP(k) = mdl0.Coefficients.pValue(2);
                PtabR2(k) = mdl0.Rsquared.Ordinary;
                clear mdl0
            end

            mdl0 = fitlm(dclim , (ad(L,yst+t:yen))');
            DtabC(k) = mdl0.Coefficients.Estimate(2);
            DtabP(k) = mdl0.Coefficients.pValue(2);
            DtabR2(k) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(aclim , (aall(L,yst+t:yen))');
            AtabC(k) = mdl0.Coefficients.Estimate(2);
            AtabP(k) = mdl0.Coefficients.pValue(2);
            AtabR2(k) = mdl0.Rsquared.Ordinary;
            clear mdl0

        end % time lag

        %% Find lag with max R2 for that driver
        maxC = max(abs(FtabR2));
        if(~isnan(maxC))
            fid = find(abs(FtabR2)==maxC);
            LFtab(i,j) = FtabC(fid);
            LFsig(i,j) = FtabP(fid);
            LFlag(i,j) = yr(fid);
        end
        clear maxC

        maxC = max(abs(PtabR2));
        if(~isnan(maxC))
            fid = find(abs(PtabR2)==maxC);
            LPtab(i,j) = PtabC(fid);
            LPsig(i,j) = PtabP(fid);
            LPlag(i,j) = yr(fid);
        end
        clear maxC

        maxC = max(abs(DtabR2));
        if(~isnan(maxC))
            fid = find(abs(DtabR2)==maxC);
            LDtab(i,j) = DtabC(fid);
            LDsig(i,j) = DtabP(fid);
            LDlag(i,j) = yr(fid);
        end
        clear maxC

        maxC = max(abs(AtabR2));
        if(~isnan(maxC))
            fid = find(abs(AtabR2)==maxC);
            LAtab(i,j) = AtabC(fid);
            LAsig(i,j) = AtabP(fid);
            LAlag(i,j) = yr(fid);
        end
        clear maxC

    end % driver

end % LME

LFtab(:,7) = lid;
LPtab(:,7) = lid;
LDtab(:,7) = lid;
LAtab(:,7) = lid;

LFsig(:,7) = lid;
LPsig(:,7) = lid;
LDsig(:,7) = lid;
LAsig(:,7) = lid;

LFlag(:,7) = lid;
LPlag(:,7) = lid;
LDlag(:,7) = lid;
LAlag(:,7) = lid;

%%
cname = {'TP','TB','Det','ZmLoss','Biom','Prod','LME'};
Ftab = array2table(LFtab,"VariableNames",cname);
Ptab = array2table(LPtab,"VariableNames",cname);
Dtab = array2table(LDtab,"VariableNames",cname);
Atab = array2table(LAtab,"VariableNames",cname);

writetable(Ftab,[rpath,'LME_linreg_maxR2lag_driver_cpue_F.csv'],...
    'Delimiter',',');
writetable(Ptab,[rpath,'LME_linreg_maxR2lag_driver_cpue_P.csv'],...
    'Delimiter',',');
writetable(Dtab,[rpath,'LME_linreg_maxR2lag_driver_cpue_D.csv'],...
    'Delimiter',',');
writetable(Atab,[rpath,'LME_linreg_maxR2lag_driver_cpue_A.csv'],...
    'Delimiter',',');

save([rpath,'LME_linreg_maxR2lag_driver_cpue.mat'],'lid',...
    'LFtab','LFlag','LPtab','LPlag','LDtab','LDlag','LAtab','LAlag',...
    'LFsig','LPsig','LDsig','LAsig','tanom','cname');


