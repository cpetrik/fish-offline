% Calc linear regression of biomass and nu against CPUE for LMEs

clear
close all

%% grid
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% FEISTY outputs
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',cfile,'/corrs'];

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

%% Fish data
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_CPUE_1961-2010_ann_mean_anoms.mat'])

%% subset effort years
fyr = 1948:2015;
cyr = 1961:2010;
[yr,fid] = intersect(fyr,cyr);

aba    = aba(:,fid);
ana    = ana(:,fid);
abf    = abf(:,fid);
anf    = anf(:,fid);
abp    = abp(:,fid);
anp    = anp(:,fid);
abd    = abd(:,fid);
and    = and(:,fid);

%% Divide anoms by 2SD
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

%% put into a matrix & use annual nuuction
aanom(:,:,1) = aba;
aanom(:,:,2) = ana;

fanom(:,:,1) = abf;
fanom(:,:,2) = anf;

panom(:,:,1) = abp;
panom(:,:,2) = anp;

danom(:,:,1) = abd;
danom(:,:,2) = and;

tanom = {'Biom','Prod'};

%% Linear regression ---------------------------------------------------------
yst = 1;
yen = length(cyr);

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = aba(:,1);
lid = find(~isnan(AA));

%Lags
yr = 0:5;

% Drivers
tanom2=tanom';
tanom2(:,2)=tanom2(:,1);
tanom2(:,3)=tanom2(:,1);
tanom2(:,4)=tanom2(:,1);

[Ymat,Jmat] = meshgrid(yr,1:length(tanom));

FtabC = nan*ones(length(tanom),length(yr),length(lid));
FtabP = nan*ones(length(tanom),length(yr),length(lid));
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

            %             LME  time   driver                
            ats = ((aanom(i,yst:yen-t,j))') ;
            fts = ((fanom(i,yst:yen-t,j))') ;
            pts = ((panom(i,yst:yen-t,j))') ;
            dts = ((danom(i,yst:yen-t,j))') ;

            %Fish
            mdl0 = fitlm(fts , (af(i,yst+t:yen))');
            FtabC(j,k,L) = mdl0.Coefficients.Estimate(2);
            FtabP(j,k,L) = mdl0.Coefficients.pValue(2);
            FtabR2(j,k,L) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(pts , (ap(i,yst+t:yen))');
            PtabC(j,k,L) = mdl0.Coefficients.Estimate(2);
            PtabP(j,k,L) = mdl0.Coefficients.pValue(2);
            PtabR2(j,k,L) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(dts , (ad(i,yst+t:yen))');
            DtabC(j,k,L) = mdl0.Coefficients.Estimate(2);
            DtabP(j,k,L) = mdl0.Coefficients.pValue(2);
            DtabR2(j,k,L) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(ats , (aall(i,yst+t:yen))');
            AtabC(j,k,L) = mdl0.Coefficients.Estimate(2);
            AtabP(j,k,L) = mdl0.Coefficients.pValue(2);
            AtabR2(j,k,L) = mdl0.Rsquared.Ordinary;
            clear mdl0

        end % time lag

    end % driver

end %LME

%%
save([spath,'LMEs_SLR_cpue_feisty_0_5_lag.mat'],'tanom','lid','yr',...
    'FtabC','FtabP','PtabC','PtabP','DtabC','DtabP',...
    'AtabC','AtabP','FtabR2','PtabR2','DtabR2','AtabR2');

