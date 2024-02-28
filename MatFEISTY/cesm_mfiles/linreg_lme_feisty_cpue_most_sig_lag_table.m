% Calc linear regression of cpue with forcing, biomass, nu
% find most sig driver and lag

clear
close all

%% FEISTY outputs
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',cfile,'/corrs/'];

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

%% %Corr of forcing ---------------------------------------------------------
cnam = {'corr','p','R2','lag','idriver','driver'};

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
tanom2(:,5)=tanom2(:,1);
tanom2(:,6)=tanom2(:,1);

[Ymat,Jmat] = meshgrid(yr,1:length(tanom));

LFtab = nan*ones(length(lid),5);
LPtab = nan*ones(length(lid),5);
LDtab = nan*ones(length(lid),5);
LAtab = nan*ones(length(lid),5);

LFt = cell(length(lid),1);
LPt = cell(length(lid),1);
LDt = cell(length(lid),1);
LAt = cell(length(lid),1);

FtabC = nan*ones(length(tanom),length(yr));
FtabP = nan*ones(length(tanom),length(yr));
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


        yst = 1;
        yen = length(cyr);


        for k=1:length(yr) %Correlations at diff lags
            t = yr(k);

            %             LME  time   driver                
            ats = ((aanom(i,yst:yen-t,j))') ;
            fts = ((fanom(i,yst:yen-t,j))') ;
            pts = ((panom(i,yst:yen-t,j))') ;
            dts = ((danom(i,yst:yen-t,j))') ;

            %Fish
            mdl0 = fitlm(fts , (af(i,yst+t:yen))');
            FtabC(j,k) = mdl0.Coefficients.Estimate(2);
            FtabP(j,k) = mdl0.Coefficients.pValue(2);
            FtabR2(j,k) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(pts , (ap(i,yst+t:yen))');
            PtabC(j,k) = mdl0.Coefficients.Estimate(2);
            PtabP(j,k) = mdl0.Coefficients.pValue(2);
            PtabR2(j,k) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(dts , (ad(i,yst+t:yen))');
            DtabC(j,k) = mdl0.Coefficients.Estimate(2);
            DtabP(j,k) = mdl0.Coefficients.pValue(2);
            DtabR2(j,k) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(ats , (aall(i,yst+t:yen))');
            AtabC(j,k) = mdl0.Coefficients.Estimate(2);
            AtabP(j,k) = mdl0.Coefficients.pValue(2);
            AtabR2(j,k) = mdl0.Rsquared.Ordinary;
            clear mdl0

        end % time lag

    end % driver

    %% find max coeff, max R2, or min p-val? - R2 just in case standarization didn't work right
    maxC = max(abs(AtabR2(:)));
    pid = find(abs(AtabR2(:))==maxC);
    LAtab(L,1) = AtabC(pid);
    LAtab(L,2) = AtabP(pid);
    LAtab(L,3) = AtabR2(pid);
    LAtab(L,4) = Ymat(pid);
    LAtab(L,5) = Jmat(pid);
    LAt(L) = tanom2(pid);
    clear pid maxC

    maxC = max(abs(FtabR2(:)));
    pid = find(abs(FtabR2(:))==maxC);
    LFtab(L,1) = FtabC(pid);
    LFtab(L,2) = FtabP(pid);
    LFtab(L,3) = FtabR2(pid);
    LFtab(L,4) = Ymat(pid);
    LFtab(L,5) = Jmat(pid);
    LFt(L) = tanom2(pid);
    clear pid maxC

    maxC = max(abs(PtabR2(:)));
    if(L~=61)
        pid = find(abs(PtabR2(:))==maxC);
        LPtab(L,1) = PtabC(pid);
        LPtab(L,2) = PtabP(pid);
        LPtab(L,3) = PtabR2(pid);
        LPtab(L,4) = Ymat(pid);
        LPtab(L,5) = Jmat(pid);
        LPt(L) = tanom2(pid);
    end
    clear pid maxC

    maxC = max(abs(DtabR2(:)));
    pid = find(abs(DtabR2(:))==maxC);
    LDtab(L,1) = DtabC(pid);
    LDtab(L,2) = DtabP(pid);
    LDtab(L,3) = DtabR2(pid);
    LDtab(L,4) = Ymat(pid);
    LDtab(L,5) = Jmat(pid);
    LDt(L) = tanom2(pid);
    clear pid maxC

end %LME

%%
lname = cellstr(num2str(lid));
% cnam, lname, tanom2
Atab1 = array2table(LAtab,"RowNames",lname);
Atab1(:,6) = LAt;
Atab1.Properties.VariableNames = cnam;

Ftab1 = array2table(LFtab,"RowNames",lname);
Ftab1(:,6) = LFt;
Ftab1.Properties.VariableNames = cnam;

Ptab1 = array2table(LPtab,"RowNames",lname);
Ptab1(:,6) = LPt;
Ptab1.Properties.VariableNames = cnam;

Dtab1 = array2table(LDtab,"RowNames",lname);
Dtab1(:,6) = LDt;
Dtab1.Properties.VariableNames = cnam;


%%
dpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';

writetable(Atab1,[spath,'LMEs_SLR_cpue_feisty_maxR2_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ftab1,[spath,'LMEs_SLR_cpue_feisty_maxR2_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'LMEs_SLR_cpue_feisty_maxR2_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'LMEs_SLR_cpue_feisty_maxR2_D.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'LMEs_SLR_cpue_feisty_maxR2s.mat'],...
    'LFtab','LPtab','LDtab','LAtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','lid');

writetable(Atab1,[dpath,'LMEs_SLR_cpue_feisty_maxR2_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ftab1,[dpath,'LMEs_SLR_cpue_feisty_maxR2_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[dpath,'LMEs_SLR_cpue_feisty_maxR2_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[dpath,'LMEs_SLR_cpue_feisty_maxR2_D.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([dpath,'LMEs_SLR_cpue_feisty_maxR2s.mat'],...
    'LFtab','LPtab','LDtab','LAtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','lid');

