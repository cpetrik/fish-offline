% Calc linear regression of cpue with sst, chl, forcing, biom, nu
% find most sig driver and lag

clear
close all

%% FOSI input forcing

cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adety','atb','atp','azlosy','azoo');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% FEISTY outputs
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

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
eyr = 1961:2010;
[yr,fid] = intersect(fyr,eyr);

% put into a matrix & use annual production
manom(:,:,1) = atp(:,fid);
manom(:,:,2) = atb(:,fid);
manom(:,:,3) = adety(:,fid);
manom(:,:,4) = azlosy(:,fid);
manom(:,:,5) = aba(:,fid);
manom(:,:,6) = ana(:,fid);

%% Divide anoms by 2SD
sdm = std(manom,0,2);
sdm2 = repmat(sdm,1,50,1);
manom = manom ./ (2*sdm2);

%Fish
af = af ./ (2*std(af,0,2));
ap = ap ./ (2*std(ap,0,2));
ad = ad ./ (2*std(ad,0,2));
aall = aall ./ (2*std(aall,0,2));

%% Drivers from satellite obs
load([fpath 'lme_satellite_sst_chl_ann_mean_anoms.mat'])

%% match years
[~,cid] = intersect(eyr,cyr);
[~,tid] = intersect(eyr,tyr);

sst = asst(:,1:length(tid)); 
chl = achl(:,1:length(cid));
manom(:,tid,7) = sst ./ (2*std(sst,0,2));
manom(:,cid,8) = chl ./ (2*std(chl,0,2));

tanom = {'TP','TB','Det','ZmLoss','Biom','Prod','SST','chl'};

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
tanom2(:,7)=tanom2(:,1);
tanom2(:,8)=tanom2(:,1);
tanom2(:,9)=tanom2(:,1);
tanom2(:,10)=tanom2(:,1);
tanom2(:,11)=tanom2(:,1);
tanom2(:,12)=tanom2(:,1);
tanom2(:,13)=tanom2(:,1);
tanom2(:,14)=tanom2(:,1);
tanom2(:,15)=tanom2(:,1);

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

        if j==7
            yst = tid(1);
            yen = tid(end);

        elseif j==8
            yst = cid(1);
            yen = cid(end);
        else
            yst = 1;
            yen = length(eyr);
        end

        for k=1:length(yr) %Lin reg at diff lags
            t = yr(k);

            %               LME time    driver         
            sclim = ((manom(i,yst:yen-t,j))') ;

            %Fish
            mdl0 = fitlm(sclim, (af(i,yst+t:yen))');
            FtabC(j,k) = mdl0.Coefficients.Estimate(2);
            FtabP(j,k) = mdl0.Coefficients.pValue(2);
            FtabR2(j,k) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(sclim, (ap(i,yst+t:yen))');
            PtabC(j,k) = mdl0.Coefficients.Estimate(2);
            PtabP(j,k) = mdl0.Coefficients.pValue(2);
            PtabR2(j,k) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(sclim, (ad(i,yst+t:yen))');
            DtabC(j,k) = mdl0.Coefficients.Estimate(2);
            DtabP(j,k) = mdl0.Coefficients.pValue(2);
            DtabR2(j,k) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(sclim, (aall(i,yst+t:yen))');
            AtabC(j,k) = mdl0.Coefficients.Estimate(2);
            AtabP(j,k) = mdl0.Coefficients.pValue(2);
            AtabR2(j,k) = mdl0.Rsquared.Ordinary;
            clear mdl0

        end % time lag

    end % driver
    
    %% find max(R2)
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

writetable(Atab1,[spath,'LMEs_SLR_cpue_sat_driver_feisty_norec_maxR2_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ftab1,[spath,'LMEs_SLR_cpue_sat_driver_feisty_norec_maxR2_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'LMEs_SLR_cpue_sat_driver_feisty_norec_maxR2_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'LMEs_SLR_cpue_sat_driver_feisty_norec_maxR2_D.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'LMEs_SLR_cpue_sat_driver_feisty_norec_maxR2s.mat'],...
    'LFtab','LPtab','LDtab','LAtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','lid');

writetable(Atab1,[dpath,'LMEs_SLR_cpue_sat_driver_feisty_norec_maxR2_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ftab1,[dpath,'LMEs_SLR_cpue_sat_driver_feisty_norec_maxR2_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[dpath,'LMEs_SLR_cpue_sat_driver_feisty_norec_maxR2_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[dpath,'LMEs_SLR_cpue_sat_driver_feisty_norec_maxR2_D.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([dpath,'LMEs_SLR_cpue_sat_driver_feisty_norec_maxR2s.mat'],...
    'LFtab','LPtab','LDtab','LAtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','lid');

