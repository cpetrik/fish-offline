% Calc linear regression of cpue with forcing
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
load([fpath 'FEISTY_FOSI_',mod,'_lme_ann_mean_anoms.mat'],...
    'af','ap','ad','aa','ab');

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
ab = ab ./ (2*std(ab,0,2));

%% %Corr of forcing ---------------------------------------------------------
cnam = {'corr','p','R2','lag','idriver','driver'};

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
tanom2(:,7)=tanom2(:,1);
tanom2(:,8)=tanom2(:,1);

[Ymat,Jmat] = meshgrid(yr,1:length(tanom));

LFtab = nan*ones(length(lid),5);
LPtab = nan*ones(length(lid),5);
LDtab = nan*ones(length(lid),5);
LAtab = nan*ones(length(lid),5);
LBtab = nan*ones(length(lid),5);

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
FtabR2 = FtabC;
PtabR2 = FtabC;
DtabR2 = FtabC;
AtabR2 = FtabC;
BtabR2 = FtabC;

%%
yst = 1;
yen = 68;
        
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = num2str(i);

    for j = 1:length(tanom)

        %input forcing
        driver = tanom{j};
        
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

            mdl0 = fitlm(sclim, (aa(i,yst+t:yen))');
            AtabC(j,k) = mdl0.Coefficients.Estimate(2);
            AtabP(j,k) = mdl0.Coefficients.pValue(2);
            AtabR2(j,k) = mdl0.Rsquared.Ordinary;
            clear mdl0

            mdl0 = fitlm(sclim, (ab(i,yst+t:yen))');
            BtabC(j,k) = mdl0.Coefficients.Estimate(2);
            BtabP(j,k) = mdl0.Coefficients.pValue(2);
            BtabR2(j,k) = mdl0.Rsquared.Ordinary;
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
    %if(i~=64)
        pid = find(abs(PtabR2(:))==maxC);
        LPtab(L,1) = PtabC(pid);
        LPtab(L,2) = PtabP(pid);
        LPtab(L,3) = PtabR2(pid);
        LPtab(L,4) = Ymat(pid);
        LPtab(L,5) = Jmat(pid);
        LPt(L) = tanom2(pid);
    %end
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

    maxC = max(abs(BtabR2(:)));
    pid = find(abs(BtabR2(:))==maxC);
    LBtab(L,1) = BtabC(pid);
    LBtab(L,2) = BtabP(pid);
    LBtab(L,3) = BtabR2(pid);
    LBtab(L,4) = Ymat(pid);
    LBtab(L,5) = Jmat(pid);
    LBt(L) = tanom2(pid);
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

Btab1 = array2table(LBtab,"RowNames",lname);
Btab1(:,6) = LBt;
Btab1.Properties.VariableNames = cnam;

%%
dpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';

writetable(Atab1,[spath,'LMEs_SLR_biom_driver_maxR2_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ftab1,[spath,'LMEs_SLR_biom_driver_maxR2_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'LMEs_SLR_biom_driver_maxR2_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'LMEs_SLR_biom_driver_maxR2_D.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Btab1,[spath,'LMEs_SLR_biom_driver_maxR2_B.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'LMEs_SLR_biom_driver_maxR2s.mat'],...
    'LFtab','LPtab','LDtab','LAtab','LBtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','Btab1','lid');

writetable(Atab1,[dpath,'LMEs_SLR_biom_driver_maxR2_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ftab1,[dpath,'LMEs_SLR_biom_driver_maxR2_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[dpath,'LMEs_SLR_biom_driver_maxR2_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[dpath,'LMEs_SLR_biom_driver_maxR2_D.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Btab1,[dpath,'LMEs_SLR_biom_driver_maxR2_B.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([dpath,'LMEs_SLR_biom_driver_maxR2s.mat'],...
    'LFtab','LPtab','LDtab','LAtab','LBtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','Btab1','lid');

