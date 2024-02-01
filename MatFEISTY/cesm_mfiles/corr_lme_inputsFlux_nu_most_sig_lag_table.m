% Calc corr of forcing-fish 
% Nu instead of biomass
% find most sig driver and lag

clear
close all

%% FOSI input forcing

%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adety','atb','atp','azlosy','azoo');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

% put into a matrix & use annual nuuction
manom(:,:,1) = atp;
manom(:,:,2) = atb;
manom(:,:,3) = adety;
manom(:,:,5) = azlosy;

tanom = {'TP','TB','Det','ZmLoss'};

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',cfile,'/corrs'];

mod = 'v15_All_fish03';

% Anoms with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms.mat'],...
    'aa','ad','af','ap')

%% %Corr of forcing ---------------------------------------------------------
yst = 1;
yen = 68;

cnam = {'corr','p','lag','idriver','driver'};

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

LFtab = nan*ones(length(lid),4);
LPtab = nan*ones(length(lid),4);
LDtab = nan*ones(length(lid),4);
LAtab = nan*ones(length(lid),4);

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

%%
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = num2str(i);

    for j = 1:length(tanom)

        %input forcing
        driver = tanom{j};

        for k=1:length(yr) %Correlations at diff lags
            t = yr(k);

            %               LME     time      driver         
            sclim = ((manom(i,yst:yen-t,j))') ;

            %Fish
            [rp,pp] = corrcoef(sclim , (af(i,yst+t:yen))');
            FtabC(j,k) = rp(1,2);
            FtabP(j,k) = pp(1,2);
            clear rp pp

            [rp,pp] = corrcoef(sclim , (ap(i,yst+t:yen))');
            PtabC(j,k) = rp(1,2);
            PtabP(j,k) = pp(1,2);
            clear rp pp

            [rp,pp] = corrcoef(sclim , (ad(i,yst+t:yen))');
            DtabC(j,k) = rp(1,2);
            DtabP(j,k) = pp(1,2);
            clear rp pp

            [rp,pp] = corrcoef(sclim , (aa(i,yst+t:yen))');
            AtabC(j,k) = rp(1,2);
            AtabP(j,k) = pp(1,2);
            clear rp pp

           
        end % time lag

    end % driver
    %%
    %save([spath,ilme,'_corr_drivers_0_5_lag_nu.mat'])

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


%%
writetable(Ftab1,[spath,'LMEs_corr_driver_maxcorr_Fnu.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'LMEs_corr_driver_maxcorr_Pnu.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'LMEs_corr_driver_maxcorr_Dnu.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Atab1,[spath,'LMEs_corr_driver_maxcorr_Anu.csv'],...
    'Delimiter',',','WriteRowNames',true);


save([spath,'LMEs_nu_corr_driver_maxcorrs.mat'],...
    'LFtab','LPtab','LDtab','LAtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','lid');

