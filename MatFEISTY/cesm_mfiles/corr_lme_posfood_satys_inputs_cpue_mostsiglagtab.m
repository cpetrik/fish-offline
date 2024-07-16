% Calc corr of cpue with forcings
% find most sig driver and lag
% min yrs as sat chl

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

%% FEISTY outputs
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',cfile,'/corrs'];

mod = 'v15_All_fish03';

%% Fish data
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_CPUE_1961-2010_ann_mean_anoms.mat'])

%% subset effort years
fyr = 1948:2015;
eyr = 1961:2010;
[yr,fid] = intersect(fyr,eyr);

adety  = adety(:,fid);
atb    = atb(:,fid);
atp    = atp(:,fid);
azlosy = azlosy(:,fid);

% put into a matrix & use annual nuuction
manom(:,:,1) = atp;
manom(:,:,2) = atb;
manom(:,:,3) = adety;
manom(:,:,4) = azlosy;

%% Drivers from satellite obs
load([fpath 'lme_satellite_sst_chl_ann_mean_anoms.mat'])

%% match years
[~,cid] = intersect(eyr,cyr);
[~,tid] = intersect(eyr,tyr);

manom(:,tid,5) = asst(:,1:length(tid));
manom(:,cid,6) = achl(:,1:length(cid));

tanom = {'TP','TB','Det','ZmLoss','SST','chl'};

%% restrict analysis to only years with satellite chl data
manom = manom(:,cid,:);
aall = aall(:,cid);
af = af(:,cid);
ap = ap(:,cid);
ad = ad(:,cid);

%% %Corr of forcing ---------------------------------------------------------
cnam = {'corr','p','lag','idriver','driver'};

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = azlosy(:,1);
lid = find(~isnan(AA));

%Lags
yr = 0:4;  %reduce lags 0:4

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
yst = 1;
yen = length(cid);

for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = num2str(i);

    for j = 1:length(tanom)

        %input forcing
        driver = tanom{j};
        %Lags based on driver
        if(j>4 && j<7)
            yr = 0:3;
        else
            yr = 0:4;  %reduce lags 0:4
        end

        for k=1:length(yr) %Correlations at diff lags
            t = yr(k);

            %               LME     time      driver
            sclim = ((manom(i,yst:yen-t,j))') ;

            %Fish
            [rp,pp] = corrcoef(sclim , (aall(i,yst+t:yen))');
            AtabC(j,k) = rp(1,2);
            AtabP(j,k) = pp(1,2);
            clear rp pp

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

        end % time lag

    end % driver
    
    %% force prey corrs to be pos or zero (3,4,6)
    AtabC(3,AtabC(3,:)<0) = 0;
    AtabC(4,AtabC(4,:)<0) = 0;
    AtabC(6,AtabC(6,:)<0) = 0;
    
    FtabC(3,FtabC(3,:)<0) = 0;
    FtabC(4,FtabC(4,:)<0) = 0;
    FtabC(6,FtabC(6,:)<0) = 0;
    
    PtabC(3,PtabC(3,:)<0) = 0;
    PtabC(4,PtabC(4,:)<0) = 0;
    PtabC(6,PtabC(6,:)<0) = 0;
    
    DtabC(3,DtabC(3,:)<0) = 0;
    DtabC(4,DtabC(4,:)<0) = 0;
    DtabC(6,DtabC(6,:)<0) = 0;
   

    %%
    maxC = max(abs(AtabC(:)));
    pid = find(abs(AtabC(:))==maxC);
    LAtab(L,1) = AtabC(pid);
    LAtab(L,2) = AtabP(pid);
    LAtab(L,3) = Ymat(pid);
    LAtab(L,4) = Jmat(pid);
    LAt(L) = tanom2(pid);
    clear pid maxC

    maxC = max(abs(FtabC(:)));
    pid = find(abs(FtabC(:))==maxC);
    LFtab(L,1) = FtabC(pid);
    LFtab(L,2) = FtabP(pid);
    LFtab(L,3) = Ymat(pid);
    LFtab(L,4) = Jmat(pid);
    LFt(L) = tanom2(pid);
    clear pid maxC

    maxC = max(abs(PtabC(:)));
    if(~isnan(maxC))
        pid = find(abs(PtabC(:))==maxC);
        LPtab(L,1) = PtabC(pid);
        LPtab(L,2) = PtabP(pid);
        LPtab(L,3) = Ymat(pid);
        LPtab(L,4) = Jmat(pid);
        LPt(L) = tanom2(pid);
    end
    clear pid maxC

    maxC = max(abs(DtabC(:)));
    pid = find(abs(DtabC(:))==maxC);
    LDtab(L,1) = DtabC(pid);
    LDtab(L,2) = DtabP(pid);
    LDtab(L,3) = Ymat(pid);
    LDtab(L,4) = Jmat(pid);
    LDt(L) = tanom2(pid);
    clear pid maxC

end %LME

%%
lname = cellstr(num2str(lid));
% cnam, lname, tanom2
Atab1 = array2table(LAtab,"RowNames",lname);
Atab1(:,5) = LAt;
Atab1.Properties.VariableNames = cnam;

Ftab1 = array2table(LFtab,"RowNames",lname);
Ftab1(:,5) = LFt;
Ftab1.Properties.VariableNames = cnam;

Ptab1 = array2table(LPtab,"RowNames",lname);
Ptab1(:,5) = LPt;
Ptab1.Properties.VariableNames = cnam;

Dtab1 = array2table(LDtab,"RowNames",lname);
Dtab1(:,5) = LDt;
Dtab1.Properties.VariableNames = cnam;


%%
writetable(Atab1,[spath,'LMEs_corr_cpue_satyrs_driver_maxcorr_posfood_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ftab1,[spath,'LMEs_corr_cpue_satyrs_driver_maxcorr_posfood_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'LMEs_corr_cpue_satyrs_driver_maxcorr_posfood_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'LMEs_corr_cpue_satyrs_driver_maxcorr_posfood_D.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'LMEs_corr_cpue_satyrs_driver_maxcorr_posfoods.mat'],...
    'LFtab','LPtab','LDtab','LAtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','lid');

