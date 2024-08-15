% Find patterns in forcing-fish correlations
% SAUp CPUE
% Do not divide by 2SD

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

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
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

%Gamma
load([fpath 'FEISTY_FOSI_',mod,'_lme_gam_rec_ann_mean_anoms.mat'],...
    'agf','agp','agd','aga');

%% Fishing data
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_CPUE_1961-2010_ann_mean_anoms.mat'])

%%
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
% AA = aa(:,1);
% lid = find(~isnan(AA));
lid=1:66;

%%
fyr = 1948:2015;
cyr = 1961:2010;
[yr,fid] = intersect(fyr,cyr);

% put into a matrix & use annual nuuction
manom(:,:,1) = atp;
manom(:,:,2) = atb;
manom(:,:,3) = adety;
manom(:,:,4) = azlosy;
manom(:,:,5) = aba;
manom(:,:,6) = ana;
manom(:,:,7) = aga;

tanom = {'TP','TB','Det','ZmLoss','Biom','Prod','Rec'};

manom = manom(:,fid,:);

%% % LM of forcing ---------------------------------------------------------
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

for L = 1:length(lid) %LME
    
    for j = 1:length(tanom) %drivers
        driver = tanom{j};

        n = 0;

        for k=1:length(yr) %Linear regression at diff lags
            t = yr(k);

            %               LME  time   driver
            sclim = ((manom(L,yst:yen-t,j))');

            %Fish
            n = n+1;
            [rp,pp] = corrcoef(sclim , (af(L,yst+t:yen))');
            FtabC(k) = rp(1,2);
            FtabP(k) = pp(1,2);
            clear rp pp

            n = n+1;
            [rp,pp] = corrcoef(sclim , (ap(L,yst+t:yen))');
            PtabC(k) = rp(1,2);
            PtabP(k) = pp(1,2);
            clear rp pp

            n = n+1;
            [rp,pp] = corrcoef(sclim , (ad(L,yst+t:yen))');
            DtabC(k) = rp(1,2);
            DtabP(k) = pp(1,2);
            clear rp pp

            n = n+1;
            [rp,pp] = corrcoef(sclim , (aall(L,yst+t:yen))');
            AtabC(k) = rp(1,2);
            AtabP(k) = pp(1,2);
            clear rp pp

        end % time lag

        %% Find lag with max corr for that driver
        maxC = max(abs(FtabC));
        if(~isnan(maxC))
            fid = find(abs(FtabC)==maxC);
            LFtab(L,j) = FtabC(fid);
            LFsig(L,j) = FtabP(fid);
            LFlag(L,j) = yr(fid);
        end
        clear maxC

        maxC = max(abs(PtabC));
        if(~isnan(maxC))
            fid = find(abs(PtabC)==maxC);
            LPtab(L,j) = PtabC(fid);
            LPsig(L,j) = PtabP(fid);
            LPlag(L,j) = yr(fid);
        end
        clear maxC

        maxC = max(abs(DtabC));
        if(~isnan(maxC))
            fid = find(abs(DtabC)==maxC);
            LDtab(L,j) = DtabC(fid);
            LDsig(L,j) = DtabP(fid);
            LDlag(L,j) = yr(fid);
        end
        clear maxC

        maxC = max(abs(AtabC));
        if(~isnan(maxC))
            fid = find(abs(AtabC)==maxC);
            LAtab(L,j) = AtabC(fid);
            LAsig(L,j) = AtabP(fid);
            LAlag(L,j) = yr(fid);
        end
        clear maxC

    end % driver

end % LME

LFtab(:,8) = lid;
LPtab(:,8) = lid;
LDtab(:,8) = lid;
LAtab(:,8) = lid;

LFsig(:,8) = lid;
LPsig(:,8) = lid;
LDsig(:,8) = lid;
LAsig(:,8) = lid;

LFlag(:,8) = lid;
LPlag(:,8) = lid;
LDlag(:,8) = lid;
LAlag(:,8) = lid;

%%
cname = {'TP','TB','Det','ZmLoss','Biom','Prod','Rec','LME'};
Ftab = array2table(LFtab,"VariableNames",cname);
Ptab = array2table(LPtab,"VariableNames",cname);
Dtab = array2table(LDtab,"VariableNames",cname);
Atab = array2table(LAtab,"VariableNames",cname);

writetable(Ftab,[rpath,'LME_corr_maxlag_driver_cpue_F.csv'],...
    'Delimiter',',');
writetable(Ptab,[rpath,'LME_corr_maxlag_driver_cpue_P.csv'],...
    'Delimiter',',');
writetable(Dtab,[rpath,'LME_corr_maxlag_driver_cpue_D.csv'],...
    'Delimiter',',');
writetable(Atab,[rpath,'LME_corr_maxlag_driver_cpue_A.csv'],...
    'Delimiter',',');

save([rpath,'LME_corr_maxlag_driver_cpue.mat'],'lid',...
    'LFtab','LFlag','LPtab','LPlag','LDtab','LDlag','LAtab','LAlag',...
    'LFsig','LPsig','LDsig','LAsig','tanom','cname');


