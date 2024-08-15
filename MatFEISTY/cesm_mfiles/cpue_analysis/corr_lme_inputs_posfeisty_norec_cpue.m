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

tanom = {'TP','TB','Det','ZmLoss','Biom','Prod'};

manom = manom(:,fid,:);

fanom = manom;
fanom(:,:,5) = abf(:,fid);
fanom(:,:,6) = anf(:,fid);

panom = manom;
panom(:,:,5) = abp(:,fid);
panom(:,:,6) = anp(:,fid);

danom = manom;
danom(:,:,5) = abd(:,fid);
danom(:,:,6) = and(:,fid);

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

       
        for k=1:length(yr) %Linear regression at diff lags
            t = yr(k);

            %               LME  time   driver
            aclim = ((manom(L,yst:yen-t,j))');
            fclim = ((fanom(L,yst:yen-t,j))');
            pclim = ((panom(L,yst:yen-t,j))');
            dclim = ((danom(L,yst:yen-t,j))');

            %Fish
            [rp,pp] = corrcoef(fclim , (af(L,yst+t:yen))');
            FtabC(k) = rp(1,2);
            FtabP(k) = pp(1,2);
            clear rp pp

            [rp,pp] = corrcoef(pclim , (ap(L,yst+t:yen))');
            PtabC(k) = rp(1,2);
            PtabP(k) = pp(1,2);
            clear rp pp

            [rp,pp] = corrcoef(dclim , (ad(L,yst+t:yen))');
            DtabC(k) = rp(1,2);
            DtabP(k) = pp(1,2);
            clear rp pp

            [rp,pp] = corrcoef(aclim , (aall(L,yst+t:yen))');
            AtabC(k) = rp(1,2);
            AtabP(k) = pp(1,2);
            clear rp pp

        end % time lag

        %% Find lag with max corr for that driver
        % force prey & fish corrs to be pos or zero (3,4,5,6,7)
        clear maxC
        if j>2
            maxC = max((FtabC));
            maxC = max(0,maxC);
            fid = find((FtabC)==maxC);
            if(~isempty(fid))
                LFtab(L,j) = FtabC(fid);
                LFsig(L,j) = FtabP(fid);
                LFlag(L,j) = yr(fid);
            end
            clear maxC

            maxC = max((PtabC));
            maxC = max(0,maxC);
            fid = find((PtabC)==maxC);
            if(~isempty(fid))
                LPtab(L,j) = PtabC(fid);
                LPsig(L,j) = PtabP(fid);
                LPlag(L,j) = yr(fid);
            end
            clear maxC

            maxC = max((DtabC));
            maxC = max(0,maxC);
            fid = find((DtabC)==maxC);
            if(~isempty(fid))
                LDtab(L,j) = DtabC(fid);
                LDsig(L,j) = DtabP(fid);
                LDlag(L,j) = yr(fid);
            end
            clear maxC

            maxC = max((AtabC));
            maxC = max(0,maxC);
            fid = find((AtabC)==maxC);
            if(~isempty(fid))
                LAtab(L,j) = AtabC(fid);
                LAsig(L,j) = AtabP(fid);
                LAlag(L,j) = yr(fid);
            end
            clear maxC

        else
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
        end %if


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

writetable(Ftab,[rpath,'LME_corr_maxlag_driver_posfeisty_norec_cpue_F.csv'],...
    'Delimiter',',');
writetable(Ptab,[rpath,'LME_corr_maxlag_driver_posfeisty_norec_cpue_P.csv'],...
    'Delimiter',',');
writetable(Dtab,[rpath,'LME_corr_maxlag_driver_posfeisty_norec_cpue_D.csv'],...
    'Delimiter',',');
writetable(Atab,[rpath,'LME_corr_maxlag_driver_posfeisty_norec_cpue_A.csv'],...
    'Delimiter',',');

save([rpath,'LME_corr_maxlag_driver_posfeisty_norec_cpue.mat'],'lid',...
    'LFtab','LFlag','LPtab','LPlag','LDtab','LDlag','LAtab','LAlag',...
    'LFsig','LPsig','LDsig','LAsig','tanom','cname');


