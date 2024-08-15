% Find patterns in forcing-fish correlations 
% biophys & effort
% Simulated prod (nu)
% Do not divide by 2SD

clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/corrs/";

%% FOSI input forcing

%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adety','atb','atp','azlosy');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% Fishing effort
epath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';
load([epath 'FishMIP_Phase3a_LME_effort_1948-2010_ann_mean_anoms.mat']);

fyr = 1948:2015;
cyr = 1948:2010;
[yr,fid] = intersect(fyr,cyr);

% put into a matrix & use annual production
manom(:,:,1) = atp(:,fid);
manom(:,:,2) = atb(:,fid);
manom(:,:,3) = adety(:,fid);
manom(:,:,4) = azlosy(:,fid);
manom(:,:,5) = aall;
manom(:,:,6) = af;
manom(:,:,7) = ap;
manom(:,:,8) = ad;

tanom = {'TP','TB','Det','ZmLoss','Eff'};

clear ad af ap

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
rpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];

mod = 'v15_obsfish';

% Anoms with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'_lme_catch_ann_mean_anoms.mat'],...
    'af','ap','ad','aa');

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
% AA = aa(:,1);
% lid = find(~isnan(AA));
lid=1:66;

%% % LM of forcing ---------------------------------------------------------
%Loop over drivers and responses
yr = 0:5;
yst = 1;
yen = length(fid);

LFtab = nan*ones(length(lid),length(tanom));
LPtab = nan*ones(length(lid),length(tanom));
LDtab = nan*ones(length(lid),length(tanom));
LAtab = nan*ones(length(lid),length(tanom));

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

        for k=1:length(yr) %Pearson corr at diff lags
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
            [rp,pp] = corrcoef(sclim , (aa(L,yst+t:yen))');
            AtabC(k) = rp(1,2);
            AtabP(k) = pp(1,2);
            clear rp pp

        end % time lag

        %% Find lag with max corr for that driver
        maxC = max(abs(FtabC));
        if(~isnan(maxC))
            fid = find(abs(FtabC)==maxC);
            if (FtabP(fid)<=0.05)
                LFtab(L,j) = FtabC(fid);
                LFlag(L,j) = yr(fid);
            end
        end
        clear maxC

        maxC = max(abs(PtabC));
        if(~isnan(maxC))
            fid = find(abs(PtabC)==maxC);
            if (PtabP(fid)<=0.05)
                LPtab(L,j) = PtabC(fid);
                LPlag(L,j) = yr(fid);
            end
        end
        clear maxC

        maxC = max(abs(DtabC));
        if(~isnan(maxC))
            fid = find(abs(DtabC)==maxC);
            if (DtabP(fid)<=0.05)
                LDtab(L,j) = DtabC(fid);
                LDlag(L,j) = yr(fid);
            end
        end
        clear maxC

        maxC = max(abs(AtabC));
        if(~isnan(maxC))
            fid = find(abs(AtabC)==maxC);
            if (AtabP(fid)<=0.05)
                LAtab(L,j) = AtabC(fid);
                LAlag(L,j) = yr(fid);
            end
        end
        clear maxC

    end % driver

end % LME

LFtab(:,length(tanom)+1) = lid;
LPtab(:,length(tanom)+1) = lid;
LDtab(:,length(tanom)+1) = lid;
LAtab(:,length(tanom)+1) = lid;

%%
cname = {'TP','TB','Det','ZmLoss','Eff','LME'};
Ftab = array2table(LFtab,"VariableNames",cname);
Ptab = array2table(LPtab,"VariableNames",cname);
Dtab = array2table(LDtab,"VariableNames",cname);
Atab = array2table(LAtab,"VariableNames",cname);

writetable(Ftab,[rpath,'LME_sig_corr_maxlag_driver_effort_nu_obsfish_F.csv'],...
    'Delimiter',',');
writetable(Ptab,[rpath,'LME_sig_corr_maxlag_driver_effort_nu_obsfish_P.csv'],...
    'Delimiter',',');
writetable(Dtab,[rpath,'LME_sig_corr_maxlag_driver_effort_nu_obsfish_D.csv'],...
    'Delimiter',',');
writetable(Atab,[rpath,'LME_sig_corr_maxlag_driver_effort_nu_obsfish_A.csv'],...
    'Delimiter',',');

save([rpath,'LME_sig_corr_maxlag_driver_effort_nu_obsfish.mat'],'lid',...
    'LFtab','LFlag','LPtab','LPlag','LDtab','LDlag','LAtab','LAlag',...
    'tanom');


