% Find patterns in forcing-fish correlations

clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";

% LMEs
lid = [54,1:2,65,10,3,5:7]; %ADD 65 = Aleutian Islands
lname = {'CHK','EBS','GAK','AI','HI','CCE','GMX','SE','NE'};

%% FOSI input forcing

%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adet','atb','atp','azlos','azoo');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

% put into a matrix
manom(:,:,1) = atp;
manom(:,:,2) = atb;
manom(:,:,3) = adet;
manom(:,:,4) = azoo;
manom(:,:,5) = azlos;

tanom = {'TP','TB','Det','Zmeso','ZmLoss'};

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

% Anoms with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'_lme_ann_mean_anoms.mat'],...
    'aa','ab','ad','af','ap','as','am','al')

%% % LM of forcing ---------------------------------------------------------
cnam = {'lag','coef','p','lagMC','maxcoef','pMC'};
ename = {'S','M','L','F','P','D','A','B'};

%Loop over drivers and responses
yr = 0:5;
yst = 1;
yen = 68;
lags = repmat(yr,length(ename),1);
lags = lags';

%%
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = lname{L};

    for j = 1:length(tanom)

        %input forcing
        driver = tanom{j};
        LtabC = nan*ones(length(ename),length(yr));
        LtabP = nan*ones(length(ename),length(yr));

        for k=1:length(yr) %Linear regression at diff lags
            t = yr(k);

            %               LME     time      driver                LME     time      driver
            sclim = ((manom(i,yst:yen-t,j))') ./ (2*std((manom(i,yst:yen-t,j))));

            %Fish
            mdl0 = fitlm(sclim , (as(i,yst+t:yen))');
            LtabC(1,k) = mdl0.Coefficients.Estimate(2);
            LtabP(1,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (am(i,yst+t:yen))');
            LtabC(2,k) = mdl0.Coefficients.Estimate(2);
            LtabP(2,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (al(i,yst+t:yen))');
            LtabC(3,k) = mdl0.Coefficients.Estimate(2);
            LtabP(3,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (af(i,yst+t:yen))');
            LtabC(4,k) = mdl0.Coefficients.Estimate(2);
            LtabP(4,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (ap(i,yst+t:yen))');
            LtabC(5,k) = mdl0.Coefficients.Estimate(2);
            LtabP(5,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (ad(i,yst+t:yen))');
            LtabC(6,k) = mdl0.Coefficients.Estimate(2);
            LtabP(6,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (aa(i,yst+t:yen))');
            LtabC(7,k) = mdl0.Coefficients.Estimate(2);
            LtabP(7,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (ab(i,yst+t:yen))');
            LtabC(8,k) = mdl0.Coefficients.Estimate(2);
            LtabP(8,k) = mdl0.Coefficients.pValue(2);
            clear mdl0
        end % time lag

        %%
        LtabC = LtabC';
        LtabP = LtabP';

        Stab = nan*ones(length(ename),6);
        Stab(:,1) = zeros(length(ename),1);
        Stab(:,2) = LtabC(1,:)';
        Stab(:,3) = LtabP(1,:)';
        
        Ptab = nan*ones(length(ename),6);
        Ptab(:,1) = zeros(length(ename),1);
        Ptab(:,2) = LtabC(1,:)';
        Ptab(:,3) = LtabP(1,:)';

        maxC = max(abs(LtabC));
        minP = min(LtabP);
        cid = nan*ones(8,1);
        pid = nan*ones(8,1);
        for z=1:8
            cid(z) = find(abs(LtabC(:,z))==maxC(z));
            pid(z) = find(LtabP(:,z)==minP(z));

            Stab(z,4) = yr(cid(z));
            Stab(z,5) = LtabC(cid(z),z);
            Stab(z,6) = LtabP(cid(z),z);

            Ptab(z,4) = yr(pid(z));
            Ptab(z,5) = LtabC(pid(z),z);
            Ptab(z,6) = LtabP(pid(z),z);
        end

        %%
        Atab1 = array2table(Stab,"VariableNames",cnam,...
            "RowNames",ename);
        Atab2 = array2table(Ptab,"VariableNames",cnam,...
            "RowNames",ename);

        writetable(Atab1,[spath,ilme,'_regress_',driver,...
            '_div2SD_maxcoef_lag.csv'],'Delimiter',',','WriteRowNames',true);
        
        writetable(Atab2,[spath,ilme,'_regress_',driver,...
            '_div2SD_minp_lag.csv'],'Delimiter',',','WriteRowNames',true);

        clear Atab1 Atab2

    end % LME

end %driver
