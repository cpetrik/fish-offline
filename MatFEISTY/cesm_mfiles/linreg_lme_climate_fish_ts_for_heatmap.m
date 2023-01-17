% Find patterns in forcing-fish correlations

clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";

%% % ------------------------------------------------------------
% load data
apath = '/Users/cpetrik/Dropbox/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'Climate_anomalies_annual_means.mat']);

tanom = canom;
clear canom

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

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

% Anoms with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'_lme_ann_mean_anoms.mat'],...
    'aa','ab','ad','af','ap','as','am','al')

%% LM of forcing ---------------------------------------------------------
cnam = {'Type','Lag','coef','p'};

%Loop over climate indices and responses
yr = 0:5;

for j = 1:11
    %climate
    driver = tanom{j};

    for L = 1:length(lid)
        %LME
        i = lid(L);
        ilme = lname{L};
        Ctab = cell(78,4);
    
        n = 0;

        for k=1:length(yr) %Linear regression at diff lags
            t = yr(k);
            
            sclim = ((manom(j,yst(j):yen(j)-t))') ./ (2*std((manom(j,yst(j):yen(j)-t))));

            %Fish
            n = n+1;
            mdl0 = fitlm(sclim , (as(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'S';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0
            
            n = n+1;
            mdl0 = fitlm(sclim , (am(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'M';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (al(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'L';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (af(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'F';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (ap(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'P';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (ad(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'D';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (aa(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'A';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (ab(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'B';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            % Inputs / Forcing
            n = n+1;
            mdl0 = fitlm(sclim , (atp(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'Tp';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (atb(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'Tb';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (adet(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'Det';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (azoo(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'Zoo';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (azlos(i,yst(j)+t:yen(j)))');
            Ctab{n,1} = 'Zlos';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0
        end % Lag
 
        %%
        Atab = array2table(Ctab,"VariableNames",cnam);

        writetable(Atab,[fpath,ilme,'_regress_',driver,...
            '_div2SD_melt_mat.csv'],'Delimiter',',');

        clear Atab

    end % LME

end %Climate
