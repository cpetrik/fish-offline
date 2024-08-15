% Find patterns in forcing-fish correlations
% Nu instead of biomass

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
    'adet','adety','atb','atp','azlos','azlosy','azoo');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

% put into a matrix & use annual nuuction
manom(:,:,1) = atp;
manom(:,:,2) = atb;
manom(:,:,3) = adety;
manom(:,:,4) = azoo;
manom(:,:,5) = azlosy;

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

% Anoms with linear trend removed (per yr instead of day)
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms.mat'],...
    'aa','ad','af','ap')

%% % LM of forcing ---------------------------------------------------------
% row = Type, col = Lag, dim 3 = driver
tanom = {'TP','TB','Det','Zmeso','ZmLoss'};
tfish = {'D','A','P','F'};
tyr = {'0','1','2','3','4','5'};

%Loop over drivers and responses
yr = 0:5;
yst = 1;
yen = 68;

save([fpath,'FOSI_nu_regress_drivers_div2SD.mat'],'tanom','tfish',...
    'tyr','yr');

%%
for L = 1:length(lid)
    %LME
    i = lid(L);
    ilme = lname{L};

    Cmat = nan*ones(4,length(yr),5);
    Pmat = nan*ones(4,length(yr),5);

    for j = 1:5 %forcing vars
        driver = tanom{j};

        for k=1:length(yr) %Linear regression at diff lags
            t = yr(k);

            %               LME time    driver                LME time     driver
            sclim = ((manom(i,yst:yen-t,j))') ./ (2*std((manom(i,yst:yen-t,j))));

            %Fish
             %Fish
            mdl0 = fitlm(sclim , (af(i,yst+t:yen))');
            Cmat(4,k,j) = mdl0.Coefficients.Estimate(2);
            Pmat(4,k,j) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (ap(i,yst+t:yen))');
            Cmat(3,k,j) = mdl0.Coefficients.Estimate(2);
            Pmat(3,k,j) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (ad(i,yst+t:yen))');
            Cmat(1,k,j) = mdl0.Coefficients.Estimate(2);
            Pmat(1,k,j) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm(sclim , (aa(i,yst+t:yen))');
            Cmat(2,k,j) = mdl0.Coefficients.Estimate(2);
            Pmat(2,k,j) = mdl0.Coefficients.pValue(2);
            clear mdl0

        end % time lag

    end %driver

    %%
    eval([ 'Cmat_nosize_' ilme '= Cmat;']); 
    eval([ 'Pmat_nosize_' ilme '= Pmat;']); 

    eval(['save(''' [fpath 'FOSI_nu_regress_drivers_div2SD.mat'] ''',''' ['Cmat_nosize_' ilme] ''',''-append'')']);
    eval(['save(''' [fpath 'FOSI_nu_regress_drivers_div2SD.mat'] ''',''' ['Pmat_nosize_' ilme] ''',''-append'')']);
    
end % LME'




