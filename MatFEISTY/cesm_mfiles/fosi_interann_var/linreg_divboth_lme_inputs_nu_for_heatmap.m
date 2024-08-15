% Find patterns in forcing-fish correlations
% Nu instead of biomass
% divide both ts by 2SD


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

tanom = {'TP','TB','Det','Zmeso','ZmLoss'};

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
rpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

% Anoms with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms.mat'],...
    'aa','ad','af','ap')

%% Divide anoms by 2SD
sdm = std(manom,0,2);
sdm2 = repmat(sdm,1,68,1);
manom = manom ./ (2*sdm2);

%Fish
af = af ./ (2*std(af,0,2));
ap = ap ./ (2*std(ap,0,2));
ad = ad ./ (2*std(ad,0,2));
aa = aa ./ (2*std(aa,0,2));

%% % LM of forcing ---------------------------------------------------------
cnam = {'Type','Lag','coef','p'};

%Loop over drivers and responses
yr = 0:5;
yst = 1; 
yen = 68;

for j = 1:5
    %fish vars
    driver = tanom{j};

    for L = 1:length(lid)
        %LME
        i = lid(L);
        ilme = lname{L};
        Ctab = cell(48,4);

        n = 0;

        for k=1:length(yr) %Linear regression at diff lags
            t = yr(k);

            %               LME  time   driver         
            sclim = ((manom(i,yst:yen-t,j))') ;

            %Fish
            n = n+1;
            mdl0 = fitlm(sclim , (af(i,yst+t:yen))');
            Ctab{n,1} = 'F';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (ap(i,yst+t:yen))');
            Ctab{n,1} = 'P';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (ad(i,yst+t:yen))');
            Ctab{n,1} = 'D';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

            n = n+1;
            mdl0 = fitlm(sclim , (aa(i,yst+t:yen))');
            Ctab{n,1} = 'A';
            Ctab{n,2} = t;
            Ctab{n,3} = mdl0.Coefficients.Estimate(2);
            Ctab{n,4} = mdl0.Coefficients.pValue(2);
            clear mdl0

        end % time lag

        %%
        Atab = array2table(Ctab,"VariableNames",cnam);

        writetable(Atab,[rpath,ilme,'_regress_',driver,...
            '_ALLdiv2SD_melt_mat_nu.csv'],'Delimiter',',');

        clear Atab

    end % LME

end %driver
