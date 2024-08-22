% Calc corr of cpue with biomass, nu
% calc only once, then put together with others
% min yrs as sat chl 1982-2015

clear
close all

%% FOSI input forcing

%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% FEISTY outputs
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regress_cpue/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',...
    cfile,'/corrs_cpue'];

mod = 'v15_All_fish03';

% Anoms with linear trend removed
% Biomass
load([fpath 'FEISTY_FOSI_',mod,'_lme_biom_ann_mean_anoms_1982_2010_2015.mat'],'eyr',...
    'aba15','abd15','abf15','abp15');

% Nu
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms_1982_2010_2015.mat'],...
    'ana15','and15','anf15','anp15');

%% Fish data
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_CPUE_1982-2015_ann_mean_anoms.mat'],...
    'aa_cpue82','af_cpue82','ap_cpue82','ad_cpue82')

%% put into a matrix & use annual nuuction
aanom(:,:,1) = aba15;
aanom(:,:,2) = ana15;

fanom(:,:,1) = abf15;
fanom(:,:,2) = anf15;

panom(:,:,1) = abp15;
panom(:,:,2) = anp15;

danom(:,:,1) = abd15;
danom(:,:,2) = and15;

tanom = {'Biom','Prod'};

%% %Corr of forcing ---------------------------------------------------------
cnam = {'corr','p','lag','idriver','driver'};

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = aba15(:,1);
lid = find(~isnan(AA));

%Lags
yr = 0:2;  %reduce lags 

% Drivers
[Ymat,Jmat] = meshgrid(yr,1:length(tanom));

FtabC = nan*ones(length(lid),length(tanom),length(yr));
FtabP = nan*ones(length(lid),length(tanom),length(yr));
PtabC = FtabC;
PtabP = FtabC;
DtabC = FtabC;
DtabP = FtabC;
AtabC = FtabC;
AtabP = FtabC;

%%
yst = 1;
yen = length(eyr);

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
            %Fish
            [rp,pp] = corrcoef(((aanom(i,yst:yen-t,j))') , (aa_cpue82(i,yst+t:yen))');
            AtabC(L,j,k) = rp(1,2);
            AtabP(L,j,k) = pp(1,2);
            clear rp pp

            [rp,pp] = corrcoef(((fanom(i,yst:yen-t,j))') , (af_cpue82(i,yst+t:yen))');
            FtabC(L,j,k) = rp(1,2);
            FtabP(L,j,k) = pp(1,2);
            clear rp pp

            [rp,pp] = corrcoef(((panom(i,yst:yen-t,j))') , (ap_cpue82(i,yst+t:yen))');
            PtabC(L,j,k) = rp(1,2);
            PtabP(L,j,k) = pp(1,2);
            clear rp pp

            [rp,pp] = corrcoef(((danom(i,yst:yen-t,j))') , (ad_cpue82(i,yst+t:yen))');
            DtabC(L,j,k) = rp(1,2);
            DtabP(L,j,k) = pp(1,2);
            clear rp pp

        end % time lag
    end % driver
end %LME

%%
save([spath,'LMEs_corr_cpue_sstyrs15_feisty_lags.mat'],...
    'FtabC','PtabC','DtabC','AtabC',...
    'FtabP','PtabP','DtabP','AtabP','lid','cnam','tanom');



