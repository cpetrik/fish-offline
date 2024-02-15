% Calc corr of catch with forcing, biomass, nu, gamma
% find most sig driver and lag

clear
close all

%% FOSI input forcing
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adety','atb','atp','azlosy','azoo');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% FEISTY outputs
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',cfile,'/corrs'];

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

% Test
% fpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';
%'FEISTY_FOSI_v15_All_fish03_lme_nu_ann_mean_anoms';
%'LME_fosi_fished_v15_All_fish03_Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%% Drivers from satellite obs
load([fpath 'lme_satellite_sst_chl_ann_mean_anoms.mat'])

%% match years
fyr = 1948:2015;
[~,cid] = intersect(fyr,cyr);
[~,tid] = intersect(fyr,tyr);

%% put into a matrix & use annual
manom(:,:,1) = atp;
manom(:,:,2) = atb;
manom(:,:,3) = adety;
manom(:,:,4) = azlosy;
manom(:,:,5) = aba;
manom(:,:,6) = ana;
manom(:,:,7) = aga;
manom(:,tid,8) = asst(:,1:length(tid));
manom(:,cid,9) = achl(:,1:length(cid));

tanom = {'TP','TB','Det','ZmLoss','Biom','Prod','Rec','SST','chl'};

%% Fish data
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

% Test
%ypath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/Phase3/fishing/';

% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_Catch_1948-2015_ann_mean_anoms.mat'],'aall')

%% %Corr of forcing ---------------------------------------------------------
cnam = {'corr','p','lag','idriver','driver'};

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = aba(:,1);
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
tanom2(:,7)=tanom2(:,1);
tanom2(:,8)=tanom2(:,1);
tanom2(:,9)=tanom2(:,1);
tanom2(:,10)=tanom2(:,1);

[Ymat,Jmat] = meshgrid(yr,1:length(tanom));

LAtab = nan*ones(length(lid),4);

LAt = cell(length(lid),1);

AtabC = nan*ones(length(tanom),length(yr));
AtabP = nan*ones(length(tanom),length(yr));

%% Test each driver with individ length of time series - also want to test min ts length for all
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = num2str(i);

    for j = 1:length(tanom)

        %input forcing
        driver = tanom{j};

        if j==8
            yst = tid(1);
            yen = tid(end);

        elseif j==9
            yst = cid(1);
            yen = cid(end);
        else
            yst = 1;
            yen = 68;
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

        end % time lag

    end % driver
    
    %%
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
Atab1 = array2table(LAtab,"RowNames",lname);
Atab1(:,5) = LAt;
Atab1.Properties.VariableNames = cnam;

%%
writetable(Atab1,[spath,'LMEs_corr_catch_sat_driver_feisty_maxcorr_A.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'LMEs_corr_catch_sat_driver_feisty_maxcorrs.mat'],...
    'LAtab','Atab1','lid');

