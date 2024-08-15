% Calc corr of forcing-fish gam for LMEs

clear
close all

%% grid
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
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

cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% put into a matrix & use annual production
aanom(:,:,1) = aba;
aanom(:,:,2) = ana;
aanom(:,:,3) = aga;

fanom(:,:,1) = abf;
fanom(:,:,2) = anf;
fanom(:,:,3) = agf;

panom(:,:,1) = abp;
panom(:,:,2) = anp;
panom(:,:,3) = agp;

danom(:,:,1) = abd;
danom(:,:,2) = and;
danom(:,:,3) = agd;

tanom = {'Biom','Prod','Rec'};

%% Fish data
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_Catch_1948-2015_ann_mean_anoms.mat'])

%% % Corr of forcing ---------------------------------------------------------
yst = 1;
yen = 68;

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = ana(:,1);
lid = find(~isnan(AA));

%Lags
yr = 0:5;

% Drivers
tanom2=tanom';
tanom2(:,2)=tanom2(:,1);
tanom2(:,3)=tanom2(:,1);
tanom2(:,4)=tanom2(:,1);

[Ymat,Jmat] = meshgrid(yr,1:length(tanom));

FtabC = nan*ones(length(tanom),length(yr),length(lid));
FtabP = nan*ones(length(tanom),length(yr),length(lid));
PtabC = FtabC;
PtabP = FtabC;
DtabC = FtabC;
DtabP = FtabC;
AtabC = FtabC;
AtabP = FtabC;

%%
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = num2str(i);

    for j = 1:length(tanom)

        %input forcing
        driver = tanom{j};

        for k=1:length(yr) %Correlation at diff lags
            t = yr(k);

            %             LME  time   driver                
            ats = ((aanom(i,yst:yen-t,j))') ;
            fts = ((fanom(i,yst:yen-t,j))') ;
            pts = ((panom(i,yst:yen-t,j))') ;
            dts = ((danom(i,yst:yen-t,j))') ;

            %Fish
            [rp,pp] = corrcoef(fts , (af(i,yst+t:yen))');
            FtabC(j,k,L) = rp(1,2);
            FtabP(j,k,L) = pp(1,2);
            clear rp pp

            [rp,pp] = corrcoef(pts , (ap(i,yst+t:yen))');
            PtabC(j,k,L) = rp(1,2);
            PtabP(j,k,L) = pp(1,2);
            clear rp pp

            [rp,pp] = corrcoef(dts , (ad(i,yst+t:yen))');
            DtabC(j,k,L) = rp(1,2);
            DtabP(j,k,L) = pp(1,2);
            clear rp pp

            [rp,pp] = corrcoef(ats , (aall(i,yst+t:yen))');
            AtabC(j,k,L) = rp(1,2);
            AtabP(j,k,L) = pp(1,2);
            clear rp pp

        end % time lag

    end % driver

end %LME

%%
save([spath,'LMEs_corr_catch_drivers_feisty_0_5_lag.mat'],'tanom','lid','yr',...
    'FtabC','FtabP','PtabC','PtabP','DtabC','DtabP',...
    'AtabC','AtabP');

