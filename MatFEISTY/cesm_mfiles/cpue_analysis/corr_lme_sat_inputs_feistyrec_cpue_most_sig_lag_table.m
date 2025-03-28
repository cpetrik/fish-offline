% Calc corr of cpue with forcing, biomass, nu, recruit
% find most sig driver and lag

clear
close all

%% FOSI input forcing

%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adety','atb','atp','azlosy','azoo');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% FEISTY outputs
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
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

%Rec
load([fpath 'FEISTY_FOSI_',mod,'_lme_gam_rec_ann_mean_anoms.mat'],...
    'arf','arp','ard','ara');

%% Fish data
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_CPUE_1961-2010_ann_mean_anoms.mat'],'aall')

%% subset effort years
fyr = 1948:2015;
eyr = 1961:2010;
[yr,fid] = intersect(fyr,eyr);

adety  = adety(:,fid);
atb    = atb(:,fid);
atp    = atp(:,fid);
azlosy = azlosy(:,fid);
aba    = aba(:,fid);
ana    = ana(:,fid);
ara    = ara(:,fid);

% put into a matrix & use annual nuuction
manom(:,:,1) = atp;
manom(:,:,2) = atb;
manom(:,:,3) = adety;
manom(:,:,4) = azlosy;
manom(:,:,5) = aba;
manom(:,:,6) = ana;
manom(:,:,7) = ara;

%% Drivers from satellite obs
load([fpath 'lme_satellite_sst_chl_ann_mean_anoms.mat'])

%% match years
[~,cid] = intersect(eyr,cyr);
[~,tid] = intersect(eyr,tyr);

manom(:,tid,8) = asst(:,1:length(tid));
manom(:,cid,9) = achl(:,1:length(cid));

tanom = {'TP','TB','Det','ZmLoss','Biom','Prod','Rec','SST','chl'};

%% %Corr of forcing ---------------------------------------------------------
cnam = {'corr','p','lag','idriver','driver'};

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = aba(:,1);
lid = find(~isnan(AA));

%Lags
yr = 0:3;  %reduce lags 0:5

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

%%
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
            yen = length(eyr);
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
    %save([spath,ilme,'_corr_drivers_0_5_lag_nu.mat'])

    %%
    maxC = max(abs(AtabC(:)));
    pid = find(abs(AtabC(:))==maxC);
    LAtab(L,1) = AtabC(pid);
    LAtab(L,2) = AtabP(pid);
    LAtab(L,3) = Ymat(pid);
    LAtab(L,4) = Jmat(pid);
    LAt(L) = tanom2(pid);

    clear pid maxC

end %LME

%%
lname = cellstr(num2str(lid));
% cnam, lname, tanom2
Atab1 = array2table(LAtab,"RowNames",lname);
Atab1(:,5) = LAt;
Atab1.Properties.VariableNames = cnam;


%%
writetable(Atab1,[spath,'LMEs_corr_cpue_sat_driver_feistyrec_maxcorr_A.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'LMEs_corr_cpue_sat_driver_feistyrec_maxcorrs.mat'],...
    'LAtab','Atab1','lid');

