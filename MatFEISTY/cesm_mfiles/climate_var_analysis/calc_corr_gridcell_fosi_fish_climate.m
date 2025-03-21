% Map correlations of FEISTY FOSI
% w/climate indices

clear all
close all

%% Climate indices
apath = '/Users/cpetrik/Dropbox/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'Climate_anomalies_annual_means.mat']);

tanom = canom;
clear canom

%% FOSI input forcing
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat'],'Adet','Atb','Atp','Azlos','Azoo');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
ppath = [pp cfile '/corrs/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

sims = {'v15_All_fish03_';'v15_climatol_';'v15_varFood_';'v15_varTemp_'};
mod = sims{1};

load([fpath 'FEISTY_FOSI_',mod,'ann_mean_anoms.mat'],'aa','ab','ad','af','ap','as','am','al') % Anoms with linear trend removed

%% Vectorize inputs
nt = length(yanom);

atp   = nan*ones(length(ID),nt);
atb   = atp;
adet  = atp;
azoo  = atp;
azlos = atp;

for i=1:nt
    temp = Atp(:,:,i);
    atp(:,i) = temp(ID);
    clear temp

    temp = Atb(:,:,i);
    atb(:,i) = temp(ID);
    clear temp

    temp = Adet(:,:,i);
    adet(:,i) = temp(ID);
    clear temp

    temp = Azoo(:,:,i);
    azoo(:,i) = temp(ID);
    clear temp

    temp = Azlos(:,:,i);
    azlos(:,i) = temp(ID);
    clear temp
end

%%
clear Adet Atb Atp Azlos Azoo

%% rep climate anoms
canom = repmat(manom(1,:),length(ID),1);
canom(:,:,2) = repmat(manom(2,:),length(ID),1);
canom(:,:,3) = repmat(manom(3,:),length(ID),1);
canom(:,:,4) = repmat(manom(4,:),length(ID),1);
canom(:,:,5) = repmat(manom(5,:),length(ID),1);
canom(:,:,6) = repmat(manom(6,:),length(ID),1);
canom(:,:,7) = repmat(manom(7,:),length(ID),1);
canom(:,:,8) = repmat(manom(8,:),length(ID),1);
canom(:,:,9) = repmat(manom(9,:),length(ID),1);
canom(:,:,10) = repmat(manom(10,:),length(ID),1);
canom(:,:,11) = repmat(manom(11,:),length(ID),1);

%% correlations by grid cell
yr = 0:4; %lags

rS = nan*ones(length(ID),11,length(yr));
rM = rS;
rL = rS;
rF = rS;
rP = rS;
rD = rS;
rV = rS;
rB = rS;
pS = rS;
pM = pS;
pL = pS;
pF = pS;
pP = pS;
pD = pS;
pV = pS;
pB = pS;

rTp  = rS;
rTb  = rS;
rDet = rS;
rZb  = rS;
rZl  = rS;
pTp  = pS;
pTb  = pS;
pDet = pS;
pZb  = pS;
pZl  = pS;

for j=1:11
    for k=1:length(yr)
        %% Corr at diff lags
        t = yr(k);
        
        %Fish
        [rho,pval] = (corr((as(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((as(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((as(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rS(:,j,k) = [rS1;rS2;rS3];
        pS(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3
        
        [rho,pval] = (corr((am(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((am(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((am(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rM(:,j,k) = [rS1;rS2;rS3];
        pM(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3
        
        [rho,pval] = (corr((al(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((al(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((al(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rL(:,j,k) = [rS1;rS2;rS3];
        pL(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3
        
        [rho,pval] = (corr((af(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((af(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((af(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rF(:,j,k) = [rS1;rS2;rS3];
        pF(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3
        
        [rho,pval] = (corr((ap(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((ap(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((ap(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rP(:,j,k) = [rS1;rS2;rS3];
        pP(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3
        
        [rho,pval] = (corr((ad(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((ad(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((ad(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rD(:,j,k) = [rS1;rS2;rS3];
        pD(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3
        
        [rho,pval] = (corr((aa(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((aa(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((aa(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rV(:,j,k) = [rS1;rS2;rS3];
        pV(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3
        
        [rho,pval] = (corr((ab(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((ab(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((ab(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rB(:,j,k) = [rS1;rS2;rS3];
        pB(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3

        % Inputs / Forcing
        [rho,pval] = (corr((atp(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((atp(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((atp(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rTp(:,j,k) = [rS1;rS2;rS3];
        pTp(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3
        
        [rho,pval] = (corr((atb(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((atb(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((atb(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rTb(:,j,k) = [rS1;rS2;rS3];
        pTb(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3
        
        [rho,pval] = (corr((adet(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((adet(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((adet(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rDet(:,j,k) = [rS1;rS2;rS3];
        pDet(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3
        
        [rho,pval] = (corr((azoo(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((azoo(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((azoo(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rZb(:,j,k) = [rS1;rS2;rS3];
        pZb(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3
        
        [rho,pval] = (corr((azlos(1:28604,yst(j)+t:yen(j)))',(canom(1:28604,yst(j):yen(j)-t,1))'));
        rS1 = diag(rho);
        p1  = diag(pval);
        [rho,pval] = (corr((azlos(28605:57208,yst(j)+t:yen(j)))',(canom(28605:57208,yst(j):yen(j)-t,1))'));
        rS2 = diag(rho);
        p2  = diag(pval);
        [rho,pval] = (corr((azlos(57209:end,yst(j)+t:yen(j)))',(canom(57209:end,yst(j):yen(j)-t,1))'));
        rS3 = diag(rho);
        p3  = diag(pval);
        rZl(:,j,k) = [rS1;rS2;rS3];
        pZl(:,j,k) = [p1;p2;p3];
        clear rS1 rS2 rS3 p1 p2 p3
        
    end
end

%stopped on j=3, k=1
save([fpath 'grid_corrs_climate_inputs_fish_FOSI_fished_',mod,'.mat'],...
    'rS','rM','rL','rF','rP','rD','rV','rB','rTp','rTb','rDet','rZb','rZl',...
    'pS','pM','pL','pF','pP','pD','pV','pB','pTp','pTb','pDet','pZb','pZl',...
    'tanom','yr');

