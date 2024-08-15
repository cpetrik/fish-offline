% Map coeffs of linear regression of FEISTY FOSI
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

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

load([fpath 'FEISTY_FOSI_',mod,'_ann_mean_anoms.mat'],'aa','ab','ad','af','ap','as','am','al') % Anoms with linear trend removed

%% Vectorize inputs
nt = length(yanom);
nid = length(ID);

atp   = nan*ones(nid,nt);
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
canom = repmat(manom(1,:),nid,1);
canom(:,:,2) = repmat(manom(2,:),nid,1);
canom(:,:,3) = repmat(manom(3,:),nid,1);
canom(:,:,4) = repmat(manom(4,:),nid,1);
canom(:,:,5) = repmat(manom(5,:),nid,1);
canom(:,:,6) = repmat(manom(6,:),nid,1);
canom(:,:,7) = repmat(manom(7,:),nid,1);
canom(:,:,8) = repmat(manom(8,:),nid,1);
canom(:,:,9) = repmat(manom(9,:),nid,1);
canom(:,:,10) = repmat(manom(10,:),nid,1);
canom(:,:,11) = repmat(manom(11,:),nid,1);

%% correlations by grid cell
yr = 0:4; %lags

rS = nan*ones(nid,11,length(yr));
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

tic
for j=1:11 %Diff climate indices
    for k=1:length(yr) %Linear regression at diff lags

        %%
        t = yr(k);

        for i=1:nid

            %Fish
            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (as(i,yst(j)+t:yen(j)))');
            rS(i,j,k) = mdl0.Coefficients.Estimate(2);
            pS(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (am(i,yst(j)+t:yen(j)))');
            rM(i,j,k) = mdl0.Coefficients.Estimate(2);
            pM(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (al(i,yst(j)+t:yen(j)))');
            rL(i,j,k) = mdl0.Coefficients.Estimate(2);
            pL(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (af(i,yst(j)+t:yen(j)))');
            rF(i,j,k) = mdl0.Coefficients.Estimate(2);
            pF(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (ap(i,yst(j)+t:yen(j)))');
            rP(i,j,k) = mdl0.Coefficients.Estimate(2);
            pP(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (ad(i,yst(j)+t:yen(j)))');
            rD(i,j,k) = mdl0.Coefficients.Estimate(2);
            pD(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (aa(i,yst(j)+t:yen(j)))');
            rV(i,j,k) = mdl0.Coefficients.Estimate(2);
            pV(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (ab(i,yst(j)+t:yen(j)))');
            rB(i,j,k) = mdl0.Coefficients.Estimate(2);
            pB(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            % Inputs / Forcing
            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (atp(i,yst(j)+t:yen(j)))');
            rTp(i,j,k) = mdl0.Coefficients.Estimate(2);
            pTp(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (atb(i,yst(j)+t:yen(j)))');
            rTb(i,j,k) = mdl0.Coefficients.Estimate(2);
            pTb(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (adet(i,yst(j)+t:yen(j)))');
            rDet(i,j,k) = mdl0.Coefficients.Estimate(2);
            pDet(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (azoo(i,yst(j)+t:yen(j)))');
            rZb(i,j,k) = mdl0.Coefficients.Estimate(2);
            pZb(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

            mdl0 = fitlm((canom(i,yst(j):yen(j)-t,j))' , (azlos(i,yst(j)+t:yen(j)))');
            rZl(i,j,k) = mdl0.Coefficients.Estimate(2);
            pZl(i,j,k) = mdl0.Coefficients.pValue(2);
            clear mdl0

        end
    end
end
toc %73.3647 hrs

save([fpath 'grid_linreg_climate_inputs_fish_FOSI_fished_',mod,'.mat'],...
    'rS','rM','rL','rF','rP','rD','rV','rB','rTp','rTb','rDet','rZb','rZl',...
    'pS','pM','pL','pF','pP','pD','pV','pB','pTp','pTb','pDet','pZb','pZl',...
    'tanom','yr');
