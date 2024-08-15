% Linear regression FEISTY inputs LME means with climate anoms
% CESM FOSI

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/';

apath = '/Users/cpetrik/Dropbox/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'Climate_anomalies_annual_means.mat']);

%% Map data
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%% FOSI input forcing
% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat']);
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat']);

%% Loop over all LMEs and all Climate
lid = [54,1:2,10,3,5:7];
lname = {'CHK','EBS','GAK','HI','CCE','GMX','SE','NE'};
fyr = 1948:2015;

colororder({'k','b'}); close all

% Try looking at LR with 0-5 yr lag
yr = 0:5;
np=0;
nb=0;
nd=0;
nm=0;
nl=0;

for i=1:length(lid)
    for j=1:length(canom)
        lme = lid(i);
        ltex = lname{i};
        ctex = canom{j};
        
        for k=1:6
            t = yr(k);
            %% linear regression at diff lags

mdlp = fitlm(manom(j,yst(j):yen(j)-t) , atp(lme,yst(j)+t:yen(j)));

            [rp,pp] = corrcoef(atp(lme,yst(j)+t:yen(j)),manom(j,yst(j):yen(j)-t));
            rP(i,j,k) = rp(1,2); pP(i,j,k) = pp(1,2);
            
            [rb,pb] = corrcoef(atb(lme,yst(j)+t:yen(j)),manom(j,yst(j):yen(j)-t));
            rB(i,j,k) = rb(1,2); pB(i,j,k) = pb(1,2);
            
            [rd,pd] = corrcoef(adet(lme,yst(j)+t:yen(j)),manom(j,yst(j):yen(j)-t));
            rD(i,j,k) = rd(1,2); pD(i,j,k) = pd(1,2);
            
            [rm,pm] = corrcoef(azoo(lme,yst(j)+t:yen(j)),manom(j,yst(j):yen(j)-t));
            rM(i,j,k) = rm(1,2); pM(i,j,k) = pm(1,2);
            
            [rl,pl] = corrcoef(azlos(lme,yst(j)+t:yen(j)),manom(j,yst(j):yen(j)-t));
            rL(i,j,k) = rl(1,2); pL(i,j,k) = pl(1,2);
            
            if (pp(1,2)<=0.05)
                np=np+1;
                sigP{np,1} = ltex;
                sigP{np,2} = ctex;
                sigP{np,3} = t;
                sigP{np,4} = rp(1,2);
                sigP{np,5} = pp(1,2);
            end
            if (pb(1,2)<=0.05)
                nb=nb+1;
                sigB{nb,1} = ltex;
                sigB{nb,2} = ctex;
                sigB{nb,3} = t;
                sigB{nb,4} = rb(1,2);
                sigB{nb,5} = pb(1,2);
            end
            if (pd(1,2)<=0.05)
                nd=nd+1;
                sigD{nd,1} = ltex;
                sigD{nd,2} = ctex;
                sigD{nd,3} = t;
                sigD{nd,4} = rd(1,2);
                sigD{nd,5} = pd(1,2);
            end
            if (pm(1,2)<=0.05)
                nm=nm+1;
                sigM{nm,1} = ltex;
                sigM{nm,2} = ctex;
                sigM{nm,3} = t;
                sigM{nm,4} = rm(1,2);
                sigM{nm,5} = pm(1,2);
            end
            if (pl(1,2)<=0.05)
                nl=nl+1;
                sigL{nl,1} = ltex;
                sigL{nl,2} = ctex;
                sigL{nl,3} = t;
                sigL{nl,4} = rl(1,2);
                sigL{nl,5} = pl(1,2);
            end
            
        end
    end
end

%%
% writecell(sigM,[cpath 'LME_fosi_inputs_sigLzooC.csv']);
% writecell(sigL,[cpath 'LME_fosi_inputs_sigZooLoss.csv']);
% writecell(sigP,[cpath 'LME_fosi_inputs_sigTP.csv']);
% writecell(sigD,[cpath 'LME_fosi_inputs_sigDet.csv']);
% writecell(sigB,[cpath 'LME_fosi_inputs_sigTB.csv']);








