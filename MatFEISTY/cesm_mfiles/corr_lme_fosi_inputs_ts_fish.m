% Corr LME biomass of FEISTY with input means
% CESM FOSI

clear all
close all

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
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
%load([cpath 'lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'])
spath='/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/FOSI/';
load([spath 'LME_fosi_input_anomalies_annual.mat']);

% put inputs in matrix
manom = nan*ones(5,66,68);
manom(1,:,:) = lme_tpa;
manom(2,:,:) = lme_tba;
manom(3,:,:) = lme_deta;
manom(4,:,:) = lme_mza;
manom(5,:,:) = lme_losa;

yanom = 1948:2015;

canom = {'Tp','Tb','Det','LZbiom','LZloss'};

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

dpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

sims = {'v14_All_fish03_';'v14_climatol_';'v14_varFood_';'v14_varTemp_'};
mod = sims{1};

%load([dpath 'LME_fosi_fished_',mod,cfile '.mat']);
fpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
load([fpath 'LME_fosi_fished_',mod,'anomalies_annual.mat'])


%% Loop over all LMEs and all inputs
lid = [54,1:2,10,3,5:7];
lname = {'CHK','EBS','GAK','HI','CCE','GMX','SE','NE'};
fyr = 1948:2015;
close all
colororder({'k','b'})

for i=1:length(lid)
    for j=1:length(canom)
        lme = lid(i);
        ltex = lname{i};
        ctex = canom{j};
        
        %% Cross corr 
        [cS,lagsS] = xcorr(lme_mSa(lme,:),squeeze(manom(j,lme,:)),15);
        [cM,lagsM] = xcorr(lme_mMa(lme,:),squeeze(manom(j,lme,:)),15);
        [cL,lagsL] = xcorr(lme_mLa(lme,:),squeeze(manom(j,lme,:)),15);
        [cF,lagsF] = xcorr(lme_mFa(lme,:),squeeze(manom(j,lme,:)),15);
        [cP,lagsP] = xcorr(lme_mPa(lme,:),squeeze(manom(j,lme,:)),15);
        [cD,lagsD] = xcorr(lme_mDa(lme,:),squeeze(manom(j,lme,:)),15);
        [cA,lagsA] = xcorr(lme_mAa(lme,:),squeeze(manom(j,lme,:)),15);
        [cB,lagsB] = xcorr(lme_mba(lme,:),squeeze(manom(j,lme,:)),15);
        
        %%
        figure(1)
        clf
        subplot(3,3,1)
        yyaxis left
        plot(yanom,squeeze(manom(j,lme,:)));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_mSa(lme,:));
        xlim([fyr(1) fyr(end)])
        title('Small')
        
        subplot(3,3,2)
        yyaxis left
        plot(yanom,squeeze(manom(j,lme,:)));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_mMa(lme,:));
        xlim([fyr(1) fyr(end)])
        str = {[ctex,' ', ltex], ' Medium'};
        title(str)
        
        subplot(3,3,3)
        yyaxis left
        plot(yanom,squeeze(manom(j,lme,:)));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_mLa(lme,:));
        xlim([fyr(1) fyr(end)])
        title('Large')
        
        subplot(3,3,4)
        yyaxis left
        plot(yanom,squeeze(manom(j,lme,:)));
        xlim([fyr(1) fyr(end)])
        ylabel(ctex)
        yyaxis right
        plot(fyr,lme_mFa(lme,:));
        xlim([fyr(1) fyr(end)])
        title('Forage')
        
        subplot(3,3,5)
        yyaxis left
        plot(yanom,squeeze(manom(j,lme,:)));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_mPa(lme,:));
        xlim([fyr(1) fyr(end)])
        title('Lg Pelagic')
        
        subplot(3,3,6)
        yyaxis left
        plot(yanom,squeeze(manom(j,lme,:)));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_mDa(lme,:));
        xlim([fyr(1) fyr(end)])
        title('Demersal')
        ylabel('Mean biomass (log_1_0 MT)')
        
        subplot(3,3,7)
        yyaxis left
        plot(yanom,squeeze(manom(j,lme,:)));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_mAa(lme,:));
        xlim([fyr(1) fyr(end)])
        xlabel('Time (y)')
        title('All Fish')
        
        subplot(3,3,9)
        yyaxis left
        plot(yanom,squeeze(manom(j,lme,:)));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_mba(lme,:));
        xlim([fyr(1) fyr(end)])
        title('Benthos')
        stamp('')
        print('-dpng',[ppath 'FOSI_',mod,ctex,'_',ltex,'_ts_corr.png'])
        
        %%
        figure(2)
        clf
        subplot(3,3,1)
        stem(lagsS,cS,'k')
        xlim([0 9])
        title('Small')
        
        subplot(3,3,2)
        stem(lagsM,cM,'k')
        xlim([0 9])
        str = {[ctex,' ', ltex], ' Medium'};
        title(str)
        
        subplot(3,3,3)
        stem(lagsL,cL,'k')
        xlim([0 9])
        title('Large')
        
        subplot(3,3,4)
        stem(lagsF,cF,'k')
        xlim([0 9])
        title('Forage')
        
        subplot(3,3,5)
        stem(lagsP,cP,'k')
        xlim([0 9])
        title('Lg Pelagic')
        
        subplot(3,3,6)
        stem(lagsD,cD,'k')
        xlim([0 9])
        title('Demersal')
        
        subplot(3,3,7)
        stem(lagsA,cA,'k')
        xlim([0 9])
        title('All Fish')
        
        subplot(3,3,9)
        stem(lagsB,cB,'k')
        xlim([0 9])
        title('Benthos')
        stamp('')
        print('-dpng',[ppath 'FOSI_',mod,ctex,'_',ltex,'_ts_crosscorr.png'])
        
    end
end

%% Try looking at corrs with 0-5 yr lag
yr = 0:5;
ns=0;
nm=0;
nl=0;
nf=0;
np=0;
nd=0;
na=0;
nb=0;

time = 1:length(fyr);
for i=1:length(lid)
    for j=1:length(canom)
        lme = lid(i);
        ltex = lname{i};
        ctex = canom{j};
        
        for k=1:6
            t = yr(k);
            %% Corr at diff lags
            [rs,ps] = corrcoef(lme_mSa(lme,(1+t):end),manom(j,lme,1:(end-t)));
            rS(i,j,k) = rs(1,2); pS(i,j,k) = ps(1,2);
            
            [rm,pm] = corrcoef(lme_mMa(lme,(1+t):end),manom(j,lme,1:(end-t)));
            rM(i,j,k) = rm(1,2); pM(i,j,k) = pm(1,2);
            
            [rl,pl] = corrcoef(lme_mLa(lme,(1+t):end),manom(j,lme,1:(end-t)));
            rL(i,j,k) = rl(1,2); pL(i,j,k) = pl(1,2);
            
            [rf,pf] = corrcoef(lme_mFa(lme,(1+t):end),manom(j,lme,1:(end-t)));
            rF(i,j,k) = rf(1,2); pF(i,j,k) = pf(1,2);
            
            [rp,pp] = corrcoef(lme_mPa(lme,(1+t):end),manom(j,lme,1:(end-t)));
            rP(i,j,k) = rp(1,2); pP(i,j,k) = pp(1,2);
            
            [rd,pd] = corrcoef(lme_mDa(lme,(1+t):end),manom(j,lme,1:(end-t)));
            rD(i,j,k) = rd(1,2); pD(i,j,k) = pd(1,2);
            
            [ra,pa] = corrcoef(lme_mAa(lme,(1+t):end),manom(j,lme,1:(end-t)));
            rA(i,j,k) = ra(1,2); pA(i,j,k) = pa(1,2);
            
            [rb,pb] = corrcoef(lme_mba(lme,(1+t):end),manom(j,lme,1:(end-t)));
            rB(i,j,k) = rb(1,2); pB(i,j,k) = pb(1,2);
            
            if (ps(1,2)<=0.05)
                ns=ns+1;
                sigS{ns,1} = ltex;
                sigS{ns,2} = ctex;
                sigS{ns,3} = t;
                sigS{ns,4} = rs(1,2);
                sigS{ns,5} = ps(1,2);
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
            if (pf(1,2)<=0.05)
                nf=nf+1;
                sigF{nf,1} = ltex;
                sigF{nf,2} = ctex;
                sigF{nf,3} = t;
                sigF{nf,4} = rf(1,2);
                sigF{nf,5} = pf(1,2);
            end
            if (pp(1,2)<=0.05)
                np=np+1;
                sigP{np,1} = ltex;
                sigP{np,2} = ctex;
                sigP{np,3} = t;
                sigP{np,4} = rp(1,2);
                sigP{np,5} = pp(1,2);
            end
            if (pd(1,2)<=0.05)
                nd=nd+1;
                sigD{nd,1} = ltex;
                sigD{nd,2} = ctex;
                sigD{nd,3} = t;
                sigD{nd,4} = rd(1,2);
                sigD{nd,5} = pd(1,2);
            end
            if (pa(1,2)<=0.05)
                na=na+1;
                sigA{na,1} = ltex;
                sigA{na,2} = ctex;
                sigA{na,3} = t;
                sigA{na,4} = ra(1,2);
                sigA{na,5} = pa(1,2);
            end
            if (pb(1,2)<=0.05)
                nb=nb+1;
                sigB{nb,1} = ltex;
                sigB{nb,2} = ctex;
                sigB{nb,3} = t;
                sigB{nb,4} = rb(1,2);
                sigB{nb,5} = pb(1,2);
            end
            
        end
    end
end

%%
writecell(sigS,[dpath 'LME_fosi_fished_',mod,'sigS_inputs.csv']);
writecell(sigM,[dpath 'LME_fosi_fished_',mod,'sigM_inputs.csv']);
writecell(sigL,[dpath 'LME_fosi_fished_',mod,'sigL_inputs.csv']);
writecell(sigF,[dpath 'LME_fosi_fished_',mod,'sigF_inputs.csv']);
writecell(sigP,[dpath 'LME_fosi_fished_',mod,'sigP_inputs.csv']);
writecell(sigD,[dpath 'LME_fosi_fished_',mod,'sigD_inputs.csv']);
writecell(sigA,[dpath 'LME_fosi_fished_',mod,'sigA_inputs.csv']);
writecell(sigB,[dpath 'LME_fosi_fished_',mod,'sigB_inputs.csv']);











