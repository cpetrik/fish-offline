% Corr FEISTY inputs LME means with climate anoms
% CESM FOSI

clear all
close all

apath = '/Users/cpetrik/Dropbox/Princeton/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'Climate_anomalies_annual_means.mat']);

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
load([fpath 'FOSI_lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'])

ppath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

%% Annual means
[nid,nt] = size(lme_mtp);
nyr = nt/12;
st=1:12:nt;
en=12:12:nt;

for n=1:length(st)
    % mean 
    lme_amtp(:,n)=nanmean(lme_mtp(:,st(n):en(n)),2);
    lme_amtb(:,n)=nanmean(lme_mtb(:,st(n):en(n)),2);
    lme_amdet(:,n)=nanmean(lme_mdet(:,st(n):en(n)),2);
    lme_ammz(:,n)=nanmean(lme_mmz(:,st(n):en(n)),2);
    lme_amloss(:,n)=nanmean(lme_mloss(:,st(n):en(n)),2);

end

%% Calc anomalies
%means
lme_tpa = lme_amtp - nanmean(lme_amtp,2);
lme_tba = lme_amtb - nanmean(lme_amtb,2);
lme_deta = lme_amdet - nanmean(lme_amdet,2);
lme_mza = lme_ammz - nanmean(lme_ammz,2);
lme_losa = lme_amloss - nanmean(lme_amloss,2);

%% SAVE
spath='/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/FOSI/';
if (~isfolder(spath))
    mkdir(spath)
end
save([spath 'LME_fosi_input_anomalies_annual.mat'],...
    'lme_tpa','lme_tba','lme_deta','lme_mza','lme_losa');
save([cpath 'LME_fosi_input_anomalies_annual.mat'],...
    'lme_tpa','lme_tba','lme_deta','lme_mza','lme_losa');

%% Loop over all LMEs and all Climate
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
        
        %% Cross corr - FIX TO BE SAME DATES
        [cP,lagsP] = xcorr(lme_tpa(lme,yst(j):yen(j)),manom(j,yst(j):yen(j)),15);
        [cB,lagsB] = xcorr(lme_tba(lme,yst(j):yen(j)),manom(j,yst(j):yen(j)),15);
        [cD,lagsD] = xcorr(lme_deta(lme,yst(j):yen(j)),manom(j,yst(j):yen(j)),15);
        [cM,lagsM] = xcorr(lme_mza(lme,yst(j):yen(j)),manom(j,yst(j):yen(j)),15);
        [cL,lagsL] = xcorr(lme_losa(lme,yst(j):yen(j)),manom(j,yst(j):yen(j)),15);
        
        %%
        figure(1)
        clf
        subplot(2,3,1)
        yyaxis left
        plot(yanom,manom(j,:));
        xlim([fyr(1) fyr(end)])
        ylabel(ctex)
        yyaxis right
        plot(fyr,lme_tpa(lme,:));
        xlim([fyr(1) fyr(end)])
        title('Tpelagic')
        
        subplot(2,3,2)
        yyaxis left
        plot(yanom,manom(j,:));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_losa(lme,:));
        xlim([fyr(1) fyr(end)])
        str = {[ctex,' ', ltex], ' LzooC'};
        title(str)
        
        subplot(2,3,3)
        yyaxis left
        plot(yanom,manom(j,:));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_mza(lme,:));
        xlim([fyr(1) fyr(end)])
        title('Lzoo loss')
        
        subplot(2,3,4)
        yyaxis left
        plot(yanom,manom(j,:));
        xlim([fyr(1) fyr(end)])
        ylabel(ctex)
        yyaxis right
        plot(fyr,lme_tba(lme,:));
        xlim([fyr(1) fyr(end)])
        title('Tbottom')
        
        subplot(2,3,5)
        yyaxis left
        plot(yanom,manom(j,:));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_deta(lme,:));
        xlim([fyr(1) fyr(end)])
        title('Detritus')
        stamp('')
        print('-dpng',[ppath 'FOSI_inputs_',ctex,'_',ltex,'_ts_corr.png'])
        
        %%
        figure(2)
        clf
        subplot(2,3,1)
        stem(lagsP,cP,'k')
        xlim([0 9])
        title('Tpelagic')
        
        subplot(2,3,2)
        stem(lagsM,cM,'k')
        xlim([0 9])
        str = {[ctex,' ', ltex], ' LzooC'};
        title(str)
        
        subplot(2,3,3)
        stem(lagsL,cL,'k')
        xlim([0 9])
        title('Lzoo loss')
        
        subplot(2,3,4)
        stem(lagsB,cB,'k')
        xlim([0 9])
        title('Tbottom')
        
        subplot(2,3,5)
        stem(lagsD,cD,'k')
        xlim([0 9])
        title('Detritus')
        stamp('')
        print('-dpng',[ppath 'FOSI_inputs_',ctex,'_',ltex,'_ts_crosscorr.png'])
        
    end
end

%% Try looking at corrs with 0-5 yr lag
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
            %% Corr at diff lags
            [rp,pp] = corrcoef(lme_tpa(lme,yst(j)+t:yen(j)),manom(j,yst(j):yen(j)-t));
            rP(i,j,k) = rp(1,2); pP(i,j,k) = pp(1,2);
            
            [rb,pb] = corrcoef(lme_tba(lme,yst(j)+t:yen(j)),manom(j,yst(j):yen(j)-t));
            rB(i,j,k) = rb(1,2); pB(i,j,k) = pb(1,2);
            
            [rd,pd] = corrcoef(lme_deta(lme,yst(j)+t:yen(j)),manom(j,yst(j):yen(j)-t));
            rD(i,j,k) = rd(1,2); pD(i,j,k) = pd(1,2);
            
            [rm,pm] = corrcoef(lme_mza(lme,yst(j)+t:yen(j)),manom(j,yst(j):yen(j)-t));
            rM(i,j,k) = rm(1,2); pM(i,j,k) = pm(1,2);
            
            [rl,pl] = corrcoef(lme_losa(lme,yst(j)+t:yen(j)),manom(j,yst(j):yen(j)-t));
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
writecell(sigM,[cpath 'LME_fosi_inputs_sigLzooC.csv']);
writecell(sigL,[cpath 'LME_fosi_inputs_sigZooLoss.csv']);
writecell(sigP,[cpath 'LME_fosi_inputs_sigTP.csv']);
writecell(sigD,[cpath 'LME_fosi_inputs_sigDet.csv']);
writecell(sigB,[cpath 'LME_fosi_inputs_sigTB.csv']);








