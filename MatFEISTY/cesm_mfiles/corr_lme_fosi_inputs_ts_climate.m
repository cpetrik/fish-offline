% Corr FEISTY inputs LME means with climate anoms
% CESM FOSI

clear all
close all

apath = '/Users/cpetrik/Dropbox/Princeton/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'climate_anomalies.mat'])

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
load([cpath 'lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'])

ppath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

%% Remove N/A values
AMO(abs(AMO)>9) = nan;
AO(abs(AO)>9) = nan;
MEI(abs(MEI)>9) = nan;
NAO(abs(NAO)>9) = nan;
Nino12(abs(Nino12)>9) = nan;
Nino34(abs(Nino34)>9) = nan;
Nino3(abs(Nino3)>9) = nan;
Nino4(abs(Nino4)>9) = nan;
NOI(abs(NOI)>9) = nan;
PDO(abs(PDO)>9) = nan;
SOI(abs(SOI)>9) = nan;

%% Climate anom annual means

mAMO = nanmean(AMO,2);
mAO = nanmean(AO,2);
mMEI = nanmean(MEI,2);
mNAO = nanmean(NAO,2);
mNino12 = nanmean(Nino12,2);
mNino34 = nanmean(Nino34,2);
mNino3 = nanmean(Nino3,2);
mNino4 = nanmean(Nino4,2);
mNOI = nanmean(NOI,2);
mPDO = nanmean(PDO,2);
mSOI = nanmean(SOI,2);

%% Isolate years of interest 1948-2015
fyr = 1948:2015;
mAMO = mAMO(AMOyr>=1948 & AMOyr<=2015);
mAO = mAO(AOyr>=1948 & AOyr<=2015);
mMEI = mMEI(MEIyr>=1948 & MEIyr<=2015);
mNAO = mNAO(NAOyr>=1948 & NAOyr<=2015);
mNino12 = mNino12(Nino12yr>=1948 & Nino12yr<=2015);
mNino34 = mNino34(Nino34yr>=1948 & Nino34yr<=2015);
mNino3 = mNino3(Nino3yr>=1948 & Nino3yr<=2015);
mNino4 = mNino4(Nino4yr>=1948 & Nino4yr<=2015);
mNOI = mNOI(NOIyr>=1948 & NOIyr<=2015);
mPDO = mPDO(PDOyr>=1948 & PDOyr<=2015);
mSOI = mSOI(SOIyr>=1948 & SOIyr<=2015);

yAMO = AMOyr(AMOyr>=1948 & AMOyr<=2015);
yAO = AOyr(AOyr>=1948 & AOyr<=2015);
yMEI = MEIyr(MEIyr>=1948 & MEIyr<=2015);
yNAO = NAOyr(NAOyr>=1948 & NAOyr<=2015);
yNino12 = Nino12yr(Nino12yr>=1948 & Nino12yr<=2015);
yNino34 = Nino34yr(Nino34yr>=1948 & Nino34yr<=2015);
yNino3 = Nino3yr(Nino3yr>=1948 & Nino3yr<=2015);
yNino4 = Nino4yr(Nino4yr>=1948 & Nino4yr<=2015);
yNOI = NOIyr(NOIyr>=1948 & NOIyr<=2015);
yPDO = PDOyr(PDOyr>=1948 & PDOyr<=2015);
ySOI = SOIyr(SOIyr>=1948 & SOIyr<=2015);

%% put in matrix
manom = nan*ones(11,68);
manom(1,:) = mAMO;
manom(2,3:end) = mAO;
manom(3,32:end) = mMEI;
manom(4,:) = mNAO;
manom(5,:) = mNino12;
manom(6,:) = mNino34;
manom(7,:) = mNino3;
manom(8,:) = mNino4;
manom(9,1:60) = mNOI;
manom(10,:) = mPDO;
manom(11,:) = mSOI;

yanom = 1948:2015;

canom = {'AMO','AO','MEI','NAO','Nino12','Nino34','Nino3','Nino4','NOI',...
    'PDO','SOI'};

%% id start and end years
yst = nan*ones(11,1);
yen = nan*ones(11,1);
for k=1:length(canom)
    nn = find(~isnan(manom(k,:)));
    yst(k) = nn(1);
    yen(k) = nn(end);
end


%% Calc anomalies
%means
%lme_tp_fosi','lme_tb_fosi','lme_det_fosi','lme_mz_fosi','lme_mzloss_fosi
lme_tpa = lme_tp_fosi - nanmean(lme_tp_fosi,2);
lme_tba = lme_tb_fosi - nanmean(lme_tb_fosi,2);
lme_deta = lme_det_fosi - nanmean(lme_det_fosi,2);
lme_mza = lme_mz_fosi - nanmean(lme_mz_fosi,2);
lme_losa = lme_mzloss_fosi - nanmean(lme_mzloss_fosi,2);

%% Loop over all LMEs and all Climate
close all
colororder({'k','b'})

for i=1%:length(lid)
    for j=2%1:length(canom)
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
        yyaxis right
        plot(fyr,lme_mSa(3,:));
        xlim([fyr(1) fyr(end)])
        title('Tpelagic')
        
        subplot(2,3,2)
        yyaxis left
        plot(yanom,manom(j,:));
        xlim([fyr(1) fyr(end)])
        ylabel('PDO')
        yyaxis right
        plot(fyr,lme_mFa(3,:));
        xlim([fyr(1) fyr(end)])
        title([ctex,' ', ltex ' LzooC'])
        
        subplot(2,3,3)
        yyaxis left
        plot(yanom,manom(j,:));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_mPa(3,:));
        xlim([fyr(1) fyr(end)])
        title('Lzoo loss')
        
        subplot(2,3,4)
        yyaxis left
        plot(yanom,manom(j,:));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_mMa(3,:));
        xlim([fyr(1) fyr(end)])
        title('Tbottom')
        
        subplot(2,3,5)
        yyaxis left
        plot(yanom,manom(j,:));
        xlim([fyr(1) fyr(end)])
        yyaxis right
        plot(fyr,lme_mLa(3,:));
        xlim([fyr(1) fyr(end)])
        title('Detritus')
        stamp('')
        print('-dpng',[ppath 'FOSI_inputs_',ctex,'_',ltex,'_ts_corr.png'])
        
        %%
        figure(2)
        clf
        subplot(2,3,1)
        stem(lagsP,cP,'k')
        xlim([0 15])
        title('Tpelagic')
        
        subplot(2,3,2)
        stem(lagsM,cM,'k')
        xlim([0 15])
        title([ctex,' ', ltex ' LzooC'])
        
        subplot(2,3,3)
        stem(lagsL,cL,'k')
        xlim([0 15])
        title('Lzoo loss')
        
        subplot(2,3,4)
        stem(lagsB,cB,'k')
        xlim([0 15])
        title('Tbottom')
        
        subplot(2,3,5)
        stem(lagsD,cD,'k')
        xlim([0 15])
        title('Detritus')
        stamp('')
        print('-dpng',[ppath 'FOSI_inputs_',ctex,'_',ltex,'_ts_crosscorr.png'])
        
    end
end