% Corr LME biomass of FEISTY with climate anoms
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

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
dpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

sims = {'v14_All_fish03_';'v14_climatol_';'v14_varFood_';'v14_varTemp_'};
mod = sims{1};

load([dpath 'LME_fosi_fished_',mod,cfile '.mat']);

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

%% Pac
figure
plot(yPDO,mPDO,'k'); hold on;
plot(yNOI,mNOI,'r'); hold on;
plot(yMEI,mMEI,'b'); hold on;
plot(yNino12,mNino12,'m'); hold on;
plot(yNino34,mNino34,'c'); hold on;

%% Atl
figure
plot(yAMO,mAMO,'k'); hold on;
plot(yNAO,mNAO,'b'); hold on;

%% Group types and sizes
lme_mFb = lme_msfb + lme_mmfb;
lme_mPb = lme_mspb + lme_mmpb + lme_mlpb;
lme_mDb = lme_msdb + lme_mmdb + lme_mldb;
lme_mSb = lme_msfb + lme_mspb + lme_msdb;
lme_mMb = lme_mmfb + lme_mmpb + lme_mmdb;
lme_mLb = lme_mlpb + lme_mldb;
lme_mAb = lme_mFb + lme_mPb + lme_mDb;

lme_sFb = lme_ssfb + lme_smfb;
lme_sPb = lme_sspb + lme_smpb + lme_slpb;
lme_sDb = lme_ssdb + lme_smdb + lme_sldb;
lme_sSb = lme_ssfb + lme_sspb + lme_ssdb;
lme_sMb = lme_smfb + lme_smpb + lme_smdb;
lme_sLb = lme_slpb + lme_sldb;
lme_sAb = lme_sFb + lme_sPb + lme_sDb;

%% Calc fish anomalies
%mean biomass
lme_msfa = lme_msfb - nanmean(lme_msfb,2);
lme_mspa = lme_mspb - nanmean(lme_mspb,2);
lme_msda = lme_msdb - nanmean(lme_msdb,2);
lme_mmfa = lme_mmfb - nanmean(lme_mmfb,2);
lme_mmpa = lme_mmpb - nanmean(lme_mmpb,2);
lme_mmda = lme_mmdb - nanmean(lme_mmdb,2);
lme_mlpa = lme_mlpb - nanmean(lme_mlpb,2);
lme_mlda = lme_mldb - nanmean(lme_mldb,2);
lme_mba = lme_mbb - nanmean(lme_mbb,2);
lme_mFa = lme_mFb - nanmean(lme_mFb,2);
lme_mPa = lme_mPb - nanmean(lme_mPb,2);
lme_mDa = lme_mDb - nanmean(lme_mDb,2);
lme_mSa = lme_mSb - nanmean(lme_mSb,2);
lme_mMa = lme_mMb - nanmean(lme_mMb,2);
lme_mLa = lme_mLb - nanmean(lme_mLb,2);
lme_mAa = lme_mAb - nanmean(lme_mAb,2);

%total biomass
lme_ssfa = lme_ssfb - nanmean(lme_ssfb,2);
lme_sspa = lme_sspb - nanmean(lme_sspb,2);
lme_ssda = lme_ssdb - nanmean(lme_ssdb,2);
lme_smfa = lme_smfb - nanmean(lme_smfb,2);
lme_smpa = lme_smpb - nanmean(lme_smpb,2);
lme_smda = lme_smdb - nanmean(lme_smdb,2);
lme_slpa = lme_slpb - nanmean(lme_slpb,2);
lme_slda = lme_sldb - nanmean(lme_sldb,2);
lme_sba = lme_sbb - nanmean(lme_sbb,2);
lme_sFa = lme_sFb - nanmean(lme_sFb,2);
lme_sPa = lme_sPb - nanmean(lme_sPb,2);
lme_sDa = lme_sDb - nanmean(lme_sDb,2);
lme_sSa = lme_sSb - nanmean(lme_sSb,2);
lme_sMa = lme_sMb - nanmean(lme_sMb,2);
lme_sLa = lme_sLb - nanmean(lme_sLb,2);
lme_sAa = lme_sAb - nanmean(lme_sAb,2);

%% Test PDO in CCE
% US LMEs only
lid = [54,1:2,10,3,5:7];
lname = {'CHK','EBS','GAK','HI','CCE','GMX','SE','NE'};

lme = 3;
[c,lags] = xcorr(mPDO,lme_mFb(lme,:),20);
figure
stem(lags,c,'k')
[R,P] = corrcoef(mPDO,lme_mFb(lme,:))

%% CCE - PDO ==========================================================
close all
colororder({'k','b'})

lme = 3;

figure(1)
subplot(3,3,1)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(3,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(3,:));
xlim([fyr(1) fyr(end)])
title(['PDO CCE'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(3,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
ylabel('PDO')
yyaxis right
plot(fyr,lme_mFa(3,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(3,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(3,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(3,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(3,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_CCE_PDO_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mPDO,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mPDO,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mPDO,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mPDO,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mPDO,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mPDO,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mPDO,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mPDO,20);

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['PDO CCE'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_CCE_PDO_ts_crosscorr.png'])

%% CCE - MEI
figure(3)
subplot(3,3,1)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(3,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(3,:));
xlim([fyr(1) fyr(end)])
title(['MEI CCE'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(3,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
ylabel('MEI')
yyaxis right
plot(fyr,lme_mFa(3,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(3,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(3,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(3,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(3,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_CCE_MEI_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mMEI,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mMEI,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mMEI,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mMEI,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mMEI,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mMEI,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mMEI,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mMEI,20);

figure(4)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['MEI CCE'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_CCE_MEI_ts_crosscorr.png'])

%% CCE - NOI
figure(5)
subplot(3,3,1)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(3,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(3,:));
xlim([fyr(1) fyr(end)])
title(['NOI CCE'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(3,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
ylabel('NOI')
yyaxis right
plot(fyr,lme_mFa(3,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(3,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(3,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(3,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(3,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_CCE_NOI_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mNOI,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mNOI,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mNOI,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mNOI,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mNOI,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mNOI,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mNOI,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mNOI,20);

figure(6)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['NOI CCE'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_CCE_NOI_ts_crosscorr.png'])

%% CCE - ENSO 1-2
figure(7)
subplot(3,3,1)
yyaxis left
plot(yNino12,mNino12);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(3,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yNino12,mNino12);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(3,:));
xlim([fyr(1) fyr(end)])
title(['Nino12 CCE'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yNino12,mNino12);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(3,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yNino12,mNino12);
xlim([fyr(1) fyr(end)])
ylabel('Nino12')
yyaxis right
plot(fyr,lme_mFa(3,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yNino12,mNino12);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(3,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yNino12,mNino12);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(3,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yNino12,mNino12);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(3,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yNino12,mNino12);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(3,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_CCE_Nino12_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mNino12,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mNino12,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mNino12,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mNino12,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mNino12,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mNino12,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mNino12,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mNino12,20);

figure(8)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['Nino12 CCE'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_CCE_Nino12_ts_crosscorr.png'])

%% CCE - ENSO 3-4
figure(9)
subplot(3,3,1)
yyaxis left
plot(yNino34,mNino34);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(3,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yNino34,mNino34);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(3,:));
xlim([fyr(1) fyr(end)])
title(['Nino34 CCE'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yNino34,mNino34);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(3,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yNino34,mNino34);
xlim([fyr(1) fyr(end)])
ylabel('Nino34')
yyaxis right
plot(fyr,lme_mFa(3,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yNino34,mNino34);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(3,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yNino34,mNino34);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(3,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yNino34,mNino34);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(3,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yNino34,mNino34);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(3,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_CCE_Nino34_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mNino34,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mNino34,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mNino34,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mNino34,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mNino34,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mNino34,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mNino34,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mNino34,20);

figure(10)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['Nino34 CCE'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_CCE_Nino34_ts_crosscorr.png'])

%% EBS - PDO ==========================================================
lme = 1;

figure(1)
subplot(3,3,1)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(1,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(1,:));
xlim([fyr(1) fyr(end)])
title(['PDO EBS'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(1,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
ylabel('PDO')
yyaxis right
plot(fyr,lme_mFa(1,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(1,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(1,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(1,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(1,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_EBS_PDO_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mPDO,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mPDO,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mPDO,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mPDO,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mPDO,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mPDO,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mPDO,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mPDO,20);

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['PDO EBS'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_EBS_PDO_ts_crosscorr.png'])

%% EBS - MEI
figure(1)
subplot(3,3,1)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(1,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(1,:));
xlim([fyr(1) fyr(end)])
title(['MEI EBS'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(1,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
ylabel('MEI')
yyaxis right
plot(fyr,lme_mFa(1,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(1,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(1,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(1,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(1,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_EBS_MEI_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mMEI,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mMEI,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mMEI,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mMEI,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mMEI,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mMEI,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mMEI,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mMEI,20);

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['MEI EBS'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_EBS_MEI_ts_crosscorr.png'])

%% EBS - NOI
figure(1)
subplot(3,3,1)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(1,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(1,:));
xlim([fyr(1) fyr(end)])
title(['NOI EBS'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(1,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
ylabel('NOI')
yyaxis right
plot(fyr,lme_mFa(1,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(1,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(1,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(1,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(1,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_EBS_NOI_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mNOI,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mNOI,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mNOI,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mNOI,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mNOI,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mNOI,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mNOI,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mNOI,20);

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['NOI EBS'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_EBS_NOI_ts_crosscorr.png'])

%% EBS - AO
figure(1)
subplot(3,3,1)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(1,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(1,:));
xlim([fyr(1) fyr(end)])
title(['AO EBS'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(1,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
ylabel('AO')
yyaxis right
plot(fyr,lme_mFa(1,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(1,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(1,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(1,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(1,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_EBS_AO_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mAO,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mAO,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mAO,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mAO,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mAO,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mAO,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mAO,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mAO,20);

stem(lags,c,'k')

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['AO EBS'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_EBS_AO_ts_crosscorr.png'])

%% GAK - PDO ==========================================================
lme = 2;
figure(1)
subplot(3,3,1)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(2,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(2,:));
xlim([fyr(1) fyr(end)])
title(['PDO GAK'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(2,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
ylabel('PDO')
yyaxis right
plot(fyr,lme_mFa(2,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(2,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(2,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(2,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(2,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_PDO_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mPDO,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mPDO,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mPDO,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mPDO,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mPDO,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mPDO,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mPDO,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mPDO,20);

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['PDO GAK'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_PDO_ts_crosscorr.png'])

%% GAK - MEI
figure(1)
subplot(3,3,1)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(2,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(2,:));
xlim([fyr(1) fyr(end)])
title(['MEI GAK'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(2,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
ylabel('MEI')
yyaxis right
plot(fyr,lme_mFa(2,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(2,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(2,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(2,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(2,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_MEI_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mMEI,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mMEI,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mMEI,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mMEI,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mMEI,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mMEI,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mMEI,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mMEI,20);

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['MEI GAK'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_MEI_ts_crosscorr.png'])

%% GAK - NOI
figure(1)
subplot(3,3,1)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(2,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(2,:));
xlim([fyr(1) fyr(end)])
title(['NOI GAK'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(2,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
ylabel('NOI')
yyaxis right
plot(fyr,lme_mFa(2,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(2,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(2,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(2,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(2,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_NOI_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mNOI,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mNOI,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mNOI,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mNOI,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mNOI,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mNOI,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mNOI,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mNOI,20);

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['NOI GAK'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_NOI_ts_crosscorr.png'])

%% GAK - AO
figure(1)
subplot(3,3,1)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(2,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(2,:));
xlim([fyr(1) fyr(end)])
title(['AO GAK'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(2,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
ylabel('AO')
yyaxis right
plot(fyr,lme_mFa(2,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(2,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(2,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(2,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(2,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_AO_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mAO,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mAO,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mAO,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mAO,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mAO,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mAO,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mAO,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mAO,20);

stem(lags,c,'k')

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['AO GAK'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_AO_ts_crosscorr.png'])

%% GAK - PDO ==========================================================
lme = 2;
figure(1)
subplot(3,3,1)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(2,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(2,:));
xlim([fyr(1) fyr(end)])
title(['PDO GAK'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(2,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
ylabel('PDO')
yyaxis right
plot(fyr,lme_mFa(2,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(2,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(2,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(2,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yPDO,mPDO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(2,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_PDO_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mPDO,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mPDO,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mPDO,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mPDO,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mPDO,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mPDO,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mPDO,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mPDO,20);

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['PDO GAK'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_PDO_ts_crosscorr.png'])

%% GAK - MEI
figure(1)
subplot(3,3,1)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(2,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(2,:));
xlim([fyr(1) fyr(end)])
title(['MEI GAK'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(2,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
ylabel('MEI')
yyaxis right
plot(fyr,lme_mFa(2,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(2,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(2,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(2,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yMEI,mMEI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(2,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_MEI_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mMEI,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mMEI,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mMEI,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mMEI,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mMEI,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mMEI,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mMEI,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mMEI,20);

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['MEI GAK'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_MEI_ts_crosscorr.png'])

%% GAK - NOI
figure(1)
subplot(3,3,1)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(2,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(2,:));
xlim([fyr(1) fyr(end)])
title(['NOI GAK'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(2,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
ylabel('NOI')
yyaxis right
plot(fyr,lme_mFa(2,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(2,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(2,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(2,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yNOI,mNOI);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(2,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_NOI_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mNOI,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mNOI,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mNOI,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mNOI,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mNOI,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mNOI,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mNOI,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mNOI,20);

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['NOI GAK'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_NOI_ts_crosscorr.png'])

%% GAK - AO
figure(1)
subplot(3,3,1)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mSa(2,:));
xlim([fyr(1) fyr(end)])
title('Small')

subplot(3,3,2)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mMa(2,:));
xlim([fyr(1) fyr(end)])
title(['AO GAK'; ' Medium'])

subplot(3,3,3)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mLa(2,:));
xlim([fyr(1) fyr(end)])
title('Large')

subplot(3,3,4)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
ylabel('AO')
yyaxis right
plot(fyr,lme_mFa(2,:));
xlim([fyr(1) fyr(end)])
title('Forage')

subplot(3,3,5)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mPa(2,:));
xlim([fyr(1) fyr(end)])
title('Lg Pelagic')

subplot(3,3,6)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mDa(2,:));
xlim([fyr(1) fyr(end)])
title('Demersal')
ylabel('Mean biomass (log_1_0 MT)')

subplot(3,3,7)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mAa(2,:));
xlim([fyr(1) fyr(end)])
xlabel('Time (y)')
title('All Fish')

subplot(3,3,9)
yyaxis left
plot(yAO,mAO);
xlim([fyr(1) fyr(end)])
yyaxis right
plot(fyr,lme_mba(2,:));
xlim([fyr(1) fyr(end)])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_AO_ts_corr.png'])

% Cross corr
[cS,lagsS] = xcorr(lme_mSb(lme,:),mAO,20);
[cM,lagsM] = xcorr(lme_mMb(lme,:),mAO,20);
[cL,lagsL] = xcorr(lme_mLb(lme,:),mAO,20);
[cF,lagsF] = xcorr(lme_mFb(lme,:),mAO,20);
[cP,lagsP] = xcorr(lme_mPb(lme,:),mAO,20);
[cD,lagsD] = xcorr(lme_mDb(lme,:),mAO,20);
[cA,lagsA] = xcorr(lme_mAb(lme,:),mAO,20);
[cB,lagsB] = xcorr(lme_mbb(lme,:),mAO,20);

stem(lags,c,'k')

figure(2)
subplot(3,3,1)
stem(lagsS,cS,'k')
xlim([0 15])
title('Small')

subplot(3,3,2)
stem(lagsM,cM,'k')
xlim([0 15])
title(['AO GAK'; ' Medium'])

subplot(3,3,3)
stem(lagsL,cL,'k')
xlim([0 15])
title('Large')

subplot(3,3,4)
stem(lagsF,cF,'k')
xlim([0 15])
title('Forage')

subplot(3,3,5)
stem(lagsP,cP,'k')
xlim([0 15])
title('Lg Pelagic')

subplot(3,3,6)
stem(lagsD,cD,'k')
xlim([0 15])
title('Demersal')

subplot(3,3,7)
stem(lagsA,cA,'k')
xlim([0 15])
title('All Fish')

subplot(3,3,9)
stem(lagsB,cB,'k')
xlim([0 15])
title('Benthos')
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_GAK_AO_ts_crosscorr.png'])





