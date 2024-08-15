% Plot climate indices by month to see min/max

clear 
close all

apath = '/Users/cpetrik/Dropbox/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'climate_anomalies.mat'])

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

%% AMO
absAMO = abs(AMO);
maAMO = nanmean(absAMO,1);

figure(1)
plot(1:12,maAMO) %max in summer (7,8,9)
title('AMO')

%% AO
absAO = abs(AO);
maAO = nanmean(absAO,1);

figure(2)
plot(1:12,maAO) %max winter (12,1,2)
title('AO')

%% NAO
absNAO = abs(NAO);
maNAO = nanmean(absNAO,1);

figure(3)
plot(1:12,maNAO) %max winter (12,1,2)
title('NAO')

%% ENSO
absNino34 = abs(Nino34);
maNino34 = nanmean(absNino34,1);

figure(4)
plot(1:12,maNino34) %max winter (12,1,2)
title('ENSO')

%% PDO
absPDO = abs(PDO);
maPDO = nanmean(absPDO,1);

figure(5)
plot(1:12,maPDO) %no clear pattern
title('PDO')

%% All together
figure(6)
subplot(3,2,1)
plot(1:12,maNino34) %max winter (12,1,2)
title('ENSO')

subplot(3,2,2)
plot(1:12,maPDO) %no clear pattern
title('PDO')

subplot(3,2,3)
plot(1:12,maAMO) %max in summer (7,8,9)
title('AMO')

subplot(3,2,4)
plot(1:12,maNAO) %max winter (12,1,2)
title('NAO')

subplot(3,2,5)
plot(1:12,maAO) %max winter (12,1,2)
title('AO')
print('-dpng',[apath 'climate_anom_ts_monthly_abs_value.png'])





