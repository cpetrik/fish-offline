% Corr LME biomass of FEISTY with climate anoms
% CESM FOSI

clear all
close all

apath = '/Users/cpetrik/Dropbox/Princeton/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
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

%% save
save([apath 'Climate_anomalies_annual_means.mat'],'manom','yanom','canom',...
    'yst','yen');

%% csv file for R
tmanom = manom';
tmanom(:,12) = yanom';
Tac = array2table(tmanom,'VariableNames',[canom, 'Year']);

writetable(Tac,[apath 'Climate_anomalies_annual_means.csv'],'WriteRowNames',false);










