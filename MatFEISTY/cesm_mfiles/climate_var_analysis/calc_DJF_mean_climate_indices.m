% Corr LME biomass of FEISTY with climate anoms
% CESM FOSI

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

%% Put previous Dec before Jan & Feb

wAMO = AMO(1:(end-1),12);
wAO = AO(1:(end-1),12);
wMEI = MEI(1:(end-1),12);
wNAO = NAO(1:(end-1),12);
wNino12 = Nino12(1:(end-1),12);
wNino34 = Nino34(1:(end-1),12);
wNino3 = Nino3(1:(end-1),12);
wNino4 = Nino4(1:(end-1),12);
wNOI = NOI(1:(end-1),12);
wPDO = PDO(1:(end-1),12);
wSOI = SOI(1:(end-1),12);

wAMO(:,2:3) = AMO(2:end,[1,2]);
wAO(:,2:3) = AO(2:end,[1,2]);
wMEI(:,2:3) = MEI(2:end,[1,2]);
wNAO(:,2:3) = NAO(2:end,[1,2]);
wNino12(:,2:3) = Nino12(2:end,[1,2]);
wNino34(:,2:3) = Nino34(2:end,[1,2]);
wNino3(:,2:3) = Nino3(2:end,[1,2]);
wNino4(:,2:3) = Nino4(2:end,[1,2]);
wNOI(:,2:3) = NOI(2:end,[1,2]);
wPDO(:,2:3) = PDO(2:end,[1,2]);
wSOI(:,2:3) = SOI(2:end,[1,2]);

AMOyr = AMOyr(2:end);
AOyr = AOyr(2:end);
MEIyr = MEIyr(2:end);
NAOyr = NAOyr(2:end);
Nino12yr = Nino12yr(2:end);
Nino34yr = Nino34yr(2:end);
Nino3yr = Nino3yr(2:end);
Nino4yr = Nino4yr(2:end);
NOIyr = NOIyr(2:end);
PDOyr = PDOyr(2:end);
SOIyr = SOIyr(2:end);

%% Climate anom DJF means
mAMO = nanmean(wAMO,2);
mAO = nanmean(wAO,2);
mMEI = nanmean(wMEI,2);
mNAO = nanmean(wNAO,2);
mNino12 = nanmean(wNino12,2);
mNino34 = nanmean(wNino34,2);
mNino3 = nanmean(wNino3,2);
mNino4 = nanmean(wNino4,2);
mNOI = nanmean(wNOI,2);
mPDO = nanmean(wPDO,2);
mSOI = nanmean(wSOI,2);

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
manom(2,4:end) = mAO;
manom(3,33:end) = mMEI;
manom(4,:) = mNAO;
manom(5,:) = mNino12;
manom(6,:) = mNino34;
manom(7,:) = mNino3;
manom(8,:) = mNino4;
manom(9,2:60) = mNOI;
manom(10,:) = mPDO;
manom(11,:) = mSOI;

yanom = 1948:2015;

tanom = {'AMO','AO','MEI','NAO','Nino12','Nino34','Nino3','Nino4','NOI',...
    'PDO','SOI'};

%% id start and end years
yst = nan*ones(11,1);
yen = nan*ones(11,1);
for k=1:length(tanom)
    nn = find(~isnan(manom(k,:)));
    yst(k) = nn(1);
    yen(k) = nn(end);
end

%% save
save([apath 'Climate_anomalies_DJF_means.mat'],'manom','yanom','tanom',...
    'yst','yen');

%% csv file for R
tmanom = manom';
tmanom(:,12) = yanom';
Tac = array2table(tmanom,'VariableNames',[tanom, 'Year']);

writetable(Tac,[apath 'Climate_anomalies_DJF_means.csv'],'WriteRowNames',false);










