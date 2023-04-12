% Corr LME biomass of FEISTY with climate anoms
% CESM FOSI

% AO - DJF mean
% NAO - DJF mean
% Nino3.4 - DJF mean
% AMO - JAS mean
% PDO - whole annual mean

clear 
close all

apath = '/Users/cpetrik/Dropbox/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'climate_anomalies.mat'])

%% Remove N/A values
AMO(abs(AMO)>9) = nan;
AO(abs(AO)>9) = nan;
NAO(abs(NAO)>9) = nan;
Nino34(abs(Nino34)>9) = nan;
PDO(abs(PDO)>9) = nan;

%% Put previous Dec before Jan & Feb

wAO = AO(1:(end-1),12);
wNAO = NAO(1:(end-1),12);
wNino34 = Nino34(1:(end-1),12);

wAO(:,2:3) = AO(2:end,[1,2]);
wNAO(:,2:3) = NAO(2:end,[1,2]);
wNino34(:,2:3) = Nino34(2:end,[1,2]);

AOyr = AOyr(2:end);
NAOyr = NAOyr(2:end);
Nino34yr = Nino34yr(2:end);

%% Climate anom DJF means
mAO = nanmean(wAO,2);
mNAO = nanmean(wNAO,2);
mNino34 = nanmean(wNino34,2);

%%
jas = [7:9];
mAMO = nanmean(AMO(:,jas),2);

mPDO = nanmean(PDO,2);


%% Isolate years of interest 1948-2015
fyr = 1948:2015;
mAMO = mAMO(AMOyr>=1948 & AMOyr<=2015);
mAO = mAO(AOyr>=1948 & AOyr<=2015);
mNAO = mNAO(NAOyr>=1948 & NAOyr<=2015);
mNino34 = mNino34(Nino34yr>=1948 & Nino34yr<=2015);
mPDO = mPDO(PDOyr>=1948 & PDOyr<=2015);

yAMO = AMOyr(AMOyr>=1948 & AMOyr<=2015);
yAO = AOyr(AOyr>=1948 & AOyr<=2015);
yNAO = NAOyr(NAOyr>=1948 & NAOyr<=2015);
yNino34 = Nino34yr(Nino34yr>=1948 & Nino34yr<=2015);
yPDO = PDOyr(PDOyr>=1948 & PDOyr<=2015);

%% put in matrix
manom = nan*ones(5,68);
manom(1,:) = mAMO;
manom(2,4:end) = mAO;
manom(3,:) = mNAO;
manom(4,:) = mNino34;
manom(5,:) = mPDO;

yanom = 1948:2015;

tanom = {'AMO','AO','NAO','Nino34','PDO'};

%% id start and end years
yst = nan*ones(11,1);
yen = nan*ones(11,1);
for k=1:length(tanom)
    nn = find(~isnan(manom(k,:)));
    yst(k) = nn(1);
    yen(k) = nn(end);
end

%% save
save([apath 'Climate_anomalies_seasonal_means.mat'],'manom','yanom','tanom',...
    'yst','yen');

%% csv file for R
tmanom = manom';
tmanom(:,6) = yanom';
Tac = array2table(tmanom,'VariableNames',[tanom, 'Year']);

writetable(Tac,[apath 'Climate_anomalies_seasonal_means.csv'],'WriteRowNames',false);










