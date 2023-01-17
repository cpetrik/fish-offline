% Write netcdf of climate anoms

clear 
close all

%%
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

%% Isolate years of interest 1948-2015
fyr = 1948:2015;
mAMO = AMO((AMOyr>=1948 & AMOyr<=2015),:);
mNAO = NAO((NAOyr>=1948 & NAOyr<=2015),:);
mNino12 = Nino12((Nino12yr>=1948 & Nino12yr<=2015),:);
mNino34 = Nino34((Nino34yr>=1948 & Nino34yr<=2015),:);
mNino3 = Nino3((Nino3yr>=1948 & Nino3yr<=2015),:);
mNino4 = Nino4((Nino4yr>=1948 & Nino4yr<=2015),:);
mPDO = PDO((PDOyr>=1948 & PDOyr<=2015),:);
mSOI = SOI((SOIyr>=1948 & SOIyr<=2015),:);

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

tAO = AO((AOyr>=1948 & AOyr<=2015),:);
mAO = nan*ones(68,12);
mAO(3:end,:) = tAO;

tMEI = MEI((MEIyr>=1948 & MEIyr<=2015),:);
mMEI = nan*ones(68,12);
mMEI(32:end,:) = tMEI;

tNOI = NOI((NOIyr>=1948 & NOIyr<=2015),:);
mNOI = nan*ones(68,12);
mNOI(1:60,:) = tNOI;

%% Setup netcdf path to store to
file_tfb = [apath 'climate_anomalies_monthly.nc'];

%Use Netcdf4 classic
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

mo=1:12;
ni = length(fyr);
nj = length(mo);

ncidFB = netcdf.create(file_tfb,cmode);

%%

yr_dim = netcdf.defDim(ncidFB,'nyr',ni);
mo_dim = netcdf.defDim(ncidFB,'nmo',nj);

vidlon = netcdf.defVar(ncidFB,'yr','NC_DOUBLE',yr_dim);
netcdf.putAtt(ncidFB,vidlon,'long_name','year');
netcdf.putAtt(ncidFB,vidlon,'standard_name','year');
netcdf.putAtt(ncidFB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidFB,'mo','NC_DOUBLE',mo_dim);
netcdf.putAtt(ncidFB,vidlat,'long_name','month');
netcdf.putAtt(ncidFB,vidlat,'standard_name','month');
netcdf.putAtt(ncidFB,vidlat,'axis','Y');

%AO
vidAO = netcdf.defVar(ncidFB,'AO','NC_FLOAT',[yr_dim,mo_dim]);
netcdf.putAtt(ncidFB,vidAO,'long_name','Arctic Oscillation from NOAA NCEP');

%AMO
vidAMO = netcdf.defVar(ncidFB,'AMO','NC_FLOAT',[yr_dim,mo_dim]);
netcdf.putAtt(ncidFB,vidAMO,'long_name','AMO unsmoothed from the Kaplan SST V2');

%MEI
vidMEI = netcdf.defVar(ncidFB,'MEI','NC_FLOAT',[yr_dim,mo_dim]);
netcdf.putAtt(ncidFB,vidMEI,'long_name','Multivariate ENSO Index Version 2');

%NAO
vidNAO = netcdf.defVar(ncidFB,'NAO','NC_FLOAT',[yr_dim,mo_dim]);
netcdf.putAtt(ncidFB,vidNAO,'long_name','NAO from CRU');

%Nino34
vidNino34 = netcdf.defVar(ncidFB,'Nino34','NC_FLOAT',[yr_dim,mo_dim]);
netcdf.putAtt(ncidFB,vidNino34,'long_name','Nino34 5N-5S 170W-120W HadISST Anomaly');

%NOI
vidNOI = netcdf.defVar(ncidFB,'NOI','NC_FLOAT',[yr_dim,mo_dim]);
netcdf.putAtt(ncidFB,vidNOI,'long_name','Northern Oscillation Index');

%PDO
vidPDO = netcdf.defVar(ncidFB,'PDO','NC_FLOAT',[yr_dim,mo_dim]);
netcdf.putAtt(ncidFB,vidPDO,'long_name','PDO derived from OI.v2 SST fields');

%SOI
vidSOI = netcdf.defVar(ncidFB,'SOI','NC_FLOAT',[yr_dim,mo_dim]);
netcdf.putAtt(ncidFB,vidSOI,'long_name','Southern Oscillation Index from CRU');


varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidFB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidFB,varid,'MissingValue','NA');
netcdf.putAtt(ncidFB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidFB,varid,'institution','UCSD');

netcdf.endDef(ncidFB);

netcdf.putVar(ncidFB,vidlat,mo);
netcdf.putVar(ncidFB,vidlon,fyr);
netcdf.putVar(ncidFB,vidAO,mAO);
netcdf.putVar(ncidFB,vidAMO,mAMO);
netcdf.putVar(ncidFB,vidMEI,mMEI);
netcdf.putVar(ncidFB,vidNAO,mNAO);
netcdf.putVar(ncidFB,vidNino34,mNino34);
netcdf.putVar(ncidFB,vidNOI,mNOI);
netcdf.putVar(ncidFB,vidPDO,mPDO);
netcdf.putVar(ncidFB,vidSOI,mSOI);

netcdf.close(ncidFB);

%%
ncdisp(file_tfb)



