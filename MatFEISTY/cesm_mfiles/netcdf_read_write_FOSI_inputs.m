% Write one netcdf file of FOSI 
% FEISTY inputs

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
spath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

%% Inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TLAT','TLONG','TAREA','time','yr',...
    'TEMP_150m','TEMP_150m_units','TEMP_150m_name',...
    'TEMP_bottom','TEMP_bottom_units','TEMP_bottom_name',...
    'POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units','POC_FLUX_IN_bottom_name',...
    'zooC_150m_units','zoo_loss_150m_units');
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo_totloss_allphytoC.mat'],...
    'LzooC_150m','Lzoo_loss_150m','fracL');
load([fpath 'gridspec_POP_gx1v6_noSeas.mat'],'mask');

%% 
LzooC_150m_units = zooC_150m_units;
Lzoo_loss_150m_units = zoo_loss_150m_units;
LzooC_150m_name = 'Mesozooplankton carbon biomass integrated in top 150m';
Lzoo_loss_150m_name = 'Mesozooplankton mortality rate integrated in top 150m';

LzooC_150m(LzooC_150m<0) = 0.0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0.0;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom<0) = 0.0;

LzooC_150m = single(LzooC_150m);
Lzoo_loss_150m = single(Lzoo_loss_150m);
fracL = single(fracL);

%% netcdf write
% nans to a large number
TEMP_150m(isnan(TEMP_150m)) = missing_value;
TEMP_bottom(isnan(TEMP_bottom)) = missing_value;
POC_FLUX_IN_bottom(isnan(POC_FLUX_IN_bottom)) = missing_value;
LzooC_150m(isnan(LzooC_150m)) = missing_value;
Lzoo_loss_150m(isnan(Lzoo_loss_150m)) = missing_value;
fracL(isnan(fracL)) = missing_value;

%% Quick look
pb = TEMP_150m(:,:,50);
db = TEMP_bottom(:,:,50);
cb = LzooC_150m(:,:,50);
bp30 = POC_FLUX_IN_bottom(:,:,50);

figure(1)
pcolor((pb'))
shading flat
colormap('jet')
colorbar
caxis([-2 25])
title('TP')

figure(2)
pcolor((db'))
shading flat
colormap('jet')
colorbar
caxis([-2 25])
title('TB')

figure(3)
pcolor(log10(cb'))
shading flat
colormap('jet')
colorbar
caxis([2 4])
title('LZ')

figure(4)
pcolor(log10(bp30'))
shading flat
colormap('jet')
colorbar
caxis([-5 -2])
title('POC')

%%
close all

%% Setup netcdf path to store to
file_tfb = [fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FEISTY_forcing_complete.nc'];

[ni,nj,nt] = size(fracL);

t=1:nt;
mo=t/12;
yr=mo+1948;

%%
LAT = single(TLAT);
LON = single(TLONG);

%Use Netcdf4 classic
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% 
ncidFB = netcdf.create(file_tfb,cmode);

time_dim = netcdf.defDim(ncidFB,'time',nt);
lon_dim = netcdf.defDim(ncidFB,'nlon',ni);
lat_dim = netcdf.defDim(ncidFB,'nlat',nj);

vidtFB = netcdf.defVar(ncidFB,'time','NC_DOUBLE',time_dim);
netcdf.putAtt(ncidFB,vidtFB,'long_name','time');
netcdf.putAtt(ncidFB,vidtFB,'standard_name','time');
netcdf.putAtt(ncidFB,vidtFB,'calendar','365_day');
netcdf.putAtt(ncidFB,vidtFB,'axis','T');
netcdf.putAtt(ncidFB,vidtFB,'units','year' );

vidlon = netcdf.defVar(ncidFB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidFB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'axis','Y');

vidTP = netcdf.defVar(ncidFB,'TEMP_150m','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidFB,vidTP,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidTP,'long_name',[TEMP_150m_name ' mean top 150m']);
netcdf.putAtt(ncidFB,vidTP,'units',TEMP_150m_units);
netcdf.defVarFill(ncidFB,vidTP,false,missing_value);
netcdf.putAtt(ncidFB,vidTP,'missing value',missing_value);

vidTB = netcdf.defVar(ncidFB,'TEMP_bottom','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidFB,vidTB,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidTB,'long_name',TEMP_bottom_name);
netcdf.putAtt(ncidFB,vidTB,'units',TEMP_bottom_units);
netcdf.defVarFill(ncidFB,vidTB,false,missing_value);
netcdf.putAtt(ncidFB,vidTB,'missing value',missing_value);

vidPOC = netcdf.defVar(ncidFB,'POC_FLUX_IN_bottom','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidFB,vidPOC,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidPOC,'long_name',POC_FLUX_IN_bottom_name);
netcdf.putAtt(ncidFB,vidPOC,'units',POC_FLUX_IN_bottom_units);
netcdf.defVarFill(ncidFB,vidPOC,false,missing_value);
netcdf.putAtt(ncidFB,vidPOC,'missing value',missing_value);

vidMZ = netcdf.defVar(ncidFB,'LzooC_150m','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidFB,vidMZ,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidMZ,'long_name',LzooC_150m_name);
netcdf.putAtt(ncidFB,vidMZ,'units',LzooC_150m_units);
netcdf.defVarFill(ncidFB,vidMZ,false,missing_value);
netcdf.putAtt(ncidFB,vidMZ,'missing value',missing_value);

vidMZL = netcdf.defVar(ncidFB,'Lzoo_loss_150m','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidFB,vidMZL,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidMZL,'long_name',Lzoo_loss_150m_name);
netcdf.putAtt(ncidFB,vidMZL,'units',Lzoo_loss_150m_units);
netcdf.defVarFill(ncidFB,vidMZL,false,missing_value);
netcdf.putAtt(ncidFB,vidMZL,'missing value',missing_value);

vidFrac = netcdf.defVar(ncidFB,'fracL','NC_FLOAT',[lon_dim,lat_dim,time_dim]);
netcdf.defVarChunking(ncidFB,vidFrac,'CHUNKED',[10, 10, 1]);
netcdf.putAtt(ncidFB,vidFrac,'long_name','Fraction of large plankton');
netcdf.putAtt(ncidFB,vidFrac,'units','');
netcdf.defVarFill(ncidFB,vidFrac,false,missing_value);
netcdf.putAtt(ncidFB,vidFrac,'missing value',missing_value);

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidFB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidFB,varid,'_FillValue',missing_value);
netcdf.putAtt(ncidFB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidFB,varid,'institution','UCSD');

netcdf.endDef(ncidFB);

netcdf.putVar(ncidFB,vidlat,LAT);
netcdf.putVar(ncidFB,vidlon,LON);
netcdf.putVar(ncidFB,vidtFB,yr);
netcdf.putVar(ncidFB,vidTP,TEMP_150m);
netcdf.putVar(ncidFB,vidTB,TEMP_bottom);
netcdf.putVar(ncidFB,vidPOC,POC_FLUX_IN_bottom);
netcdf.putVar(ncidFB,vidMZ,LzooC_150m);
netcdf.putVar(ncidFB,vidMZL,Lzoo_loss_150m);
netcdf.putVar(ncidFB,vidFrac,fracL);

netcdf.close(ncidFB);

%%
ncdisp(file_tfb)

