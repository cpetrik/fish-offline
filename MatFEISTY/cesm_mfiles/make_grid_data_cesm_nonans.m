% Find grid cells with non-nan values of all variables
% FOSI POP grid

clear all
close all

Cdir = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';

%% Depth, lat, lon
ncdisp([Cdir 'grid-data-POP_gx1v6.nc'])

HTunits      = 'centimeter';
HTlong_name  = 'ocean depth at T points';
AREAunits       = 'cm^2';
AREAlong_name   = 'area of T cells';

%%
ncid = netcdf.open([Cdir 'grid-data-POP_gx1v6.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

nv = 1:nvars;
vid = nv(nv~=18);
for i = vid
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
end
netcdf.close(ncid);

%% check orientation
figure
pcolor(HT)

figure
pcolor(TLAT)

figure
pcolor(TLONG)

%% create land mask
pcolor(double(REGION_MASK))

regions = unique(double(REGION_MASK));
%remove inland seas <0

mask = double(REGION_MASK);
mask(mask<=0) = NaN;

figure
pcolor(mask)

%% compare to manual
field_masked = nan(size(KMT));
field_masked(KMT(:) > 0) = 1;
field_masked(KMT(:) <= 0) = nan;

figure
pcolor(field_masked) %has inland seas

%%
clatlim=[-90 90];
clonlim=[-180 180];

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,HT)
title('CESM depth')

%%
MID = find(~isnan(mask(:)));    %85813
NID = length(MID);              %85813

WID = find(HT(:)~=0);           %86212
KID = find(KMT(:)>0);           %86212

eq1 = (WID==MID); 
sum(eq1)

eq2 = (KID==MID); 
sum(eq2)

%% compare this against vars
load([Cdir 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'time','yr');
load([Cdir 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo_totloss_allphytoC.mat'],...
    'LzooC_150m','Lzoo_loss_150m');

% nans & zeros
TEMP_150m = double(TEMP_150m);
TEMP_bottom = double(TEMP_bottom);
POC_FLUX_IN_bottom = double(POC_FLUX_IN_bottom);

TEMP_bottom(TEMP_bottom >= 9.9e+36) = nan;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;

LzooC_150m(LzooC_150m<0) = 0.0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0.0;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom<0) = 0.0;

Tp  = squeeze(TEMP_150m(:,:,1));
Tb  = squeeze(TEMP_bottom(:,:,5));
Zm  = squeeze(LzooC_150m(:,:,10));
dZm = squeeze(Lzoo_loss_150m(:,:,15));
det = squeeze(POC_FLUX_IN_bottom(:,:,20));

TID = find(~isnan(Tp(:)));
BID = find(~isnan(Tb(:)));
ZID = find(~isnan(Zm(:)));
LID = find(~isnan(dZm(:)));
DID = find(~isnan(det(:)));

whos TID BID ZID LID DID

% eq1 = (WID==TID); 
% sum(eq1)
% eq2 = (WID==BID); 
% sum(eq2)
% eq3 = (WID==ZID); 
% sum(eq3)
% 
% eq4 = (DID==BID); 
% sum(eq4)
% eq5 = (LID==ZID); 
% sum(eq5)

%% Retain only water cells
ID = MID;
GRD.ID = ID;
GRD.N = NID;
GRD.LON = TLONG(ID);
GRD.LAT = TLAT(ID);
GRD.Z   = HT(ID) * 1e-2;     %from cm to m
GRD.area = TAREA(ID) * 1e-4; %from cm2 to m2
GRD.lmask = mask(ID);

%% Save needed variables
save([Cdir 'gridspec_POP_gx1v6_noSeas.mat'],'HT','TLAT','TLONG','TAREA','mask',...
    'HTunits','HTlong_name','AREAunits','AREAlong_name');
save([Cdir 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');

