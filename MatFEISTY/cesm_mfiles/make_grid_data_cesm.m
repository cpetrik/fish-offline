% Make GRD file for FEISTY input from CESM FOSI

clear all
close all

Cdir = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';

%% Depth, lat, lon
ncdisp([Cdir 'grid-data-POP_gx1v6.nc'])

% HT
% Size:       320x384
% Dimensions: nlon,nlat
% Datatype:   double
% Attributes:
% _FillValue = NaN
HTunits      = 'centimeter';
HTlong_name  = 'ocean depth at T points';
% note       = 'this field ignores overflows, which comprise isolated KMT pop-down points'

% TAREA
% Size:       320x384
% Dimensions: nlon,nlat
% Datatype:   double
% Attributes:
% _FillValue  = NaN
AREAunits       = 'cm^2';
AREAlong_name   = 'area of T cells';
% coordinates = 'TLONG TLAT'

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

mask = double(REGION_MASK);
mask(mask==0) = NaN;

figure
pcolor(mask)

%%
clatlim=[-90 90];
clonlim=[-180 180];

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,HT)
title('CESM depth')

%%
LID = find(~isnan(mask(:))); 
WID = find(HT(:)~=0);
NID = length(WID); %86212

eq1 = (WID==LID); 
sum(eq1)

%% Retain only water cells
ID = WID;
GRD.ID = ID;
GRD.N = NID;
GRD.LON = TLONG(ID);
GRD.LAT = TLAT(ID);
GRD.Z   = HT(ID) * 1e-2;     %from cm to m
GRD.area = TAREA(ID) * 1e-4; %from cm2 to m2
GRD.lmask = mask(ID);

%% Save needed variables
save([Cdir 'gridspec_POP_gx1v6.mat'],'HT','TLAT','TLONG','TAREA','mask',...
    'HTunits','HTlong_name','AREAunits','AREAlong_name');
save([Cdir 'Data_grid_POP_gx1v6.mat'],'GRD');
