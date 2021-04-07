% Make GRD file for FEISTY input from CESM FOSI

clear all
close all

Cdir = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
fpath='/Volumes/MIP/GCM_DATA/CESM/4P4Z/';

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

%% Use bottom temp to set WID
ncid = netcdf.open([fpath 'g.e22a06.G1850ECOIAF_JRA_PHYS_DEV.TL319_g17.4p4z.004.FIESTY-forcing.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 14%1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
netcdf.close(ncid);

figure
pcolor(TEMP_bottom(:,:,20))

%%
tb = TEMP_bottom(:,:,1);
TID = find(~isnan(tb(:))); 
LID = find(~isnan(mask(:))); 
WID = find(HT(:)~=0);
NID = length(TID); %86096

eq1 = (WID==LID); 
sum(eq1)

eq2 = (WID==TID); 
sum(eq2)

%% Retain only water cells
ID = TID;
GRD.ID = ID;
GRD.N = NID;
GRD.LON = TLONG(ID);
GRD.LAT = TLAT(ID);
GRD.Z   = HT(ID) * 1e-2;     %from cm to m
GRD.area = TAREA(ID) * 1e-4; %from cm2 to m2
GRD.lmask = mask(ID);

%% Save needed variables
save([fpath 'gridspec_POP_gx1v6_4p4z.mat'],'HT','TLAT','TLONG','TAREA','mask',...
    'HTunits','HTlong_name','AREAunits','AREAlong_name');
save([fpath 'Data_grid_POP_gx1v6_4p4z.mat'],'GRD');
