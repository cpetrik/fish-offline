% Make GRD file for FEISTY input from CESM FOSI

clear 
close all

Cdir = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/4P2Z/Hist/';
spath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/4P2Z/';


%% Depth, lat, lon
ncdisp([Cdir 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.TEMP_mean_150m.20000102-20141231.nc'])

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
ncid = netcdf.open([Cdir 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.TEMP_mean_150m.20000102-20141231.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

%nv = 1:nvars;
%vid = nv(nv~=18);
%for i = vid
for i = 1:(nvars-1)
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

%% Use bottom depth to set WID
WID = find(HT(:)~=0);
NID = length(WID); %86096


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
save([spath 'gridspec_POP_4p2z.mat'],'HT','TLAT','TLONG','TAREA','mask',...
    'HTunits','HTlong_name','AREAunits','AREAlong_name');
save([spath 'Data_grid_POP_4p2z.mat'],'GRD');
