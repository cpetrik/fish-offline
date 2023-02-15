% FEISTY output at all locations

clear
close all

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/DPLE/'];

%pick year
StartYr = 1954;
%pick member
Member = 1;

harv = ['v15_Y' num2str(StartYr) '_M' num2str(Member) '_All_fish03_' ];

%% Read netcdfs ----------------------------------------------------------
% SP
ncid = netcdf.open([fpath 'DPLE_' harv 'sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.bio = biomass;
SP.prod = prod;
clear biomass prod

% SF
ncid = netcdf.open([fpath 'DPLE_' harv 'sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass;
SF.prod = prod;
clear biomass prod

% SD
ncid = netcdf.open([fpath 'DPLE_' harv 'sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
SD.prod = prod;
clear biomass prod

% MP
ncid = netcdf.open([fpath 'DPLE_' harv 'med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.prod = prod;
MP.yield = yield;
clear yield biomass prod

% MF
ncid = netcdf.open([fpath 'DPLE_' harv 'med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.prod = prod;
MF.yield = yield;
clear yield biomass prod

% MD
ncid = netcdf.open([fpath 'DPLE_' harv 'med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.prod = prod;
MD.yield = yield;
clear yield biomass prod

% LP
ncid = netcdf.open([fpath 'DPLE_' harv 'lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.prod = prod;
LP.yield = yield;
clear yield biomass prod

% LD
ncid = netcdf.open([fpath 'DPLE_' harv 'lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.prod = prod;
LD.yield = yield;
clear yield biomass prod

% Benthic material
ncid = netcdf.open([fpath 'DPLE_' harv 'bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
clear biomass

%% Calc groups
AllF = SF.bio + MF.bio;
AllP = SP.bio + MP.bio + LP.bio;
AllD = SD.bio + MD.bio + LD.bio;
AllS = SF.bio + SP.bio + SD.bio;
AllM = MF.bio + MP.bio + MD.bio;
AllL = LP.bio + LD.bio;
All = AllF + AllP + AllD;

%% Put on grid ------------------------------------------------------------
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
% 2-D POP grid
load([cpath 'gridspec_POP_gx1v6_noSeas.mat'],'TLAT','TLONG');
% GRD is matrix of vector location in 2-D grid; also has Lat, Lon, Depth
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');

[ni,nj] = size(TLONG);

[nid,nt] = size(LD.bio);

tF=NaN*ones(ni,nj);
tP=NaN*ones(ni,nj);
tD=NaN*ones(ni,nj);
tS=NaN*ones(ni,nj);
tM=NaN*ones(ni,nj);
tL=NaN*ones(ni,nj);
tA=NaN*ones(ni,nj);
tB=NaN*ones(ni,nj);

F=NaN*ones(ni,nj,nt);
P=NaN*ones(ni,nj,nt);
D=NaN*ones(ni,nj,nt);
S=NaN*ones(ni,nj,nt);
M=NaN*ones(ni,nj,nt);
L=NaN*ones(ni,nj,nt);
A=NaN*ones(ni,nj,nt);
B=NaN*ones(ni,nj,nt);

for t=1:length(nt)

    tF(GRD.ID) = AllF(:,t);
    tP(GRD.ID) = AllP(:,t);
    tD(GRD.ID) = AllD(:,t);
    tS(GRD.ID) = AllS(:,t);
    tM(GRD.ID) = AllM(:,t);
    tL(GRD.ID) = AllL(:,t);
    tA(GRD.ID) = All(:,t);
    tB(GRD.ID) = Bent.bio(:,t);

    F(:,:,t) = tF;
    P(:,:,t) = tP;
    D(:,:,t) = tD;
    S(:,:,t) = tS;
    M(:,:,t) = tM;
    L(:,:,t) = tL;
    A(:,:,t) = tA;
    B(:,:,t) = tB;

end











