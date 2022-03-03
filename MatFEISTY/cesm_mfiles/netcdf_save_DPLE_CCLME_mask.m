% FEISTY output at CCLME locations

clear all
close all

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;
glme = double(lme_mask);
glme(glme<0) = nan;
tlme = glme(ID);

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];

%pick year
StartYr = 2015;
%loop over members
submem = 1:40;
mem=1;
Member = submem(mem);
harv = ['v14_Y' num2str(StartYr) '_M' num2str(Member) '_All_fish03_' ];

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

%% Netcdf OUTPUTS =================================================
t=time;
mo=t/12;
mo=mo+StartYr;
nt = 1:t;

allB = Bent.bio;

%% Reshape to lat,lon,yr
% Extract CC LME lat-lon only
vid = find(tlme==3);
AllB = NaN*ones(ni,nj,nt);

for z=1:nt
    Zb=NaN*ones(ni,nj);
    Zb(GRD.ID(vid))=allB(vid,z);
    AllB(:,:,z) = Zb;
end

%% Extract CC LME lat-lon only to reduce size
test = AllB(:,:,1);
cccol = find(~isnan(nansum(test)));
ccrow = find(~isnan(nansum(test,2)));

lat2 = TLAT(ccrow,cccol);
lon2 = TLONG(ccrow,cccol);
%cclme = glme(ccrow,cccol);
cclme = glme;
cclme(cclme~=3) = nan;
cclme(cclme==3) = ones;

%% test
plotminlat=21; %Set these bounds for your data
plotmaxlat=50;
plotminlon=-130;
plotmaxlon=-110;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%%
test = AllB(:,:,1);

figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat2,lon2,test); hold on;
hcb = colorbar('h');
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gcf,'renderer','painters')

%% save mask

save([cpath 'CESM_DPLE_CCLME_mask.mat'],...
    'lat2','lon2','cccol','ccrow',...
    'glme','tlme','cclme');

writematrix(cclme,[cpath 'CESM_DPLE_CCLME_mask.csv'],'Delimiter',',')
writematrix(TLAT,[cpath 'CESM_DPLE_TLAT.csv'],'Delimiter',',')
writematrix(TLONG,[cpath 'CESM_DPLE_TLONG.csv'],'Delimiter',',')

%% Setup netcdf path to store to
fname1 = 'CESM_DPLE_CCLME_mask.nc';

file_tfb = [cpath fname1];

LAT = TLAT;
LON = TLONG;

%Use Netcdf4 classic
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% 
ncidFB = netcdf.create(file_tfb,cmode);

lon_dim = netcdf.defDim(ncidFB,'lon',ni);
lat_dim = netcdf.defDim(ncidFB,'lat',nj);

vidlon = netcdf.defVar(ncidFB,'lon','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'standard_name','longitude');
netcdf.putAtt(ncidFB,vidlon,'axis','X');

vidlat = netcdf.defVar(ncidFB,'lat','NC_DOUBLE',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'standard_name','latitude');
netcdf.putAtt(ncidFB,vidlat,'axis','Y');

vidbioFB = netcdf.defVar(ncidFB,'cclme_mask','NC_FLOAT',[lon_dim,lat_dim]);
netcdf.putAtt(ncidFB,vidbioFB,'long_name','California Current LME mask');
netcdf.putAtt(ncidFB,vidbioFB,'missing value','NaN');

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidFB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidFB,varid,'_FillValue','NaN');
netcdf.putAtt(ncidFB,varid,'contact','C. Petrik <cpetrik@ucsd.edu>');
netcdf.putAtt(ncidFB,varid,'institution','UCSD');

netcdf.endDef(ncidFB);

%tfb = single(tfb);

netcdf.putVar(ncidFB,vidlat,LAT);
netcdf.putVar(ncidFB,vidlon,LON);
netcdf.putVar(ncidFB,vidbioFB,cclme);

netcdf.close(ncidFB);

%%
ncdisp(file_tfb)
