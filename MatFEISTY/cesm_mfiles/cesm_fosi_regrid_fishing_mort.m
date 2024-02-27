% Regrid fishing effort to POP grid

clear
close all

%% Map
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;
WID = GRD.ID;
NID = GRD.N;

tlme = double(lme_mask);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
clatlim=[plotminlat plotmaxlat];
clonlim=[plotminlon plotmaxlon];

load coastlines;

%% ! Shift Longitudes !
TLON = TLONG;
TLON(TLONG>180) = TLONG(TLONG>180)-360;

%% ! --> use assessment estimate
alt1 = 'grid_mortality_guilds_v3'; %grid_mortality_guilds_v2, pristine_grid_mortality_guilds_v1
alt2 = '_v3'; %_v2, _v1_pristine
spath = ['/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/',alt1,'/'];
fpath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY_other/fishing_ms_ideas/fishing_effort_impl/',alt1,'/'];
load([fpath 'grid_mortality_all',alt2,'.mat'])

%% ! --> subset 1948-2015
yrall = 1841:2010;
yid = find(yrall>=1948);

fmortD = fmortD(:,yid);
fmortF = fmortF(:,yid);
fmortP = fmortP(:,yid);

%% 1/2 degree
lats = unique([LatD, LatF, LatP]);
lons = unique([LonD, LonF, LonP]);

%%
nt = length(yid);
fmD = zeros(NID,nt);
fmF = zeros(NID,nt);
fmP = zeros(NID,nt);

testD = fmortD(:,t);
F = scatteredInterpolant(LonD,LatD,testD(:));
zGrid = F(LON,LAT);

for t=1:nt
    clear testD testF testP

    testD = fmortD(:,t);
    F = scatteredInterpolant(LonD,LatD,testD(:));
    zGrid = F(LON,LAT);
    fmD(:,t) = testD(WID);

    testD = fmortD(:,t);
    F = scatteredInterpolant(LonD,LatD,testD(:));
    zGrid = F(LON,LAT);
    fmF(:,t) = testF(WID);

    testD = fmortD(:,t);
    F = scatteredInterpolant(LonD,LatD,testD(:));
    zGrid = F(LON,LAT);
    fmP(:,t) = testP(WID);
end

%%
clatlim=[-90 90];
clonlim=[-180 180];
load coastlines;
figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,testF)
caxis([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,testP)
caxis([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,testD)
caxis([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% temp scaling
tpath ='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/';
load([tpath 'gfdl-mom6-cobalt2_obsclim_mtemp_15arcmin_global_1961_2010.mat'])
vmtp = mtp(WID);
vmtb = mtb(WID);

%% scale with Fmsy and temp ! --> find temp means from FOSI !
% fm = F/Fmsy, need to mult by Fmsy ~= 0.3
%tsc = (exp(0.063*(temp-10.0));

fmF = 0.3 * fmF .* (exp(0.063*(vmtp-10.0)));
fmP = 0.3 * fmP .* (exp(0.063*(vmtp-10.0)));
fmD = 0.3 * fmD .* (exp(0.063*(vmtb-10.0)));

%%
fmD(isnan(fmD)) = 0.0;
fmF(isnan(fmF)) = 0.0;
fmP(isnan(fmP)) = 0.0;

fmD(fmD<0) = 0.0;
fmF(fmF<0) = 0.0;
fmP(fmP<0) = 0.0;

%% save
%year = 1961:2010;

% save([tpath 'gfdl-mom6-cobalt2_obsclim_15arcmin_fmort_ID_annual_1961_2010_tempSc',alt2,'.mat'],'year','WID',...
%     'fmD','fmF','fmP');
% save([fpath 'gfdl-mom6-cobalt2_obsclim_15arcmin_fmort_ID_annual_1961_2010_tempSc',alt2,'.mat'],'year','WID',...
%     'fmD','fmF','fmP');
% save([spath 'gfdl-mom6-cobalt2_obsclim_15arcmin_fmort_ID_annual_1961_2010_tempSc',alt2,'.mat'],'year','WID',...
%     'fmD','fmF','fmP');


