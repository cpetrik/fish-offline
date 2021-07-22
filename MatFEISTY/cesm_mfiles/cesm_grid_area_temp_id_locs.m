%% Area of locations

clear all
close all

ppath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/4P4Z/';
cpath = '/Volumes/MIP/GCM_DATA/CESM/4P4Z/';

load([cpath 'gridspec_POP_gx1v6_4p4z.mat']);
load([cpath 'Data_grid_POP_gx1v6_4p4z.mat']);

%% fix lon shift

geolat = TLAT;
geolon = TLONG;
test = geolon-360;
id=find(test<-180);
test(id)=test(id)+360;
geolon = test;
figure
pcolor(TLONG);
figure
pcolor(geolon)

figure
pcolor(TLAT);

%%
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,TLONG)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,geolon,geolon)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,geolon)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%%
geolat = GRD.LAT;
geolon = GRD.LON;
test = geolon-360;
id=find(test<-180);
test(id)=test(id)+360;
geolon = test;

%% Grid

t2      = GRD.ID;
t2(:,2) = geolon;
t2(:,3) = geolat;
t2(:,4) = GRD.area;
t2(:,5) = GRD.Z;


%% Georges Bank (Northeast Peak - 41 deg 43.92' N x 66 deg 32.18' W
%Southern Flank - 40 deg 57.95' N x 67 deg 18.91' W)
lon=find(t2(:,2)<=-67 & t2(:,2)>=-68);
lat=find(t2(:,3)<=41.5 & t2(:,3)>=41);
gid=intersect(lon,lat);

% W Scotian Shelf (42.4910,-65.4670)
lon=find(t2(:,2)<=-65 & t2(:,2)>=-66);
lat=find(t2(:,3)<=42.5 & t2(:,3)>=42);
wssid=intersect(lon,lat);

% C Scotian Shelf (43.5120,-62.4780)
lon=find(t2(:,2)<=-62 & t2(:,2)>=-63);
lat=find(t2(:,3)<=44 & t2(:,3)>=43.5);
cssid=intersect(lon,lat);

% E Scotian Shelf (45.0730,-59.1160)
lon=find(t2(:,2)<=-59 & t2(:,2)>=-60);
lat=find(t2(:,3)<=46 & t2(:,3)>=45);
essid=intersect(lon,lat);

% Greenland Sea (77.8710° N, 5.6501° W)
lon=find(t2(:,2)<=-5.5 & t2(:,2)>=-6);
lat=find(t2(:,3)<=78 & t2(:,3)>=77.5);
gsid=intersect(lon,lat);

% North Sea
lon=find(t2(:,2)<=4 & t2(:,2)>=3.25);
lat=find(t2(:,3)<=57 & t2(:,3)>=56.5);
nid=intersect(lon,lat);

% Norwegian Sea (68.8774° N, 3.1397° E)
lon=find(t2(:,2)<=4 & t2(:,2)>=3);
lat=find(t2(:,3)<=69 & t2(:,3)>=68);
nwid=intersect(lon,lat);

% Barents Sea (74.9884° N, 37.1064° E)
lon=find(t2(:,2)<=38 & t2(:,2)>=37);
lat=find(t2(:,3)<=75.5 & t2(:,3)>=74.5);
bsid=intersect(lon,lat);

% Eastern Bering Sea (M2 Southeastern Bering Sea (56.87°N, -164.06°W))
lon=find(t2(:,2)<=-164 & t2(:,2)>=-165);
lat=find(t2(:,3)<=57 & t2(:,3)>=56.5);
eid=intersect(lon,lat);

% Subarctic Pacific Gyre (Ocean Station Papa)
lon=find(t2(:,2)<=-145 & t2(:,2)>=-146);
lat=find(t2(:,3)<=51 & t2(:,3)>=50);
pid=intersect(lon,lat);

% HOT (22° 45'N, 158° 00'W)
lon=find(t2(:,2)<=-157.35 & t2(:,2)>=-158);
lat=find(t2(:,3)<=23 & t2(:,3)>=22.45);
hid=intersect(lon,lat);

% BATS (31 50'N, 64 10'W)
lon=find(t2(:,2)<=-63.95 & t2(:,2)>=-65);
lat=find(t2(:,3)<=31.9 & t2(:,3)>=31.1);
bid=intersect(lon,lat);

% Eastern Equatorial Pacific
lon=find(t2(:,2)<=-110.3 & t2(:,2)>=-111);
lat=find(t2(:,3)<=5.4 & t2(:,3)>=5);
qid=intersect(lon,lat);

% Peru Upwelling
lon=find(t2(:,2)<=-79 & t2(:,2)>=-80);
lat=find(t2(:,3)<=-12.7 & t2(:,3)>=-13);
uid=intersect(lon,lat);

%Subpolar W Pac station K2: 47oN, 160oE
lon=find(t2(:,2)<=161 & t2(:,2)>=160);
lat=find(t2(:,3)<=47.5 & t2(:,3)>=47);
kid=intersect(lon,lat);

%Subtropical W Pac station S1: 30oN, 145oE
lon=find(t2(:,2)<=146 & t2(:,2)>=145);
lat=find(t2(:,3)<=30.5 & t2(:,3)>=30);
sid=intersect(lon,lat);


%% Save
names={'Georges Bank','W Scotian Shelf','C Scotian Shelf','E Scotian Shelf',...
    'Greenland Sea','North Sea','Norwegian Sea','Barents Sea',...
    'Ocean Station Papa','Eastern Bering Sea','K2','S1',...
    'HOT','BATS','Eastern Equatorial Pacific','Peru Upwell'};
abbrev = {'GB','WSS','CSS','ESS',...
    'GS','NS','NwS','BS',...
    'OSP','EBS','K2','S1',...
    'HOT','BATS','EEP','PUp'};

ids(1,1)=gid;
ids(2,1)=wssid;
ids(3,1)=cssid;
ids(4,1)=essid;
ids(5,1)=gsid;
ids(6,1)=nid;
ids(7,1)=nwid;
ids(8,1)=bsid;
ids(9,1)=pid;
ids(10,1)=eid;
ids(11,1)=kid;
ids(12,1)=sid;
ids(13,1)=hid;
ids(14,1)=bid;
ids(15,1)=qid;
ids(16,1)=uid;

lons = t2(ids,2);
lats = t2(ids,3);
area = t2(ids,4);
depth = t2(ids,5);

%%
T=table(names',ids,lons,lats,area,depth,...
    'VariableNames',{'Location','ID','Lon','Lat','Area','Depth'});
writetable(T,[cpath 'cesm_4p4z_grid_id_locs_area_dep.csv'],'Delimiter',',');
save([cpath 'cesm_4p4z_grid_id_locs_area_dep.mat'],'T','depth','ids','abbrev');

%%
% plot info
[ni,nj]=size(lon);
geolon_t = TLONG;
geolat_t = TLAT;
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; 
load coastlines;     

%%
figure(1)
%m_proj('miller','lat',82);
%m_coast('patch',[.5 .5 .5],'edgecolor','none');
%m_grid;
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
for i =1:16
    %m_text(GRD.LON(ids(i)),GRD.LAT(ids(i)),abbrev{i},'Color','black','HorizontalAlignment','center');
    textm(GRD.LAT(ids(i)),GRD.LON(ids(i)),abbrev{i},'Color','black','HorizontalAlignment','center');
end
%print('-dpng',[ppath 'POEM_Clim_locations.png'])





