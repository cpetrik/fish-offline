% Fishing effort additional years 2011-2015
% Check against FFmsy ending in 2010

clear
close all

%% 
spath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
    'FEISTY_other/fishing_ms_ideas/fishing_effort_ms/fishing_for_FEISTY/grid_mortality_guilds_v32/'];
load([spath 'grid_mortality_all_grid_mortality_guilds_v32.mat'])

% added yrs for geoengineering sims
gpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/CDR_Nuclear/geoengineering_fishing_processed/';
load([gpath 'all_enom_FFmsy_2011-2015.mat'])

%% 2010 and 2011
d10 = fmortD(:,end);
f10 = fmortF(:,end);
p10 = fmortP(:,end);

did = find(DYear==2011);
d11 = demFFmsy(did);
dLon11 = DLon(did);
dLat11 = DLat(did);

fid = find(FYear==2011);
f11 = forageFFmsy(fid);
fLon11 = FLon(fid);
fLat11 = FLat(fid);

pid = find(PYear==2011);
p11 = lgpelFFmsy(pid);
pLon11 = PLon(pid);
pLat11 = PLat(pid);

%% find same locations to compare
d10all = [LatD,LonD,d10]; 
f10all = [LatF,LonF,f10]; 
p10all = [LatP,LonP,p10]; 

d11all = [dLat11,dLon11,d11]; 
f11all = [fLat11,fLon11,f11]; 
p11all = [pLat11,pLon11,p11]; 

%% manually regrid? - THIS DID NOT WORK
lats = unique([DLat; FLat; PLat]);
lons = unique([DLon; FLon; PLon]);

%[gLAT,gLON] = meshgrid(min(lats):0.5:max(lats),min(lons):0.5:max(lons));
[gLAT,gLON] = meshgrid(lats,lons);

[ni,nj] = size(gLON);
d10G = nan*ones(ni,nj);
f10G = nan*ones(ni,nj);
p10G = nan*ones(ni,nj);
d11G = nan*ones(ni,nj);
f11G = nan*ones(ni,nj);
p11G = nan*ones(ni,nj);


for i=1:ni
    for j=1:nj
        fid10 = find(f10all(:,1)==gLAT(i) & f10all(:,2)==gLON(j));
        if(~isempty(fid10))
            f10G(i,j) = f10all(fid10,3);
        end
        pid10 = find(p10all(:,1)==gLAT(i) & p10all(:,2)==gLON(j));
        if(~isempty(pid10))
            p10G(i,j) = p10all(pid10,3);
        end
        did10 = find(d10all(:,1)==gLAT(i) & d10all(:,2)==gLON(j));
        if(~isempty(did10))
            d10G(i,j) = d10all(fid10,3);
        end

        fid11 = find(f11all(:,1)==gLAT(i) & f11all(:,2)==gLON(j));
        if(~isempty(fid11))
            f11G(i,j) = f11all(fid11,3);
        end
        pid11 = find(p11all(:,1)==gLAT(i) & p11all(:,2)==gLON(j));
        if(~isempty(pid11))
            p11G(i,j) = p11all(pid11,3);
        end
        did11 = find(d11all(:,1)==gLAT(i) & d11all(:,2)==gLON(j));
        if(~isempty(did11))
            d11G(i,j) = d11all(fid11,3);
        end
    end
end

%%
clatlim=[-90 90];
clonlim=[-180 180];
load coastlines;

figure
subplot('Position',[0.05 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,f10G)
clim([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.5 0.5 0.45 0.45])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,f11G)
clim([0 0.6])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.25 0.05 0.45 0.45])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(gLAT,gLON,f11G-f10G)
cmocean('balance')
clim([-1 1])
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);