% Visualize output of FEISTY
% CESM Spinup first year of FOSI
% Time series plots and maps

clear all
close all

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
mod = 'quad_v7_All_fish03';
%mod = 'All_fish03';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([fpath 'Means_Spinup_',mod,'_lowbiom.mat']);

% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI//';
load([cpath 'gridspec_POP_gx1v6.mat']);
load([cpath 'Data_grid_POP_gx1v6.mat']);

[ni,nj]=size(TLONG);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
% plotminlon=-180;
% plotmaxlon=180;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;     

%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black

set(groot,'defaultAxesColorOrder',cm10);

%cmBP50=cbrewer('seq','BuPu',50,'PCHIP');
cmBP50=cbrewer('seq','YlGnBu',10,'PCHIP');

%% Locations of low biom
lowF = zeros(size(TLONG));
lowF(fid) = ones;

lowP = zeros(size(TLONG));
lowP(pid) = ones;

Fid = zeros(size(TLONG));
Fid(fid) = fid;

Pid = zeros(size(TLONG));
Pid(pid) = pid;

%%
figure(1)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,lowF)
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('low F')

% all D
% subplot('Position',[0 0 0.5 0.5])


% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,lowP)
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('low P')

% All
% subplot('Position',[0.5 0 0.5 0.5])

stamp('')
print('-dpng',[ppath 'Spinup_',mod,'_global_FP_lowbiom.png'])

%% grid cell id > 10e4
figure(2)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,Fid)
colormap('jet')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
title('low F')

% all D
% subplot('Position',[0 0 0.5 0.5])

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,Pid)
colormap('jet')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%caxis([0 1]);
colorbar
set(gcf,'renderer','painters')
title('low P')

%% Time series of those cells
harv = 'quad_v6_All_fish03';

% SP
ncid = netcdf.open([fpath 'Spinup_' harv '_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%
[ni,nt] = size(biomass);

SP.bio = biomass;
clear biomass

% SF
ncid = netcdf.open([fpath 'Spinup_' harv '_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass 

% MP
ncid = netcdf.open([fpath 'Spinup_' harv '_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
clear biomass

% MF
ncid = netcdf.open([fpath 'Spinup_' harv '_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
clear biomass

% LP
ncid = netcdf.open([fpath 'Spinup_' harv '_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
clear biomass

%% Take means for my own visualization
% grid cell id > 10e4
fid2 = fid(fid>10e4);
pid2 = pid(pid>10e4);
[fx,fid3] = intersect(GRD.ID,fid2);
[px,pid3] = intersect(GRD.ID,pid2);

%%
sp_tmean=mean(SP.bio(pid3,:),1);
sf_tmean=mean(SF.bio(fid3,:),1);
mp_tmean=mean(MP.bio(pid3,:),1);
mf_tmean=mean(MF.bio(fid3,:),1);
lp_tmean=mean(LP.bio(pid3,:),1);

%% Plots in time
t = 1:length(sf_tmean); %time;
y = t/12;

F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;

%% Types together

figure(2)
plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
legend('F','P')
legend('location','east')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('Spinup')
stamp(mod)
%print('-dpng',[ppath 'Spinup_',mod,'_all_types.png'])

%%
figure(3)
subplot(2,1,1)
plot(y,log10(SF.bio(fid3(20:30),:))); hold on;
xlim([y(1) y(end)])
ylim([-8 0])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('SF')

subplot(2,1,2)
plot(y,log10(SP.bio(pid3(20:30),:))); hold on;
xlim([y(1) y(end)])
ylim([-8 0])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('SP')
stamp(mod)


figure(4)
subplot(2,1,1)
plot(y,log10(MF.bio(fid3(20:30),:))); hold on;
xlim([y(1) y(end)])
ylim([-8 0])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('MF')

subplot(2,1,2)
plot(y,log10(MP.bio(pid3(20:30),:))); hold on;
xlim([y(1) y(end)])
ylim([-8 0])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('MP')
stamp(mod)
 

