% Calc powerspectra from COBALT ts for diff var
% Steps 3 & 5
% How get signif of trend? (Step 4)

% Steps
% 1. Calc seasonal climatology
% 2. Create anomaly time series by removing climatol
% 3. Calc Theil-Sen linear trend
% 4. Remove trend if significant
% 5. Calc spectral slopes

clear all
close all

%% FEISTY ---------------------------------------------------------
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

[nr,nc]=size(geolon_t);

% anomaly time series
load([fpath 'feisty_pi400_anom_ln.mat'],'yid','B_anom','F_anom','P_anom','D_anom',...
    'S_anom','M_anom','L_anom','all_anom')
nt = length(yid);

%% remove 1st 50 yrs (spinup)
% may need to remove transpose
xB  = B_anom(:,601:end)';
xF  = F_anom(:,601:end)';
xP  = P_anom(:,601:end)';
xD  = D_anom(:,601:end)';
xS  = S_anom(:,601:end)';
xM  = M_anom(:,601:end)';
xL  = L_anom(:,601:end)';
xall = all_anom(:,601:end)';

%% 
B = NaN*ones(size(xB,2),1);
F = NaN*ones(size(xF,2),1);
P = NaN*ones(size(xP,2),1);
D = NaN*ones(size(xD,2),1);
S = NaN*ones(size(xS,2),1);
M = NaN*ones(size(xM,2),1);
L = NaN*ones(size(xL,2),1);
All = NaN*ones(size(xall,2),1);

%% calc powerspec
for i = 1:size(xP, 2) 
    i
 
    xi = xB(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    B(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xF(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    F(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xP(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    P(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xD(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    D(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xS(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    S(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xM(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    M(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xL(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    L(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xall(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    All(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
end

%% put FEISTY back to the original map
nidx = GRD.ID;

B_ps = NaN*ones(nr,nc);
F_ps = NaN*ones(nr,nc);
P_ps = NaN*ones(nr,nc);
D_ps = NaN*ones(nr,nc);
S_ps = NaN*ones(nr,nc);
M_ps = NaN*ones(nr,nc);
L_ps = NaN*ones(nr,nc);
All_ps = NaN*ones(nr,nc);

B_ps(nidx)  = B;
F_ps(nidx)  = F;
P_ps(nidx)  = P;
D_ps(nidx)  = D;
S_ps(nidx)  = S;
M_ps(nidx)  = M;
L_ps(nidx)  = L;
All_ps(nidx) = All;

%%
save([fpath 'powerspec_feisty_types_pi400_ln.mat'],...
    'geolat_t','geolon_t','yid','GRD',...
    'B_ps','F_ps','P_ps','D_ps',...
    'S_ps','M_ps','L_ps','All_ps');

%% Maps
plotminlat=-90; 
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

cmapB = cmocean('balance');
cmapB = flipud(cmapB);

lat = double(geolat_t);
lon = double(geolon_t);

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
ppath = [pp cfile '/PreIndust400/'];

%% 8 plot for fish groups
f4 = figure('Units','inches','Position',[1 3 6.5 8]);
%f4.Units = 'inches';

%A 
subplot('Position',[0.025 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,F_ps)
colormap(cmapB)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-6 6]);
set(gcf,'renderer','painters')
text(0,1.75,'Forage','HorizontalAlignment','center')

%B 
subplot('Position',[0.025 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,P_ps)
colormap(cmapB)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-6 6]);
set(gcf,'renderer','painters')
text(0,1.75,'Large Pelagic','HorizontalAlignment','center')

%C 
subplot('Position',[0.025 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,D_ps)
colormap(cmapB)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-6 6]);
set(gcf,'renderer','painters')
text(0,1.75,'Demersal','HorizontalAlignment','center')

%D 
subplot('Position',[0.025 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,B_ps)
colormap(cmapB)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-6 6]);
set(gcf,'renderer','painters')
text(0,1.75,'Benthos','HorizontalAlignment','center')

%E 
subplot('Position',[0.475 0.75 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,S_ps)
colormap(cmapB)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-6 6]);
set(gcf,'renderer','painters')
text(0,1.75,'Small','HorizontalAlignment','center')

%F 
subplot('Position',[0.475 0.5 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,M_ps)
colormap(cmapB)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-6 6]);
colorbar('Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'Medium','HorizontalAlignment','center')

%G 
subplot('Position',[0.475 0.25 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,L_ps)
colormap(cmapB)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-6 6]);
set(gcf,'renderer','painters')
text(0,1.75,'Large','HorizontalAlignment','center')

%H 
subplot('Position',[0.475 0.0 0.4 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,All_ps)
colormap(cmapB)
load coastlines;                     %decent looking coastlines
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-6 6]);
set(gcf,'renderer','painters')
text(0,1.75,'All fish','HorizontalAlignment','center')

print('-dpng',[ppath 'Pre400_fish_groups_global_ps_ln.png'])

%% Histograms
edges = -6:0.25:1;
 
figure(6)
subplot(3,3,1)
histogram(F,edges)
title('F')

subplot(3,3,2)
histogram(P,edges)
title('P')

subplot(3,3,3)
histogram(D,edges)
title('D')

subplot(3,3,4)
histogram(S,edges)
title('S')

subplot(3,3,5)
histogram(M,edges)
title('M')
subplot(3,3,6)
histogram(L,edges)
title('L')

subplot(3,3,7)
histogram(All,edges)
title('All')

subplot(3,3,8)
histogram(B,edges)
title('B')
print('-dpng',[ppath 'Pre400_histogram_fish_groups_ps_ln.png'])


