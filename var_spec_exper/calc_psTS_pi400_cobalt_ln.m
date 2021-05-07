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
%clear tp_clim tb_clim mz_clim lz_clim mhp_clim lhp_clim det_clim
close all

%% COBALT --------------------------------------------------------
spath='/Volumes/MIP/GCM_DATA/ESM2M_PI/';

% anomaly time series
load([spath 'cobalt_pi400_temp_anom.mat']);
load([spath 'cobalt_pi400_det_anom_ln.mat']);
load([spath 'cobalt_pi400_zoo_anom_ln.mat']);
load([spath 'cobalt_pi400_hploss_anom_ln.mat']);
[nr,nc] = size(lon);
nt = length(yid);

%% reshape to save only land-free cells
% TP
nidxmat = reshape(1:nr*nc, nr, nc);
xTP = reshape(tp_anom, nr*nc, nt);
nnan = reshape(sum(isnan(xTP')), nr, nc);
nidx = nidxmat(nnan < 1);
xTP = xTP(nidx, :)';
% TB
xTB = reshape(tb_anom, nr*nc, nt);
xTB = xTB(nidx, :)';
% DET
xDet = reshape(det_anom, nr*nc, nt);
xDet = xDet(nidx, :)';
% ZOOMED_INT100
xMZ = reshape(mz_anom, nr*nc, nt);
xMZ = xMZ(nidx, :)';
% ZOOLRG_INT100
xLZ = reshape(lz_anom, nr*nc, nt);
xLZ = xLZ(nidx, :)';
% HPLOSS MED
xHPM = reshape(hpmz_anom, nr*nc, nt);
xHPM = xHPM(nidx, :)';
% HPLOSS LRG
xHPL = reshape(hplz_anom, nr*nc, nt);
xHPL = xHPL(nidx, :)';

%% Exclude 1st 50 yrs (spinup for FEISTY)
xTP = xTP(601:end,:);
xTB = xTB(601:end,:);
xDet = xDet(601:end,:);
xMZ = xMZ(601:end,:);
xLZ = xLZ(601:end,:);
xHPM = xHPM(601:end,:);
xHPL = xHPL(601:end,:);

clear tp_anom tb_anom mz_anom lz_anom hpmz_anom hpmz_anom det_anom

%% 
tp = NaN*ones(size(xTP,2),1);
tb = NaN*ones(size(xTB,2),1);
det = NaN*ones(size(xDet,2),1);
zm = NaN*ones(size(xMZ,2),1);
zl = NaN*ones(size(xLZ,2),1);
hpmz = NaN*ones(size(xHPM,2),1);
hplz = NaN*ones(size(xHPL,2),1);

%% calc powerspec
for i = 1:size(xTP, 2) 
    i
 
    xi = xTP(:, i);
    inan = ~isnan(xi);
    %R = log(xi(inan)); %Already log-transformed
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    %b = Theil_Sen_Regress(t,R); %279.987615 seconds.
    %int = nanmedian(R-b.*t);
    %tH = b*t + int;
    data = [t R];
    [m,b] = TheilSen(data); %0.651784 seconds.
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    tp(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xTB(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    tb(i) = b1; 
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xDet(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    det(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xMZ(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data);
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    zm(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xLZ(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    zl(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xHPM(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    hpmz(i) = b1; 
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xHPL(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    hplz(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
end


%% put COBALT back to the original map
tp_ps = nan(nr*nc,1);
tp_ps(nidx) = tp;
tp_ps = reshape(tp_ps, nr, nc);

tb_ps = nan(nr*nc,1);
tb_ps(nidx) = tb;
tb_ps = reshape(tb_ps, nr, nc);

det_ps = nan(nr*nc,1);
det_ps(nidx) = det;
det_ps = reshape(det_ps, nr, nc);

zm_ps = nan(nr*nc,1);
zm_ps(nidx) = zm;
zm_ps = reshape(zm_ps, nr, nc);

zl_ps = nan(nr*nc,1);
zl_ps(nidx) = zl;
zl_ps = reshape(zl_ps, nr, nc);

hpmz_ps = nan(nr*nc,1);
hpmz_ps(nidx) = hpmz;
hpmz_ps = reshape(hpmz_ps, nr, nc);

hplz_ps = nan(nr*nc,1);
hplz_ps(nidx) = hplz;
hplz_ps = reshape(hplz_ps, nr, nc);

%%
% tp_ps = tp;
% tb_ps = tb;
% det_ps = det;
% zm_ps = zm;
% zl_ps = zl;
% hpmz_ps = hpmz;
% hplz_ps = hplz;
% 
% save([spath 'powerspec_cobalt_pi400_ln_141.mat'],...
%     'tp_ps','tb_ps','det_ps',...
%     'zl_ps','zm_ps','hpmz_ps','hplz_ps');

save([spath 'powerspec_cobalt_pi400_ln.mat'],...
    'tp_ps','tb_ps','det_ps',...
    'zl_ps','zm_ps','hpmz_ps','hplz_ps');

%% Maps
plotminlat=-90; 
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

cmapB = cmocean('balance');
cmapB = flipud(cmapB);

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
ppath = [pp cfile '/'];

%% ALL
figure(1)
% TP
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,tp_ps)
colormap(cmapB)
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 4]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('TP')

% TB
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,tb_ps)
colormap(cmapB)
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 4]);
set(gcf,'renderer','painters')
title('TB')

% 
% subplot('Position',[0.5 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(lat,lon,log10(AllP))
% colormap(cmapB)
% load coastlines;                     
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([-4 4]);
% set(gcf,'renderer','painters')
% title('log10 mean All P (g m^-^2)')

% Det
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,det_ps)
colormap(cmapB)
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 4]);
set(gcf,'renderer','painters')
title('Det')
stamp('')
print('-dpng',[ppath 'Pre400_temp_det_global_ps_ln.png'])


%% ALL
figure(2)
% MZ
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,zm_ps)
colormap(cmapB)
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 4]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('MZ')

% LZ
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,zl_ps)
colormap(cmapB)
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 4]);
set(gcf,'renderer','painters')
title('LZ')

% MZloss
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,hpmz_ps)
colormap(cmapB)
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 4]);
set(gcf,'renderer','painters')
title('HPloss MZ')

% LZloss
subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,hplz_ps)
colormap(cmapB)
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-4 4]);
set(gcf,'renderer','painters')
title('HPloss LZ')
stamp('')
print('-dpng',[ppath 'Pre400_zoo_global_ps_ln.png'])


%% Histograms
edges = -6:0.25:1;

figure(4)
subplot(3,3,1)
histogram(tp,edges)
title('TP')

subplot(3,3,2)
histogram(tb,edges)
title('TB')

subplot(3,3,3)
histogram(det,edges)
title('Det')

subplot(3,3,4)
histogram(zm,edges)
title('MZ')

subplot(3,3,5)
histogram(hpmz,edges)
title('MZ hploss')

subplot(3,3,7)
histogram(zl,edges)
title('LZ')

subplot(3,3,8)
histogram(hplz,edges)
title('LZ hploss')
print('-dpng',[ppath 'Pre400_histogram_phys_bgc_ps_ln.png'])

