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

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

% anomaly time series
load([fpath 'cobalt_core_anom_1950_2007.mat'])

%%
[nr,nc] = size(lon);

% t = 121:480; %last 30 yrs 1978-2007
% nt = length(t);

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
xD = reshape(det_anom, nr*nc, nt);
xD = xD(nidx, :)';
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

%% calc powerspec
tp = NaN*ones(size(xTP,2),1);
tb = NaN*ones(size(xTB,2),1);
det = NaN*ones(size(xD,2),1);
zm = NaN*ones(size(xMZ,2),1);
zl = NaN*ones(size(xLZ,2),1);
hpmz = NaN*ones(size(xHPM,2),1);
hplz = NaN*ones(size(xHPL,2),1);

%%
for i = 1:size(xTP, 2) 
    i
 
    xi = xTP(:, i);
    inan = ~isnan(xi);
    %R = log(xi(inan)); %Already log-transformed
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    tp(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xTB(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    tb(i) = b1; 
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xD(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    det(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xMZ(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    zm(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xLZ(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    zl(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xHPM(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    hpmz(i) = b1; 
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xHPL(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    hplz(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
end


%% put it back to the original map
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
save([fpath 'powerspec_cobalt_core_1950_2007.mat'],...
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
caxis([-3 3]);
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
caxis([-3 3]);
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
% caxis([-3 3]);
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
caxis([-3 3]);
set(gcf,'renderer','painters')
title('Det')
stamp('')
%print('-dpng',[ppath 'CORE_temp_det_global_ps.png'])


%% ALL
figure(1)
% MZ
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(lat,lon,zm_ps)
colormap(cmapB)
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-3 3]);
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
caxis([-3 3]);
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
caxis([-3 3]);
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
caxis([-3 3]);
set(gcf,'renderer','painters')
title('HPloss LZ')
stamp('')
%print('-dpng',[ppath 'CORE_zoo_global_ps.png'])




