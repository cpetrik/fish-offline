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

%% COBALT --------------------------------------------------------
spath='/Volumes/MIP/GCM_DATA/CORE-forced/';

% anomaly time series
load([spath 'cobalt_core_anom_1950_2007.mat'])
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

%% 
tp = NaN*ones(size(xTP,2),1);
tb = NaN*ones(size(xTB,2),1);
det = NaN*ones(size(xDet,2),1);
zm = NaN*ones(size(xMZ,2),1);
zl = NaN*ones(size(xLZ,2),1);
hpmz = NaN*ones(size(xHPM,2),1);
hplz = NaN*ones(size(xHPL,2),1);


%% FEISTY ---------------------------------------------------------
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

[nr,nc]=size(geolon_t);

% anomaly time series
load([fpath 'feisty_core_anom_1950_2007.mat'])
nt = length(yid);

%% use only land-free cells
% may need to remove nans and transpose
xsf = sf_anom';
xsp = sp_anom';
xsd = sd_anom';
xmf = mf_anom';
xmp = mp_anom'; 
xmd = md_anom';
xlp = lp_anom'; 
xld = ld_anom'; 
xB  = B_anom';
xF  = F_anom';
xP  = P_anom';
xD  = D_anom';
xS  = S_anom';
xM  = M_anom';
xL  = L_anom';
xall = all_anom';

%% 
sf = NaN*ones(size(xsf,2),1);
sp = NaN*ones(size(xsp,2),1);
sd = NaN*ones(size(xsd,2),1);
mf = NaN*ones(size(xmf,2),1);
mp = NaN*ones(size(xmp,2),1);
md = NaN*ones(size(xmd,2),1);
lp = NaN*ones(size(xlp,2),1);
ld = NaN*ones(size(xld,2),1);
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
    
    xi = xDet(:, i);
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
    
    
    xi = xsf(:, i);
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
    sf(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xsp(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    sp(i) = b1; 
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xsd(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    sd(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xmf(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    mf(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xmp(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    mp(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xmd(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    md(i) = b1; 
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xlp(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    lp(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xld(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    ld(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xB(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    B(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xF(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    F(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xP(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    P(i) = b1;
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
    D(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xS(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    S(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xM(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    M(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xL(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    L(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xall(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    All(i) = b1;
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

%% put FEISTY back to the original map
nidx = GRD.ID;

sf_ps = NaN*ones(nr,nc);
sp_ps = NaN*ones(nr,nc);
sd_ps = NaN*ones(nr,nc);
mf_ps = NaN*ones(nr,nc);
mp_ps = NaN*ones(nr,nc);
md_ps = NaN*ones(nr,nc);
lp_ps = NaN*ones(nr,nc);
ld_ps = NaN*ones(nr,nc);
B_ps = NaN*ones(nr,nc);
F_ps = NaN*ones(nr,nc);
P_ps = NaN*ones(nr,nc);
D_ps = NaN*ones(nr,nc);
S_ps = NaN*ones(nr,nc);
M_ps = NaN*ones(nr,nc);
L_ps = NaN*ones(nr,nc);
All_ps = NaN*ones(nr,nc);

sf_ps(nidx) = sf;
sp_ps(nidx) = sp;
sd_ps(nidx) = sd;
mf_ps(nidx) = mf;
mp_ps(nidx) = mp;
md_ps(nidx) = md;
lp_ps(nidx) = lp;
ld_ps(nidx) = ld;
B_ps(nidx)  = B;
F_ps(nidx)  = F;
P_ps(nidx)  = P;
D_ps(nidx)  = D;
S_ps(nidx)  = S;
M_ps(nidx)  = M;
L_ps(nidx)  = L;
All_ps(nidx) = All;

%%
save([spath 'powerspec_cobalt_core_1950_2007_sqrt.mat'],...
    'tp_ps','tb_ps','det_ps',...
    'zl_ps','zm_ps','hpmz_ps','hplz_ps');

save([fpath 'powerspec_feisty_core_1950_2007_sqrt.mat'],...
    'geolat_t','geolon_t','yid','GRD',...
    'sf_ps','sp_ps','sd_ps','mf_ps','mp_ps','md_ps',...
    'lp_ps','ld_ps','B_ps','F_ps','P_ps','D_ps',...
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

%%
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
print('-dpng',[ppath 'CORE_temp_det_global_ps.png'])


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
print('-dpng',[ppath 'CORE_zoo_global_ps.png'])

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
caxis([-5 5]);
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
caxis([-5 5]);
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
caxis([-5 5]);
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
caxis([-5 5]);
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
caxis([-5 5]);
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
caxis([-5 5]);
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
caxis([-5 5]);
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
caxis([-5 5]);
set(gcf,'renderer','painters')
text(0,1.75,'All fish','HorizontalAlignment','center')

print('-dpng',[ppath 'CORE_fish_groups_global_ps.png'])

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
print('-dpng',[ppath 'CORE_histogram_phys_bgc_ps_sqrt.png'])

%%
figure(5)
subplot(3,3,1)
histogram(sf,edges)
title('SF')

subplot(3,3,2)
histogram(sp,edges)
title('SP')

subplot(3,3,3)
histogram(sd,edges)
title('SD')

subplot(3,3,4)
histogram(mf,edges)
title('MF')

subplot(3,3,5)
histogram(mp,edges)
title('MP')

subplot(3,3,6)
histogram(md,edges)
title('MD')

subplot(3,3,8)
histogram(lp,edges)
title('LP')

subplot(3,3,9)
histogram(ld,edges)
title('LD')
print('-dpng',[ppath 'CORE_histogram_fish_stages_ps_sqrt.png'])

%% 
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
print('-dpng',[ppath 'CORE_histogram_fish_groups_ps_sqrt.png'])


