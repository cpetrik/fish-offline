% Calc powerspectra from COBALT ts for diff var
% By 4 biomes
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
load([spath 'cobalt_biome_core_anom_1950_2007.mat'])
[nr,nc] = size(geolon_t);
nt = length(yid);

%% reshape to save only land-free cells
% TP
xTP = tp_anom';
% TB
xTB = tb_anom';
% DET
xDet = det_anom';
% ZOOMED_INT100
xMZ = mz_anom';
% ZOOLRG_INT100
xLZ = lz_anom';
% HPLOSS MED
xHPM = hpmz_anom';
% HPLOSS LRG
xHPL = hplz_anom';

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
load([fpath 'feisty_biome_core_anom_1950_2007.mat'])
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
b = NaN*ones(size(xB,2),1);
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


%% bar graphs
figure(1)
subplot(3,3,1)
bar(tp,'k')
title('TP')
ylim([-3 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

subplot(3,3,2)
bar(tb,'k')
title('TB')
ylim([-3 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

subplot(3,3,3)
bar(det,'k')
title('Det')
ylim([-3 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

subplot(3,3,4)
bar(zm,'k')
title('MZ')
ylim([-3 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

subplot(3,3,5)
bar(zl,'k')
title('LZ')
ylim([-3 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

subplot(3,3,7)
bar(hpmz,'k')
title('hpMZ')
ylim([-3 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

subplot(3,3,8)
bar(hplz,'k')
title('hpLZ')
ylim([-3 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

%% 
figure(2)
subplot(3,3,1)
bar(sf,'k')
title('SF')
ylim([-5 0])
set(gca,'XTickLabel','')

subplot(3,3,2)
bar(sp,'k')
title('SP')
ylim([-5 0])
set(gca,'XTickLabel','')

subplot(3,3,3)
bar(sd,'k')
title('SD')
ylim([-5 0])
set(gca,'XTickLabel','')

subplot(3,3,4)
bar(mf,'k')
title('MF')
ylim([-5 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

subplot(3,3,5)
bar(mp,'k')
title('MP')
ylim([-5 0])
set(gca,'XTickLabel','')

subplot(3,3,6)
bar(md,'k')
title('MD')
ylim([-5 0])
set(gca,'XTickLabel','')

subplot(3,3,8)
bar(lp,'k')
title('LP')
ylim([-5 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

subplot(3,3,9)
bar(ld,'k')
title('LD')
ylim([-5 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

%% 
figure(3)
subplot(3,3,1)
bar(F,'k')
title('F')
ylim([-5 0])
set(gca,'XTickLabel','')

subplot(3,3,2)
bar(P,'k')
title('P')
ylim([-5 0])
set(gca,'XTickLabel','')

subplot(3,3,3)
bar(D,'k')
title('D')
ylim([-5 0])
set(gca,'XTickLabel','')

subplot(3,3,4)
bar(S,'k')
title('S')
ylim([-5 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

subplot(3,3,5)
bar(M,'k')
title('M')
ylim([-5 0])
set(gca,'XTickLabel','')

subplot(3,3,6)
bar(L,'k')
title('L')
ylim([-5 0])
set(gca,'XTickLabel','')

subplot(3,3,7)
bar(All,'k')
title('All')
ylim([-5 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

subplot(3,3,8)
bar(B,'k')
title('B')
ylim([-5 0])
set(gca,'XTickLabel',{'LC','CS','SS','Co'})

%%
save([spath 'powerspec_biome_cobalt_core_1950_2007.mat'],...
    'tp','tb','det',...
    'zl','zm','hpmz','hplz');

save([fpath 'powerspec_biome_feisty_core_1950_2007.mat'],...
    'geolat_t','geolon_t','yid','GRD',...
    'sf','sp','sd','mf','mp','md',...
    'lp','ld','B','F','P','D',...
    'S','M','L','All');

