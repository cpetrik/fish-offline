% Calc powerspectra from COBALT ts for diff var
% Steps 3 & 5
% How get signif of trend? (Step 4)

% Steps
% 1. Calc seasonal climatology
% 2. Create anomaly time series by removing climatol - TEST TRANSFORM
% 3. Calc Theil-Sen linear trend
% 4. Remove trend if significant
% 5. Calc spectral slopes

clear all
close all

trans = 'raw';

%% COBALT --------------------------------------------------------
spath='/Volumes/MIP/GCM_DATA/CORE-forced/';

% anomaly time series
load([spath 'cobalt_core_anom_1950_2007_',trans,'.mat'])
[nr,nc] = size(lon);
%[nr,nc] = size(geolon_t);
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
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

[nr,nc]=size(geolon_t);

% anomaly time series
load([fpath 'feisty_core_anom_1950_2007_',trans,'.mat'])
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
for i = 1:9e3:48111%1:size(xP, 2) 
    i
 
    xi = xTP(:, i);
    inan = ~isnan(xi);
    %R = log(xi(inan)); %Not temp
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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
    R = log(xi(inan)); 
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


%% 
tp_ps = tp;
tb_ps = tb;
det_ps = det;
zm_ps = zm;
zl_ps = zl;
hpmz_ps = hpmz;
hplz_ps = hplz;

%
sf_ps = sf;
sp_ps = sp;
sd_ps = sd;
mf_ps = mf;
mp_ps = mp;
md_ps = md;
lp_ps = lp;
ld_ps = ld;
B_ps  = B;
F_ps  = F;
P_ps  = P;
D_ps  = D;
S_ps  = S;
M_ps  = M;
L_ps  = L;
All_ps = All;

%%
save([spath 'powerspec_cobalt_core_1950_2007_',trans,'_test.mat'],...
    'tp_ps','tb_ps','det_ps',...
    'zl_ps','zm_ps','hpmz_ps','hplz_ps');

save([fpath 'powerspec_feisty_core_1950_2007_',trans,'_test.mat'],...
    'geolat_t','geolon_t','yid','GRD',...
    'sf_ps','sp_ps','sd_ps','mf_ps','mp_ps','md_ps',...
    'lp_ps','ld_ps','B_ps','F_ps','P_ps','D_ps',...
    'S_ps','M_ps','L_ps','All_ps');

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
ppath = [pp cfile '/'];

% Histograms
edges = -6:0.5:1;

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
print('-dpng',[ppath 'CORE_histogram_phys_bgc_ps_',trans,'_test.png'])

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
print('-dpng',[ppath 'CORE_histogram_fish_stages_ps_',trans,'_test.png'])

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
print('-dpng',[ppath 'CORE_histogram_fish_groups_ps_',trans,'_test.png'])


