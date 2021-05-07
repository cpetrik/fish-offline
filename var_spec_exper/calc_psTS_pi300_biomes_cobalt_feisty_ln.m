% Calc powerspectra from COBALT & FEISTY ts for diff var
% Biome means
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
spath='/Volumes/MIP/GCM_DATA/ESM2M_PI/';

% anomaly time series
load([spath 'cobalt_pi400_biomes_temp_anom.mat']);
load([spath 'cobalt_pi400_biomes_det_anom_ln.mat']);
load([spath 'cobalt_pi400_biomes_zoo_anom_ln.mat']);
load([spath 'cobalt_pi400_biomes_hploss_anom_ln.mat']);

%% reshape to save only land-free cells
% and transpose
% Exclude 1st 100 yrs (spinup for FEISTY)
xTP = tp_anom(:,1201:end)';
xTB = tb_anom(:,1201:end)';
xDet = det_anom(:,1201:end)';
xMZ = mz_anom(:,1201:end)';
xLZ = lz_anom(:,1201:end)';
xZ = z_anom(:,1201:end)';
xHP = hp_anom(:,1201:end)';
xHPM = hpmz_anom(:,1201:end)';
xHPL = hplz_anom(:,1201:end)';

clear tp_anom tb_anom mz_anom lz_anom hpmz_anom hplz_anom det_anom
clear z_anom hp_anom

%% FEISTY ---------------------------------------------------------
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

% anomaly time series
load([fpath 'feisty_pi400_biomes_anom300_ln.mat'])
nt = length(yid);
[nr,nc]=size(geolon_t);

%% already removed 1st 100 yrs (spinup)
%transpose
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

clear sf_anom sp_anom sd_anom mf_anom mp_anom md_anom lp_anom ld_anom

%%  
tp = NaN*ones(size(xTP,2),1);
tb = NaN*ones(size(xTB,2),1);
det = NaN*ones(size(xDet,2),1);
zm = NaN*ones(size(xMZ,2),1);
zl = NaN*ones(size(xLZ,2),1);
hpmz = NaN*ones(size(xHPM,2),1);
hplz = NaN*ones(size(xHPL,2),1);
z_ps = NaN*ones(4,1);
hp_ps = NaN*ones(4,1);

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
for i = 1:size(xTP, 2) 
    i
    % COBALT
    xi = xTP(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); %0.651784 seconds.
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    tp(i) = b1;
    tp_frq(:,i) = freq1;
    tp_pow(:,i) = p1;
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
    tb_frq(:,i) = freq1;
    tb_pow(:,i) = p1;
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
    det_frq(:,i) = freq1;
    det_pow(:,i) = p1;
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
    zm_frq(:,i) = freq1;
    zm_pow(:,i) = p1;
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
    zl_frq(:,i) = freq1;
    zl_pow(:,i) = p1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xZ(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    z_ps(i) = b1;
    z_frq(:,i) = freq1;
    z_pow(:,i) = p1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    xi = xHP(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    hp_ps(i) = b1; 
    hp_frq(:,i) = freq1;
    hp_pow(:,i) = p1;
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
    hpmz_frq(:,i) = freq1;
    hpmz_pow(:,i) = p1;
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
    hplz_frq(:,i) = freq1;
    hplz_pow(:,i) = p1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
    %STAGES
    xi = xsf(:, i);
    inan = ~isnan(xi);
    %R = log(xi(inan)); %Already log-transformed
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    sf(i) = b1;
    sf_frq(:,i) = freq1;
    sf_pow(:,i) = p1;
    clear R T t b int tH dR freq1 p1 b1 int1 data
    
    xi = xsp(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    sp(i) = b1; 
    sp_frq(:,i) = freq1;
    sp_pow(:,i) = p1;
    clear R T t b int tH dR freq1 p1 b1 int1 data
    
    xi = xsd(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    sd(i) = b1;
    sd_frq(:,i) = freq1;
    sd_pow(:,i) = p1;
    clear R T t b int tH dR freq1 p1 b1 int1 data
    
    xi = xmf(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    mf(i) = b1;
    mf_frq(:,i) = freq1;
    mf_pow(:,i) = p1;
    clear R T t b int tH dR freq1 p1 b1 int1 data
    
    xi = xmp(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    mp(i) = b1;
    mp_frq(:,i) = freq1;
    mp_pow(:,i) = p1;
    clear R T t b int tH dR freq1 p1 b1 int1 data
    
    xi = xmd(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    md(i) = b1; 
    md_frq(:,i) = freq1;
    md_pow(:,i) = p1;
    clear R T t b int tH dR freq1 p1 b1 int1 data
    
    xi = xlp(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    lp(i) = b1;
    lp_frq(:,i) = freq1;
    lp_pow(:,i) = p1;
    clear R T t b int tH dR freq1 p1 b1 int1 data
    
    xi = xld(:, i);
    inan = ~isnan(xi);
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra_TS(dR,12,0);
    ld(i) = b1;
    ld_frq(:,i) = freq1;
    ld_pow(:,i) = p1;
    clear R T t b int tH dR freq1 p1 b1 int1 data
    
    %TYPES
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
    B_frq(:,i) = freq1;
    B_pow(:,i) = p1;
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
    F_frq(:,i) = freq1;
    F_pow(:,i) = p1;
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
    P_frq(:,i) = freq1;
    P_pow(:,i) = p1;
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
    D_frq(:,i) = freq1;
    D_pow(:,i) = p1;
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
    S_frq(:,i) = freq1;
    S_pow(:,i) = p1;
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
    M_frq(:,i) = freq1;
    M_pow(:,i) = p1;
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
    L_frq(:,i) = freq1;
    L_pow(:,i) = p1;
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
    All_frq(:,i) = freq1;
    All_pow(:,i) = p1;
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
save([fpath 'powerspec_pi400_cobalt_fesity_biomes_300yr_ln.mat'],...
    'tp_ps','tb_ps','det_ps','z_ps','hp_ps',...
    'zl_ps','zm_ps','hpmz_ps','hplz_ps',...
    'sf_ps','sp_ps','sd_ps','mf_ps','mp_ps','md_ps',...
    'lp_ps','ld_ps','yid','biome4_hist','biome',...
    'B_ps','F_ps','P_ps','D_ps',...
    'S_ps','M_ps','L_ps','All_ps');

save([fpath 'freq_power_pi400_cobalt_fesity_biomes_300yr_ln.mat'],...
    'tp_frq','tb_frq','det_frq','z_frq','hp_frq',...
    'zl_frq','zm_frq','hpmz_frq','hplz_frq',...
    'sf_frq','sp_frq','sd_frq','mf_frq','mp_frq','md_frq',...
    'lp_frq','ld_frq',...
    'B_frq','F_frq','P_frq','D_frq',...
    'S_frq','M_frq','L_frq','All_frq',...
    'tp_pow','tb_pow','det_pow','z_pow','hp_pow',...
    'zl_pow','zm_pow','hpmz_pow','hplz_pow',...
    'sf_pow','sp_pow','sd_pow','mf_pow','mp_pow','md_pow',...
    'lp_pow','ld_pow','yid','biome4_hist','biome',...
    'B_pow','F_pow','P_pow','D_pow',...
    'S_pow','M_pow','L_pow','All_pow');


%% Figures
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
ppath = [pp cfile '/'];

%% Bar graph?
figure(1)
subplot(3,3,1)
bar(tp_ps,'k')
title('TP')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-3.5 0])

subplot(3,3,2)
bar(tb_ps,'k')
title('TB')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-3.5 0])

subplot(3,3,3)
bar(det_ps,'k')
title('Det')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-3.5 0])

subplot(3,3,4)
bar(zm_ps,'k')
title('MZ')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-3.5 0])

subplot(3,3,7)
bar(hpmz_ps,'k')
title('MZ hploss')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-3.5 0])

subplot(3,3,5)
bar(zl_ps,'k')
title('LZ')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-3.5 0])

subplot(3,3,8)
bar(hplz_ps,'k')
title('LZ hploss')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-3.5 0])

subplot(3,3,6)
bar(z_ps,'k')
title('Z')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-3.5 0])

subplot(3,3,9)
bar(hp_ps,'k')
title('Z hploss')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-3.5 0])
print('-dpng',[ppath 'Pre300_biomes_bar_phys_bgc_ps_ln.png'])

%%
figure(2)
subplot(3,3,1)
bar(sf_ps,'k')
title('SF')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,2)
bar(sp_ps,'k')
title('SP')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,3)
bar(sd_ps,'k')
title('SD')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,4)
bar(mf_ps,'k')
title('MF')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,5)
bar(mp_ps,'k')
title('MP')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,6)
bar(md_ps,'k')
title('MD')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,8)
bar(lp_ps,'k')
title('LP')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,9)
bar(ld_ps,'k')
title('LD')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])
print('-dpng',[ppath 'Pre300_biomes_bar_fish_stages_ps_ln.png'])

%%
figure(3)
subplot(3,3,1)
bar(F_ps,'k')
title('F')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,2)
bar(P_ps,'k')
title('P')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,3)
bar(D_ps,'k')
title('D')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,4)
bar(S_ps,'k')
title('S')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,5)
bar(M_ps,'k')
title('M')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,6)
bar(L_ps,'k')
title('L')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,7)
bar(All_ps,'k')
title('All')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])

subplot(3,3,8)
bar(B_ps,'k')
title('B')
set(gca,'XTickLabel',{'LC','CS','SS','Co'},'ylim',[-5.5 0])
print('-dpng',[ppath 'Pre300_biomes_bar_fish_groups_ps_ln.png'])

%% Plot periodogram?
f=12;

for i=1:4
    figure(4)
    subplot(2,2,i)
    plot(log10(tp_frq(:,i)),log10(tp_pow(:,i)),'k')
    set( gca, 'XTick',log10(1./([ 200 100 50 20 10 5 2 1 6/12 2/12].*f)), ...
        'XTickLabel', { '200a';'100a'; '50a'; '20a'; '10a'; '5a';
        '2a'; '1a'; '6m'; '2m' }), ylim([-8 2])
    title([biome{i} ' Periodogram Using FFT'])
    xlabel('Period')
    ylabel('Power/Frequency (B/Hz)')
    grid on;
    
    figure(5)
    subplot(2,2,i)
    plot(log10(tb_frq(:,i)),log10(tb_pow(:,i)),'k')
    set( gca, 'XTick',log10(1./([ 200 100 50 20 10 5 2 1 6/12 2/12].*f)), ...
        'XTickLabel', { '200a';'100a'; '50a'; '20a'; '10a'; '5a';
        '2a'; '1a'; '6m'; '2m' }), ylim([-11 0])
    title([biome{i} ' Periodogram Using FFT'])
    xlabel('Period')
    ylabel('Power/Frequency (B/Hz)')
    grid on;
    
    figure(6)
    subplot(2,2,i)
    plot(log10(det_frq(:,i)),log10(det_pow(:,i)),'k')
    set( gca, 'XTick',log10(1./([ 200 100 50 20 10 5 2 1 6/12 2/12].*f)), ...
        'XTickLabel', { '200a';'100a'; '50a'; '20a'; '10a'; '5a';
        '2a'; '1a'; '6m'; '2m' }), ylim([-8 0])
    title([biome{i} ' Periodogram Using FFT'])
    xlabel('Period')
    ylabel('Power/Frequency (B/Hz)')
    grid on;
    
    figure(7)
    subplot(2,2,i)
    plot(log10(zm_frq(:,i)),log10(zm_pow(:,i)),'k')
    set( gca, 'XTick',log10(1./([ 200 100 50 20 10 5 2 1 6/12 2/12].*f)), ...
        'XTickLabel', { '200a';'100a'; '50a'; '20a'; '10a'; '5a';
        '2a'; '1a'; '6m'; '2m' }), ylim([-8 0])
    title([biome{i} ' Periodogram Using FFT'])
    xlabel('Period')
    ylabel('Power/Frequency (B/Hz)')
    grid on;
end
   
