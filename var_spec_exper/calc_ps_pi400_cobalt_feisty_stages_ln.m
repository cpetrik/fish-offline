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

% load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
%     'geolon_t','geolat_t');
% load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
%     'GRD');

% anomaly time series
load([fpath 'feisty_pi400_anom_ln.mat'],'geolat_t','geolon_t','yid','grid',...
    'sf_anom','sp_anom','sd_anom','mf_anom','mp_anom','md_anom',...
    'lp_anom','ld_anom')
nt = length(yid);
[nr,nc]=size(geolon_t);

%% remove 1st 50 yrs (spinup)
% may need to remove transpose
xsf = sf_anom(:,601:end)';
xsp = sp_anom(:,601:end)';
xsd = sd_anom(:,601:end)';
xmf = mf_anom(:,601:end)';
xmp = mp_anom(:,601:end)'; 
xmd = md_anom(:,601:end)';
xlp = lp_anom(:,601:end)'; 
xld = ld_anom(:,601:end)'; 

%% 
sf = NaN*ones(size(xsf,2),1);
sp = NaN*ones(size(xsp,2),1);
sd = NaN*ones(size(xsd,2),1);
mf = NaN*ones(size(xmf,2),1);
mp = NaN*ones(size(xmp,2),1);
md = NaN*ones(size(xmd,2),1);
lp = NaN*ones(size(xlp,2),1);
ld = NaN*ones(size(xld,2),1);

clear sf_anom sp_anom sd_anom mf_anom mp_anom md_anom lp_anom ld_anom

%% calc powerspec
for i = 1:size(xsp, 2) 
    i
   
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
    clear R T t b int tH dR freq1 p1 b1 int1 data
    
end

%% put FEISTY back to the original map
%nidx = GRD.ID;
nidx = grid(:,1);

sf_ps = NaN*ones(nr,nc);
sp_ps = NaN*ones(nr,nc);
sd_ps = NaN*ones(nr,nc);
mf_ps = NaN*ones(nr,nc);
mp_ps = NaN*ones(nr,nc);
md_ps = NaN*ones(nr,nc);
lp_ps = NaN*ones(nr,nc);
ld_ps = NaN*ones(nr,nc);

sf_ps(nidx) = sf;
sp_ps(nidx) = sp;
sd_ps(nidx) = sd;
mf_ps(nidx) = mf;
mp_ps(nidx) = mp;
md_ps(nidx) = md;
lp_ps(nidx) = lp;
ld_ps(nidx) = ld;

%%
save([fpath 'powerspec_feisty_stages_pi400_ln.mat'],...
    'geolat_t','geolon_t','yid','grid',...
    'sf_ps','sp_ps','sd_ps','mf_ps','mp_ps','md_ps',...
    'lp_ps','ld_ps');

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
ppath = [pp cfile '/'];


%% Histograms
edges = -6:0.25:1;

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
print('-dpng',[ppath 'Pre400_histogram_fish_stages_ps_ln.png'])



