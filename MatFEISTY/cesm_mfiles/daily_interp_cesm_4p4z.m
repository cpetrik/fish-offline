% Make mat files of interpolated time series from GFDL
% Preindust runs 1950-2100
% New vertical integrations

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/4P4Z/';
gpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';

%% Inputs
load([fpath 'g.e22a06.G1850ECOIAF_JRA_PHYS_DEV.TL319_g17.4p4z.004.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'TLAT','TLONG','TAREA','time','yr','zoo3C_150m','zoo3_loss_150m',...
    'zoo4C_150m','zoo4_loss_150m');
load([fpath 'Data_grid_POP_gx1v6_4p4z.mat'],'GRD');

%% doubles
TEMP_150m = double(TEMP_150m);
TEMP_bottom = double(TEMP_bottom);
POC_FLUX_IN_bottom = double(POC_FLUX_IN_bottom);
zoo3C_150m = double(zoo3C_150m);
zoo4C_150m = double(zoo4C_150m);
zoo3_loss_150m = double(zoo3_loss_150m);
zoo4_loss_150m = double(zoo4_loss_150m);

zoo3C_150m(zoo3C_150m<0) = 0.0;
zoo4C_150m(zoo4C_150m<0) = 0.0;
zoo3_loss_150m(zoo3_loss_150m<0) = 0.0;
zoo4_loss_150m(zoo4_loss_150m<0) = 0.0;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom<0) = 0.0;

%% Units
%poc flux: mmol/m^3 cm/s
%zoo loss: mmol/m^3/s cm
%zoo: mmolC/m^3 cm
%tp: degC
%tb: degC

% From nmolC cm-2 s-1 to g(WW) m-2 d-1
% 1e9 nmol in 1 mol C
% 1e4 cm2 in 1 m2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
% * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

% From nmolC cm-2 s-1 to g(WW) m-2
% * 1e-9 * 1e4 * 12.01 * 9.0;

%%
mos = length(time);
mstart = 1:12:mos;
mend = 12:12:mos;
nyrs = mos/12;
yrs = 1:nyrs;

Tdays=1:365;
Time=Tdays(15:30:end);

%% test that all same orientation
test1 = squeeze(double(TEMP_150m(:,:,200)));
test2 = squeeze(double(TEMP_bottom(:,:,200)));
test3 = squeeze(double(zoo3C_150m(:,:,200)));
test4 = squeeze(double(zoo4C_150m(:,:,200)));
test5 = squeeze(double(zoo3_loss_150m(:,:,200)));
test6 = squeeze(double(zoo4_loss_150m(:,:,200)));
test7 = squeeze(double(POC_FLUX_IN_bottom(:,:,200)));

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/4P4Z/';

figure
subplot(3,3,1)
pcolor(test1)
shading flat
subplot(3,3,2)
pcolor(test2)
shading flat
subplot(3,3,3)
pcolor(test7)
shading flat
subplot(3,3,4)
pcolor(test3)
shading flat
subplot(3,3,5)
pcolor(test4)
shading flat
subplot(3,3,7)
pcolor(test5)
shading flat
subplot(3,3,8)
pcolor(test6)
shading flat
print('-dpng',[pp 'cesm_4p4z.png'])

%%
% index of water cells
[ni,nj,nt] = size(TEMP_bottom);
%WID = find(~isnan(mask));  % spatial index of water cells
WID = find(~isnan(test2));  % use bottom temp instead
NID = length(WID);         % number of water cells

for y = 1:nyrs
    yr = yrs(y)
    
    Tp  = (TEMP_150m(:,:,mstart(y):mend(y)));
    Tb  = (TEMP_bottom(:,:,mstart(y):mend(y)));
    Zm  = (zoo3C_150m(:,:,mstart(y):mend(y)));
    Zl  = (zoo4C_150m(:,:,mstart(y):mend(y)));
    dZm = (zoo3_loss_150m(:,:,mstart(y):mend(y)));
    dZl = (zoo4_loss_150m(:,:,mstart(y):mend(y)));
    det = (POC_FLUX_IN_bottom(:,:,mstart(y):mend(y)));
    
    % setup FEISTY data files
    D_Tp  = nan*zeros(NID,365);
    D_Tb  = nan*zeros(NID,365);
    D_Zm  = nan*zeros(NID,365);
    D_Zl  = nan*zeros(NID,365);
    D_dZm = nan*zeros(NID,365);
    D_dZl = nan*zeros(NID,365);
    D_det = nan*zeros(NID,365);
%     D_Zm  = zeros(NID,365);
%     D_Zl  = zeros(NID,365);
%     D_dZm = zeros(NID,365);
%     D_dZl = zeros(NID,365);
%     D_det = zeros(NID,365);
    
    %% interpolate to daily resolution
    for j = 1:NID
        % indexes
        [m,n] = ind2sub([ni,nj],WID(j)); % spatial index of water cell
        
        % pelagic temperature (in Celcius)
        Y = squeeze(Tp(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Tp(j,:) = yi;
        
        % bottom temperature (in Celcius)
        Y = squeeze(Tb(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Tb(j,:) = yi;
        
        % meso zoo: nmolC cm-2 to g(WW) m-2
        % 1e9 nmol in 1 mol C
        % 1e4 cm2 in 1 m2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(Zm(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Zm(j,:) = yi * 1e-9 * 1e4 * 12.01 * 9.0;
        
        % macro zoo: nmolC cm-2 to g(WW) m-2
        Y = squeeze(Zl(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_Zl(j,:) = yi * 1e-9 * 1e4 * 12.01 * 9.0;
        
        % medium zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
        % 1e9 nmol in 1 mol C
        % 1e4 cm2 in 1 m2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(dZm(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_dZm(j,:) = yi * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
        
        % macro zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
        Y = squeeze(dZl(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_dZl(j,:) = yi * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
        
        % detrital flux to benthos: nmolC cm-2 s-1 to g(WW) m-2 d-1
        % 1e9 nmol in 1 mol C
        % 1e4 cm2 in 1 m2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(det(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_det(j,:) = yi * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
        
    end
    
    % Negative biomass or mortality loss from interp
    D_Zm(D_Zm<0) = 0.0;
    D_Zl(D_Zl<0) = 0.0;
    D_dZm(D_dZm<0) = 0.0;
    D_dZl(D_dZl<0) = 0.0;
    D_det(D_det<0) = 0.0;
    
%     D_Zm(isnan(D_Zm)) = 0.0;
%     D_Zl(isnan(D_Zl)) = 0.0;
%     D_dZm(isnan(D_dZm)) = 0.0;
%     D_dZl(isnan(D_dZl)) = 0.0;
%     D_det(isnan(D_det)) = 0.0;
    
    ESM.Tp = D_Tp;
    ESM.Tb = D_Tb;
    ESM.Zm = D_Zm;
    ESM.Zl = D_Zl;
    ESM.dZm = D_dZm;
    ESM.dZl = D_dZl;
    ESM.det = D_det;
    
    % save
    save([fpath 'Data_cesm_4p4z_daily_',num2str(yr),'.mat'], 'ESM');
    
    
end


