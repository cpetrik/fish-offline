% Make mat files of interpolated time series from CESM 
% CESM 1 degree companion run for JRA55-forced high res
% Has dedicated mesozooplankton and loss
% Add last mon y-1 and 1st mo y+1 to interp
% Remove grid cells of inland seas

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/4P2Z/';
gpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

%% Inputs
load([fpath 'g.e22.GOMIPECOIAF_JRA-1p4-2018.TL319_g17.4p2z.001.marbl-scope-FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'mesozooC_150m','mesozooC_150m_units','mesozoo_loss_150m','mesozoo_loss_150m_units',...
    'TLAT','TLONG','TAREA','time','yr');
load([gpath 'gridspec_POP_gx1v6_noSeas.mat'],'mask');

YR = yr;
clear yr

%% nans & zeros
TEMP_150m = double(TEMP_150m);
TEMP_bottom = double(TEMP_bottom);
POC_FLUX_IN_bottom = double(POC_FLUX_IN_bottom);

TEMP_bottom(TEMP_bottom >= 9.9e+36) = nan;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;
mesozooC_150m(TEMP_bottom >= 9.9e+36) = nan;
mesozoo_loss_150m(TEMP_bottom >= 9.9e+36) = nan;

mesozooC_150m(mesozooC_150m<0) = 0.0;
mesozoo_loss_150m(mesozoo_loss_150m<0) = 0.0;
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
%yrs = 1:nyrs;
yrs = ceil(YR(1)):YR(end);

Tdays=1:365;

%% test that all same orientation
test1 = squeeze(double(TEMP_150m(:,:,200)));
test2 = squeeze(double(TEMP_bottom(:,:,200)));
test3 = squeeze(double(mesozooC_150m(:,:,200)));
test4 = squeeze(double(mesozoo_loss_150m(:,:,200)));
test5 = squeeze(double(POC_FLUX_IN_bottom(:,:,200)));

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/4P2Z/';

figure
subplot(2,3,1)
pcolor(test1)
shading flat
subplot(2,3,2)
pcolor(test2)
shading flat
subplot(2,3,3)
pcolor(test5)
shading flat
subplot(2,3,4)
pcolor(test3)
shading flat
subplot(2,3,5)
pcolor(test4)
shading flat
%print('-dpng',[pp 'cesm_fosi.png'])

%%
for y = 1:nyrs
    yr = yrs(y)
    
    if y==1
        range = mstart(y):(mend(y)+1);
        Time=15:30:395;
    elseif y==nyrs
        range = (mstart(y)-1):mend(y);
        Time=-15:30:365;
    else
        range = (mstart(y)-1):(mend(y)+1);
        Time=-15:30:395;
    end
    
    Tp  = (TEMP_150m(:,:,range));
    Tb  = (TEMP_bottom(:,:,range));
    Zm  = (mesozooC_150m(:,:,range));
    dZm = (mesozoo_loss_150m(:,:,range));
    det = (POC_FLUX_IN_bottom(:,:,range));
    
    % index of water cells
    [ni,nj,nt] = size(TEMP_bottom);
    WID = find(~isnan(mask));  % spatial index of water cells
    NID = length(WID);         % number of water cells
    
    % setup FEISTY data files
    D_Tp  = nan*zeros(NID,365);
    D_Tb  = nan*zeros(NID,365);
    D_Zm  = nan*zeros(NID,365);
    D_dZm = nan*zeros(NID,365);
    D_det = nan*zeros(NID,365);
    
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
        
        % medium zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
        % 1e9 nmol in 1 mol C
        % 1e4 cm2 in 1 m2
        % 12.01 g C in 1 mol C
        % 1 g dry W in 9 g wet W (Pauly & Christiansen)
        Y = squeeze(dZm(m,n,:));
        yi = interp1(Time, Y, Tdays,'linear','extrap');
        D_dZm(j,:) = yi * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
        %Y and D_dZm not in same units!
        %D_dZm(j,:) = max(D_dZm(j,:), min(Y)); 
        
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
    D_dZm(D_dZm<0) = 0.0;
    D_det(D_det<0) = 0.0;
    
%     D_Zm(isnan(D_Zm)) = 0.0;
%     D_dZm(isnan(D_dZm)) = 0.0;
%     D_det(isnan(D_det)) = 0.0;
    
    ESM.Tp = D_Tp;
    ESM.Tb = D_Tb;
    ESM.Zm = D_Zm;
    ESM.dZm = D_dZm;
    ESM.det = D_det;
    
    % save
    save([fpath 'Data_cesm_4p2z_daily_',num2str(yr),'.mat'], 'ESM');
    
    
end


