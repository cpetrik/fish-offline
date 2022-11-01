% Function for
% DAILY INTERPOLATION, INCLUDES UNIT CONVERSION

function ESM = daily_interp_dple_member(range,Time,Tdays,...
    TEMP_150m,TEMP_bottom,POC_FLUX_IN_bottom,LzooC_150m,Lzoo_loss_150m)

% months of interest
Tp  = (TEMP_150m(:,:,range));
Tb  = (TEMP_bottom(:,:,range));
Zm  = (LzooC_150m(:,:,range));
dZm = (Lzoo_loss_150m(:,:,range));
det = (POC_FLUX_IN_bottom(:,:,range));

% index of water cells
[ni,nj,nt] = size(LzooC_150m);
WID = find(~isnan(LzooC_150m(:,:,1)));    % spatial index of water cells excludes interior seas
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
    X = squeeze(Tp(m,n,:));
    tp = interp1(Time, X, Tdays,'linear','extrap');
    D_Tp(j,:) = tp;

    % bottom temperature (in Celcius)
    X = squeeze(Tb(m,n,:));
    tb = interp1(Time, X, Tdays,'linear','extrap');
    D_Tb(j,:) = tb;

    % meso zoo: nmolC cm-2 to g(WW) m-2
    % 1e9 nmol in 1 mol C
    % 1e4 cm2 in 1 m2
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    X = squeeze(Zm(m,n,:));
    mz = interp1(Time, X, Tdays,'linear','extrap');
    D_Zm(j,:) = mz * 1e-9 * 1e4 * 12.01 * 9.0;

    % meso zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
    % 1e9 nmol in 1 mol C
    % 1e4 cm2 in 1 m2
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    X = squeeze(dZm(m,n,:));
    zl = interp1(Time, X, Tdays,'linear','extrap');
    D_dZm(j,:) = zl * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

    % detrital flux to benthos: nmolC cm-2 s-1 to g(WW) m-2 d-1
    % 1e9 nmol in 1 mol C
    % 1e4 cm2 in 1 m2
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    X = squeeze(det(m,n,:));
    de = interp1(Time, X, Tdays,'linear','extrap');
    D_det(j,:) = de * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
end %grid cells

% Negative biomass or mortality loss from interp
D_Zm(D_Zm<0) = 0.0;
D_dZm(D_dZm<0) = 0.0;
D_det(D_det<0) = 0.0;

D_Zm(isnan(D_Zm)) = 0.0;
D_dZm(isnan(D_dZm)) = 0.0;
D_det(isnan(D_det)) = 0.0;

ESM.Tp = D_Tp;
ESM.Tb = D_Tb;
ESM.Zm = D_Zm;
ESM.dZm = D_dZm;
ESM.det = D_det;
end
