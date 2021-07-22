% Use phyto size classes to split zoop into meso and micro

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';

%% Plankton inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'zooC_150m','zoo_loss_150m','diatC_150m','spC_150m');

%%
% Size:       320x384x816
% Dimensions: nlon,nlat,time
FillValue   = NaN;
missing_value = 9.969209968386869e+36;

% zooC_150m
zooC_150m_units = 'mmol/m^3 cm';
zooC_150m_name  = 'zoo biomass';

% zoo_loss_150m
zoo_loss_150m_units = 'mmol/m^3/s cm';
zoo_loss_150m_name  = 'zoo loss';

% diatC_150m
diatC_150m_units = 'mmol/m^3 cm';
diatC_150m_name  = 'diatom biomass';

% spC_150m
spC_150m_units = 'mmol/m^3 cm';
spC_150m_name  = 'SP biomass';

%% doubles
zooC_150m = double(zooC_150m); 
zoo_loss_150m = double(zoo_loss_150m); 
diatC_150m = double(diatC_150m); 
spC_150m = double(spC_150m); 

%% fraction large phyto -> frac large zoo
fracL = diatC_150m ./ (diatC_150m + spC_150m);
LzooC_150m = fracL .* zooC_150m;
Lzoo_loss_150m = fracL .* zoo_loss_150m;

%%
save([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo.mat'],...
    'fracL','LzooC_150m','Lzoo_loss_150m','zooC_150m_units','zoo_loss_150m_units');





