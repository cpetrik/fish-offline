% Use seasonal climatologies of phys & BGC to create
% experimental time-series

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

%%
load([fpath 'Data_core_cobalt_biome_climatol_daily_150yr.mat'], 'COBALT');

%% Use these as control forcing
COBALT.Tp = COBALT.Tp(:,1:365);
COBALT.Tb = COBALT.Tb(:,1:365);
COBALT.det = COBALT.det(:,1:365);
COBALT.Zm = COBALT.Zm(:,1:365);
COBALT.Zl = COBALT.Zl(:,1:365);
COBALT.dZm = COBALT.dZm(:,1:365);
COBALT.dZl = COBALT.dZl(:,1:365);

%%
COBALT.Tp = repelem(COBALT.Tp,220,1);
COBALT.Tb = repelem(COBALT.Tb,220,1);
COBALT.det = repelem(COBALT.det,220,1);
COBALT.Zm = repelem(COBALT.Zm,220,1);
COBALT.Zl = repelem(COBALT.Zl,220,1);
COBALT.dZm = repelem(COBALT.dZm,220,1);
COBALT.dZl = repelem(COBALT.dZl,220,1);

%%
save([fpath 'Data_core_cobalt_biome_climatol_daily_150yr.mat'], 'COBALT');

