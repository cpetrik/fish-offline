% FEISTY low fish biomass
% <1 fish pr grid cell

clear
close all

%% Grid data - CESM
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);
ID = GRD.ID;

%TAREA units 'cm^2' -> m2
AREA_OCN = TAREA * 1e-4;
area = AREA_OCN(ID);

%% Fish size
param.M_s = 10^((log10(0.001)+log10(0.5))/2);  %0.0224
param.M_m = 10^((log10(0.5)+log10(250))/2);    %11.1803
param.M_l = 10^((log10(250)+log10(125000))/2); %5.5902e3

%% Means

marea = mean(area);

Mone = param.M_m ./ marea;
Lone = param.M_l ./ marea;
