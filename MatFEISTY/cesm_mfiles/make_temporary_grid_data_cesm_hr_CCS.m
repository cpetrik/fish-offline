% Make GRD file for FEISTY input from CESM HR
% Cal Curr region only

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM_HR/';

%% lat, lon
load([fpath 'g.e22.G1850ECO_JRA_HR.TL319_t13.004.FIESTY-forcing_hiresJRA_CAcurr.mat'],...
    'TLAT','TLONG','time','yr','KMT')
load([fpath 'g.e22.G1850ECO_JRA_HR.TL319_t13.004.FIESTY-forcing_hiresJRA_CAcurr_meszoo_totloss_allphytoC.mat'],...
    'LzooC_100m');

%% check orientation
figure
pcolor(KMT(:,:,1))

figure
pcolor(TLAT(:,:,1))

figure
pcolor(TLONG(:,:,1))

figure
pcolor(LzooC_100m(:,:,1))

%% create land mask
mask = squeeze(LzooC_100m(:,:,1));
mask(~isnan(mask)) = 1;

figure
pcolor(mask)

%%
LID = find(~isnan(mask(:))); 
NID = length(LID); %62598

%% Retain only water cells
ID = LID;
GRD.ID = ID;
GRD.N = NID;
GRD.LON = TLONG(ID);
GRD.LAT = TLAT(ID);
% GRD.Z   = HT(ID) * 1e-2;     %from cm to m
% GRD.area = TAREA(ID) * 1e-4; %from cm2 to m2
GRD.lmask = mask(ID);

%% Save needed variables
% save([Cdir 'gridspec_POP_gx1v6.mat'],'HT','TLAT','TLONG','TAREA','mask',...
%     'HTunits','HTlong_name','AREAunits','AREAlong_name');
save([fpath 'Data_grid_g.e22.G1850ECO_JRA_HR.TL319_t13.004.CAcurr.mat'],'GRD');
