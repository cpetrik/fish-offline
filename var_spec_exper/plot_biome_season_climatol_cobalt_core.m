% Plot experimental seasonal climatologies 
% Sqrt-transform abund
% No transform environ (temp)

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% CORE-forced output
load([fpath 'cobalt_biome_core_climatol_1950_2007.mat']);
load([fpath 'cobalt_biome_core_anom_1950_2007.mat']);

%% un-sqrt-transform
det_clim = det_clim.^2;
mz_clim = mz_clim.^2;
lz_clim = lz_clim.^2;
mhp_clim = mhp_clim.^2;
lhp_clim = lhp_clim.^2;

%%
cm9=[0.5 0.5 0.5 ;...   %grey
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0 0.7 0]; %...    %green
    %0 0 0];...      %black
    
set(groot,'defaultAxesColorOrder',cm9);

%% 
figure(1)
subplot(3,3,1)
plot(1:12,tp_clim,'LineWidth',2);
xlim([1 12])
title('TP')

subplot(3,3,2)
plot(1:12,tb_clim,'LineWidth',2);
xlim([1 12])
title('TB')

subplot(3,3,3)
plot(1:12,log10(det_clim),'LineWidth',2);
xlim([1 12])
title('log_1_0 Det')

subplot(3,3,4)
plot(1:12,mz_clim,'LineWidth',2);
xlim([1 12])
title('MZ')

subplot(3,3,5)
plot(1:12,mhp_clim,'LineWidth',2);
xlim([1 12])
title('MZ HPloss')

subplot(3,3,7)
plot(1:12,lz_clim,'LineWidth',2);
xlim([1 12])
title('LZ')

subplot(3,3,8)
plot(1:12,lhp_clim,'LineWidth',2);
xlim([1 12])
title('LZ HPloss')

subplot(3,3,6)
plot(1:12,tb_clim,'LineWidth',2);
ylim([20 30])
set(gca,'XTickLabel',[],'YTickLabel',[])
legend('SLC','SECCS','SECSS','SCoast','NLC','NECCS','NECSS','NCoast')
legend('location','west')
print('-dpng',[pp 'ts_CORE_biome_climatol.png'])
