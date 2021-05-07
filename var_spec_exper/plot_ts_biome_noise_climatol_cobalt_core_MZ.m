% Plot experimental seasonal climatologies + Noise
% Sqrt-transform abund
% No transform environ (temp)

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% CORE-forced output + noise
load([fpath 'cobalt_biome_core_noise_mz_mhp.mat'],...
    'mz_clim_spl_rep_backtrans');

%% Each biome rep 220x = 20 reps, 11 slopes
b1_mz_clim = mz_clim_spl_rep_backtrans(1:220,:);
b2_mz_clim = mz_clim_spl_rep_backtrans(221:440,:);
b3_mz_clim = mz_clim_spl_rep_backtrans(441:660,:);
b4_mz_clim = mz_clim_spl_rep_backtrans(661:880,:);
b5_mz_clim = mz_clim_spl_rep_backtrans(881:1100,:);
b6_mz_clim = mz_clim_spl_rep_backtrans(1101:1320,:);
b7_mz_clim = mz_clim_spl_rep_backtrans(1321:1540,:);
b8_mz_clim = mz_clim_spl_rep_backtrans(1541:1760,:);

%% time is 150 yrs daily
[nid,nt] =size(mz_clim_spl_rep_backtrans);
time = (1:nt)/365;

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
biomes = {'SLC','SECCS','SECSS','SCoast','NLC','NECCS','NECSS','NCoast'};
alp = 0:-0.15:-1.5;
ts = 1:20:220;

%%
figure(1)
for k=1:11
    subplot(4,3,k)
    plot(time,b1_mz_clim(ts(k),:));
    xlim([148 150])
    if (k==2)
        title({biomes{1}, ['\alpha = ' num2str(alp(k))]})
    else
        title(['\alpha = ' num2str(alp(k))])
    end
end
print('-dpng',['ts_MZ_noise_biome_' biomes{1} '.png'])

%%
figure(2)
for k=1:11
    subplot(4,3,k)
    plot(time,b7_mz_clim(ts(k),:));
    xlim([148 150])
    ylim([0 1])
    if (k==2)
        title({biomes{7}, ['\alpha = ' num2str(alp(k))]})
    else
        title(['\alpha = ' num2str(alp(k))])
    end
end
print('-dpng',['ts_MZ_noise_biome_' biomes{7} '.png'])

%%
figure(3)
for k=1:11
    subplot(4,3,k)
    plot(time,b5_mz_clim(ts(k),:));
    xlim([148 150])
    ylim([0.18 0.4])
    if (k==2)
        title({biomes{5}, ['\alpha = ' num2str(alp(k))]})
    else
        title(['\alpha = ' num2str(alp(k))])
    end
end
print('-dpng',['ts_MZ_noise_biome_' biomes{5} '.png'])

