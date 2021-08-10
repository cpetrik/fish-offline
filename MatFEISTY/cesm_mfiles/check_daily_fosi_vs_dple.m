% Check that daily interp files got made correctly
% And that units are correct

clear all
close all

%% FOSI
fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([fpath 'gridspec_POP_gx1v6.mat']);
load([fpath 'Data_grid_POP_gx1v6.mat'],'GRD');
load([fpath 'Data_cesm_fosi_daily_7.mat']); %1954?

CI = ESM;
clear ESM

%% 4P4Z
dpath='/Volumes/MIP/GCM_DATA/CESM/DPLE/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/DPLE/';

load([dpath 'Data_grid_POP_gx1v6_DPLE.mat']);
load([dpath 'Data_cesm_dple_daily_M1_Y1954_1.mat']);

HI = ESM;
clear ESM

%%
zm_ci=mean(CI.Zm);
d_ci=mean(CI.det);
tp_ci=mean(CI.Tp);
tb_ci=mean(CI.Tb);
dzm_ci=mean(CI.dZm);

zm_hi=mean(HI.Zm);
d_hi=mean(HI.det);
tp_hi=mean(HI.Tp);
tb_hi=nanmean(HI.Tb);
dzm_hi=mean(HI.dZm);

%%
ntp = isnan(HI.Tp);
ntb = isnan(HI.Tb);
ndet = isnan(HI.det);
nmz = isnan(HI.Zm);
ndmz = isnan(HI.dZm);

sum(ntp(:)) 
sum(ntb(:)) %4745
sum(ndet(:))
sum(nmz(:))
sum(ndmz(:))


%% 
%Det
figure(1)
subplot(2,2,1)
plot(d_ci,'k'); hold on;
plot(d_hi,'b'); hold on;
title('yr 1 Det')
legend({'fosi','dple'})

%Temp P
subplot(2,2,3)
plot(tp_ci,'k'); hold on;
plot(tp_hi,'b'); hold on;
title('yr 1 Tp')
%legend({'fosi','dple'})

%Temp B
subplot(2,2,4)
plot(tb_ci,'k'); hold on;
plot(tb_hi,'b'); hold on;
title('yr 1 Tb')
%legend({'fosi','dple'})
stamp('')
print('-dpng',[pp 'CESM_FOSI_DPLE_M1_Y1954_daily_yr1.png'])

%% 
%MZ
figure(2)
subplot(2,1,1)
plot(zm_ci,'k'); hold on;
plot(zm_hi,'b'); hold on;
title('yr 1 MZ')
legend({'fosi LZ','dple LZ'})

subplot(2,1,2)
plot(dzm_ci,'k'); hold on;
plot(dzm_hi,'b'); hold on;
title('yr 1 LZ loss')

stamp('')
print('-dpng',[pp 'CESM_FOSI_DPLE_M1_Y1954_daily_zoo_yr1.png'])

%%
figure(3)
subplot(2,2,1)
plot(log10(zm_ci),'k'); hold on;
plot(log10(zm_hi),'b'); hold on;
title('yr 1 log10 LZ')
legend({'fosi','dple'})

subplot(2,2,3)
plot(log10(dzm_ci),'k'); hold on;
plot(log10(dzm_hi),'b'); hold on;
title('yr 1 log10 LZ loss')
legend({'fosi','dple'})

%% DPLE zoo too low
mean(CI.Zm(:))/mean(HI.Zm(:)) %8.9676



