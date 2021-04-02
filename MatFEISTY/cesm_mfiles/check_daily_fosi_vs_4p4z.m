% Check that daily interp files got made correctly
% And that units are correct

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/4P4Z/';

%% FOSI
fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([fpath 'gridspec_POP_gx1v6.mat']);
load([fpath 'Data_grid_POP_gx1v6.mat'],'GRD');
load([fpath 'Data_cesm_fosi_daily_60.mat']);

CI = ESM;
clear ESM

%% 4P4Z
zpath='/Volumes/MIP/GCM_DATA/CESM/4P4Z/';
load([zpath 'Data_cesm_4p4z_daily_60.mat']);

HI = ESM;
clear ESM

%%
zm_ci=mean(CI.Zm);
d_ci=mean(CI.det);
tp_ci=mean(CI.Tp);
tb_ci=mean(CI.Tb);
dzm_ci=mean(CI.dZm);

zm_hi=mean(HI.Zm);
zl_hi=mean(HI.Zl);
d_hi=mean(HI.det);
tp_hi=mean(HI.Tp);
tb_hi=mean(HI.Tb);
dzm_hi=mean(HI.dZm);
dzl_hi=mean(HI.dZl);

tz_hi = zm_hi + zl_hi;
dtz_hi = dzm_hi + dzl_hi;

%% 
%Det
figure(1)
subplot(2,2,1)
plot(d_ci,'k'); hold on;
plot(d_hi,'b'); hold on;
title('yr 60 Det')
legend({'fosi','4p4z'})

%Temp P
subplot(2,2,3)
plot(tp_ci,'k'); hold on;
plot(tp_hi,'b'); hold on;
title('yr 60 Tp')
%legend({'fosi','4p4z'})

%Temp B
subplot(2,2,4)
plot(tb_ci,'k'); hold on;
plot(tb_hi,'b'); hold on;
title('yr 60 Tb')
%legend({'fosi','4p4z'})
stamp('')
%print('-dpng',[pp 'CESM_daily_yr60.png'])

%% 
%MZ
figure(2)
subplot(2,3,1)
plot(zm_ci,'k'); hold on;
plot(zm_hi,'b'); hold on;
title('yr 60 MZ')
legend({'fosi LZ','4p4z Z3'})

subplot(2,3,4)
plot(dzm_ci,'k'); hold on;
plot(dzm_hi,'b'); hold on;
title('yr 60 hplMZ')

%LZ
subplot(2,3,2)
plot(zm_ci,'k'); hold on;
plot(zl_hi,'b'); hold on;
title('yr 60 LZ')
legend({'fosi LZ','4p4z Z4'})

subplot(2,3,5)
plot(dzm_ci,'k'); hold on;
plot(dzl_hi,'b'); hold on;
title('yr 60 hplLZ')
%legend({'fosi','4p4z'})

%AllZ
subplot(2,3,3)
plot(zm_ci,'k'); hold on;
plot(tz_hi,'b'); hold on;
title('yr 60 LZ')
legend({'fosi LZ','4p4z Z3+Z4'})

%Temp B
subplot(2,3,6)
plot(dzm_ci,'k'); hold on;
plot(dtz_hi,'b'); hold on;
title('yr 60 hplLZ')
%legend({'fosi','4p4z'})

stamp('')
%print('-dpng',[pp 'CESM_daily_zoo_yr60.png'])

%%
figure(3)
subplot(2,2,1)
plot(log10(zm_ci),'k'); hold on;
plot(log10(zm_hi),'b'); hold on;
title('yr 60 log10 MZ')
legend({'fosi','4p4z Z3'})

subplot(2,2,3)
plot(log10(zm_ci),'k'); hold on;
plot(log10(zl_hi),'b'); hold on;
legend({'fosi','4p4z Z4'})
title('yr 60 log10 Z')


%% 4P4Z zoo very high
mean(CI.Zm(:))/mean(HI.Zl(:)) %0.3419
mean(CI.Zm(:))/mean(HI.Zm(:)) %0.3419
mean(CI.Zm(:))/ (mean(HI.Zm(:)) + mean(HI.Zm(:))) %0.1709


