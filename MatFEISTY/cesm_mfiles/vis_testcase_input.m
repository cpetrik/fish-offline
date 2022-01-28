% Visualize testcase forcing

clear all
close all

%%
fpath='/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/';
load([fpath 'cyclic_test_forcing.mat']);

%% bathym
figure(1)
plot(X,-1*bathymetry,'.-b')
ylabel('depth')
print('-dpng',[fpath 'testcase_bathym.png'])

%% temp
figure(2)
pcolor(T_pelagic')
shading flat
ylabel('time')
xlabel('loc')
colorbar
caxis([10 18])
print('-dpng',[fpath 'testcase_TP.png'])

figure(3)
pcolor(T_bottom')
shading flat
ylabel('time')
xlabel('loc')
colorbar
caxis([3.8 4.2])
print('-dpng',[fpath 'testcase_TB.png'])

%% Prey at loc 1
figure(4)
plot(time,poc_flux_bottom(1,:),'-b')
ylabel('POC')
ylim([0.018 0.028])
print('-dpng',[fpath 'testcase_POC_loc1.png'])

figure(5)
plot(time,zooC(1,:),'-b')
ylabel('zooc biomass')
ylim([3.15 4.85])
print('-dpng',[fpath 'testcase_zooC_loc1.png'])

figure(6)
plot(time,zoo_mort(1,:),'-b')
ylabel('zoo mort')
ylim([0.7 1.7])
print('-dpng',[fpath 'testcase_zoomort_loc1.png'])




