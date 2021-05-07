% Correlations between variance of forcing and fishe

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
spath='/Volumes/MIP/GCM_DATA/CORE-forced/';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

%% 
load([spath 'cobalt_core_variance_1950_2007.mat'])
load([fpath 'fesity_core_variance_1950_2007.mat'])

%% Correlations
%r
r(1,1:16) = corr(var_vec,tp_100_var)';
r(2,:) = corr(var_vec,tb_var)';
r(3,:) = corr(var_vec,det_var)';
r(4,:) = corr(var_vec,mz_var)';
r(5,:) = corr(var_vec,lz_var)';
r(6,:) = corr(var_vec,mz_hp_var)';
r(7,:) = corr(var_vec,lz_hp_var)';

%% Figure of scatter plots
for n = 1:length(var_vec_tex)
clf
figure(1)
subplot(3,3,1)
scatter(tp_100_var,var_vec(:,n),10,tp_100_mean,'filled'); hold on;
cmocean('thermal');
colorbar('Position',[0.65 0.1 0.025 0.5],'orientation','vertical')
%text(0.5,-1.5,['r = ' sprintf('%2.2f',r(1,n))])
%axis([-3 3 -6 -1])
xlabel('TP')

subplot(3,3,2)
scatter(tb_var,var_vec(:,n),10,tp_100_mean,'filled'); hold on;
cmocean('thermal');
% text(-20,-1.5,['r = ' sprintf('%2.2f',r(2,n))])
% axis([-50 3 -6 -1])
xlabel('TB')
title(var_vec_tex{n})

subplot(3,3,3)
scatter(det_var,var_vec(:,n),10,tp_100_mean,'filled'); hold on;
cmocean('thermal');
% text(0.5,-1.5,['r = ' sprintf('%2.2f',r(3,n))])
% axis([-4 3 -6 -1])
xlabel('Det')

subplot(3,3,4)
scatter(mz_var,var_vec(:,n),10,tp_100_mean,'filled'); hold on;
cmocean('thermal');
% text(0,-1.5,['r = ' sprintf('%2.2f',r(4,n))])
% axis([-3 3 -6 -1])
xlabel('MZ')

subplot(3,3,5)
scatter(lz_var,var_vec(:,n),10,tp_100_mean,'filled'); hold on;
cmocean('thermal');
% text(-20,-1.5,['r = ' sprintf('%2.2f',r(5,n))])
% axis([-50 3 -6 -1])
xlabel('LZ')

subplot(3,3,7)
scatter(mz_hp_var,var_vec(:,n),10,tp_100_mean,'filled'); hold on;
cmocean('thermal');
% text(2,-1.5,['r = ' sprintf('%2.2f',r(6,n))])
% axis([0.5 3 -7 -1])
xlabel('MZ HPloss')

subplot(3,3,8)
scatter(lz_hp_var,var_vec(:,n),10,tp_100_mean,'filled'); hold on;
cmocean('thermal');
% text(-10,-1.5,['r = ' sprintf('%2.2f',r(7,n))])
% axis([-30 3 -6 -1])
xlabel('LZ HPloss')

stamp('')
print('-dpng',['Corr_var_CORE_' var_vec_tex{n} '.png'])
end


