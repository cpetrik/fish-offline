% MARBL zoo loss fn
% in model orig units

clear all
close all

figp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

%% Paths

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6.mat'],'GRD');

%% FEISTY Inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units',...
    'TLAT','TLONG','TAREA','time','yr');
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo.mat'],...
    'LzooC_150m','Lzoo_loss_150m');

%% nans & zeros
TEMP_150m = double(TEMP_150m);
LzooC_150m(LzooC_150m<0) = 0.0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0.0;

%% Temp fn
T = -2:0.5:35;
Q10 = 2.0;
T0Kelv = 273.15;
Tref = 30.0;
tfn = Q10 .^ ( ( (T + T0Kelv) - (Tref + T0Kelv) ) / 10.0 );

figure(1)
plot(T,tfn,'LineWidth',2);
title('Tfunc')
xlabel('temp')
ylabel('mult on rate')
xlim([-2 35])
set(gca,'YTick',[0:0.25:1.5],'YTickLabel',[0:0.25:1.5])
print('-dpng',[figp 'zoo_loss_tfunc.png'])

%% f_loss_thres
Z = (0:0.5:300); %depth in cm
Z = Z*1e2;
thres_z1 = 100.0e2; 
thres_z2 = 200.0e2;

f_loss_thres = nan(length(Z),1);

for k=1:length(Z)
    z = Z(k);
if (z > thres_z1)
    if (z < thres_z2)
        f_loss_thres(k) = (thres_z2 - z) ./ (thres_z2 - thres_z1);
    else
        f_loss_thres(k) = 0.0;
    end
else
    f_loss_thres(k) = 1.0;
end
end

figure(2)
plot(f_loss_thres,(-1*Z*1e-2),'LineWidth',2);
xlim([-0.1 1.1])
title('f loss thres')
xlabel('f loss thres')
ylabel('depth (m)')
print('-dpng',[figp 'zoo_loss_f_loss_thres_by_depth.png'])

zid = find(Z==15000);
mean(f_loss_thres(1:zid))

avg_f_loss_thres = 0.9167;

%% Zprime
loss_thres_zoo = 0.2; 
C_loss_thres = avg_f_loss_thres * loss_thres_zoo;

zooC_150m = quantile(LzooC_150m(:),[0.01 0.05:0.05:0.95 0.99]);

Zprime = max((zooC_150m - C_loss_thres),0.0);

%% Z quad
spd = 86400; %seconds per day
dps = 1 ./ spd; %days per second
parm_z_mort = 0.08 * dps;
parm_z_mort2 = 0.42 * dps;
tfn = tfn';

zmort = parm_z_mort .* tfn;
zmort2 = parm_z_mort2 .* tfn;
zoo_loss = (zmort2 .* Zprime.^1.4) + (zmort .* Zprime);
zoo_quad = zoo_loss - (zmort .* Zprime);
zoo_lin = zoo_loss - (zmort2 .* Zprime.^1.4);

%%
figure(3)
subplot(2,2,1)
plot(T,zoo_quad(:,[1 11 21]),'LineWidth',2);
title('Effect of temp')
xlabel('temp')
ylabel('zoo quad loss')
xlim([-2 35])
legend({'low','mid','high'},'location','northwest')

subplot(2,2,2)
plot(zooC_150m,zoo_quad([1 38 75],:),'LineWidth',2);
title('Effect of biomass')
xlabel('biomass')
ylabel('zoo quad loss')
%xlim([-2 35])

subplot(2,2,3)
plot(T,zoo_loss(:,[1 11 21]),'LineWidth',2);
xlabel('temp')
ylabel('zoo tot loss')
xlim([-2 35])
legend({'low','mid','high'},'location','northwest')

subplot(2,2,4)
plot(zooC_150m,zoo_loss([1 38 75],:),'LineWidth',2); hold on;
xlabel('biomass')
ylabel('zoo tot loss')
print('-dpng',[figp 'zoo_quad_loss_theor.png'])

%%
figure(4)
subplot(2,2,1)
plot(T,zoo_quad(:,7),'r','LineWidth',2); hold on;
plot(T,zoo_loss(:,7),'k--','LineWidth',2); hold on;
plot(T,zoo_lin(:,7),'b','LineWidth',2);
title('Effect of temp')
xlabel('temp')
ylabel('biomass = 1.9e3')
xlim([-2 35])
legend({'quad','tot','linear'},'location','northwest')

subplot(2,2,2)
plot(zooC_150m,zoo_quad(25,:),'r','LineWidth',2); hold on;
plot(zooC_150m,zoo_loss(25,:),'k--','LineWidth',2); hold on;
plot(zooC_150m,zoo_lin(25,:),'b','LineWidth',2);
title('Effect of biomass')
xlabel('biomass')
ylabel('T = 10')
%legend({'low','mid','high'},'location','northwest')

subplot(2,2,3)
plot(T,zoo_quad(:,14),'r','LineWidth',2); hold on;
plot(T,zoo_loss(:,14),'k--','LineWidth',2); hold on;
plot(T,zoo_lin(:,14),'b','LineWidth',2); hold on;
xlabel('temp')
ylabel('biomass = 2.8e3')
xlim([-2 35])
%legend({'low','mid','high'},'location','northwest')

subplot(2,2,4)
plot(zooC_150m,zoo_quad(50,:),'r','LineWidth',2); hold on;
plot(zooC_150m,zoo_loss(50,:),'k--','LineWidth',2); hold on;
plot(zooC_150m,zoo_lin(50,:),'b','LineWidth',2); hold on;
xlabel('biomass')
ylabel('T = 22.5')
print('-dpng',[figp 'zoo_quad_loss_theor2.png'])

%% Apply

%Lzoo_quad = Lzoo_loss_150m - (tfn .* parm_z_mort .* Zprime);
%Lzoo_quad = max(Lzoo_quad_150m,0);

