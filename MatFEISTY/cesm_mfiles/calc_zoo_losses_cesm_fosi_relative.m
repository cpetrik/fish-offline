% CESM FOSI output
% Back calc zoo quad mort = loss - linear mort
% Compare fraction of total that is linear vs. quad

clear 
close all

%% Paths

fpath='/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
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

%% Units
%zoo loss: mmol/m^3/s cm
%zoo: mmolC/m^3 cm
%tp: degC

% meso zoo: nmolC cm-2 to g(WW) m-2
LzooC_150m = LzooC_150m * 1e-9 * 1e4 * 12.01 * 9.0;

% meso zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
Lzoo_loss_150m = Lzoo_loss_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

%% Temp fn
T = min(TEMP_150m(:)):max(TEMP_150m(:));
tfn = exp(-4000 .* ( (1./(T+273.15)) - (1./303.15) ));

figure(1)
plot(T,tfn);

Tfn = exp(-4000 .* ( (1./(TEMP_150m+273.15)) - (1./303.15) ));

%% Zprime
Zprime = max((LzooC_150m - 0.01),0);
Lzoo_lin_150m = (Tfn .* 0.1 .* Zprime);

%% Z quad
Lzoo_quad_150m = Lzoo_loss_150m - (Tfn .* 0.1 .* Zprime);
Lzoo_quad_150m = max(Lzoo_quad_150m,0);

%%
quad_tot = Lzoo_quad_150m ./ (Lzoo_loss_150m+eps); %mean 0.4058
quad_lin = Lzoo_quad_150m ./ (Lzoo_lin_150m+eps);  %mean 2.5230e+05
% times where quad is 1e7 times more than lin

%% space means
sZbio = nanmean(LzooC_150m,3);
sZtot = nanmean(Lzoo_loss_150m,3);
sZlin = nanmean(Lzoo_lin_150m,3);
sZquad = nanmean(Lzoo_quad_150m,3);
sZqt = nanmean(quad_tot,3);
sZql = nanmean(quad_lin,3);

%% time means
[ni,nj,nt] = size(Lzoo_loss_150m);
tZbio = nanmean(reshape(LzooC_150m,ni*nj,nt));
tZtot = nanmean(reshape(Lzoo_loss_150m,ni*nj,nt));
tZlin = nanmean(reshape(Lzoo_lin_150m,ni*nj,nt));
tZquad = nanmean(reshape(Lzoo_quad_150m,ni*nj,nt));
tZqt = nanmean(reshape(quad_tot,ni*nj,nt));
tZql = nanmean(reshape(quad_lin,ni*nj,nt));

%%
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

figure(4)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(sZtot))
cmocean('tempo')
caxis([-2 1])
colorbar%('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo tot loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(sZbio))
cmocean('tempo')
caxis([0 2])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo biom','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(sZlin))
cmocean('tempo')
caxis([-2 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo linear loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(sZquad))
cmocean('tempo')
caxis([-2 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'log_1_0 Zoo quad loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(sZql))
cmocean('dense')
caxis([0 3])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Quad / Linear','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(sZqt))
cmocean('matter')
caxis([0 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'Quad / Total','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CESM_FOSI_mean_zoo_loss_comp_quad_lin_tot.png'])

%% time series
figure(5)
subplot(2,1,1)
plot(yr,log10(tZtot),'k'); hold on;
plot(yr,log10(tZlin),'r'); hold on;
plot(yr,log10(tZquad),'b'); hold on;
plot(yr,log10(tZbio),'color',[0 0.7 0]); hold on;
legend('total','linear','quad','biom')
ylabel('log_1_0')

subplot(2,1,2)
plot(yr,(tZqt),'k'); hold on;
plot(yr,(tZql),'r'); hold on;
legend('quad/total','quad/linear')
xlabel('Time (yrs)')
ylim([0 1.5])


