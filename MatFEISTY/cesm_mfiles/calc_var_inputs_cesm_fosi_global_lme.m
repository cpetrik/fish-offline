% CESM FOSI output
% calc interann variability by grid cell and lme

clear all
close all

%% Paths
fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6.mat'],'GRD');
load([fpath 'LME-mask-POP_gx1v6.mat']);

tlme = double(lme_mask);
tlme(tlme<0) = nan;

%% FEISTY Inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'TLAT','TLONG','TAREA','time','yr');
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo.mat'],...
    'LzooC_150m','Lzoo_loss_150m');

%% nans & zeros
TEMP_150m = double(TEMP_150m);
TEMP_bottom = double(TEMP_bottom);
POC_FLUX_IN_bottom = double(POC_FLUX_IN_bottom);

TEMP_bottom(TEMP_bottom >= 9.9e+36) = nan;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;

LzooC_150m(LzooC_150m<0) = 0.0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0.0;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom<0) = 0.0;

%% Units
%poc flux: mmol/m^3 cm/s
%zoo loss: mmol/m^3/s cm
%zoo: mmolC/m^3 cm
%tp: degC
%tb: degC

% meso zoo: nmolC cm-2 to g(WW) m-2
% 1e9 nmol in 1 mol C
% 1e4 cm2 in 1 m2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
LzooC_150m = LzooC_150m * 1e-9 * 1e4 * 12.01 * 9.0;

% meso zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
% detrital flux to benthos: nmolC cm-2 s-1 to g(WW) m-2 d-1
% 1e9 nmol in 1 mol C
% 1e4 cm2 in 1 m2
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
Lzoo_loss_150m = Lzoo_loss_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
POC_FLUX_IN_bottom = POC_FLUX_IN_bottom * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

%% annual means
nt = length(time);
st=1:12:nt;
en=12:12:nt;
nyr = nt/12;
[ni,nj] = size(TLONG);

tp = nan*ones(ni,nj,nyr);
tb = nan*ones(ni,nj,nyr);
det = nan*ones(ni,nj,nyr);
zoo = nan*ones(ni,nj,nyr);
zlos = nan*ones(ni,nj,nyr);
for n=1:nyr
    tp(:,:,n)=nanmean(TEMP_150m(:,:,st(n):en(n)),3);
    tb(:,:,n)=nanmean(TEMP_bottom(:,:,st(n):en(n)),3);
    det(:,:,n)=nanmean(POC_FLUX_IN_bottom(:,:,st(n):en(n)),3);
    zoo(:,:,n)=nanmean(LzooC_150m(:,:,st(n):en(n)),3);
    zlos(:,:,n)=nanmean(Lzoo_loss_150m(:,:,st(n):en(n)),3);
end

%% anomalies
atp = tp - nanmean(tp,3);
atb = tb - nanmean(tb,3);
adet = det - nanmean(det,3);
azoo = zoo - nanmean(zoo,3);
azlos = zlos - nanmean(zlos,3);

save([fpath 'CESM_FOSI_interann_mean_forcings_anom.mat'],...
    'tp','tb','det','zoo','zlos',...
    'atp','atb','adet','azoo','azlos');


%% var by grid cell
vtp = var(atp,0,3,'omitnan');
vtb = var(atb,0,3,'omitnan');
vdet = var(adet,0,3,'omitnan');
vzoo = var(azoo,0,3,'omitnan');
vzlos = var(azlos,0,3,'omitnan');

%% var by lme
atp2 = reshape(atp,ni*nj,nyr);
atb2 = reshape(atb,ni*nj,nyr);
adet2 = reshape(adet,ni*nj,nyr);
amz = reshape(azoo,ni*nj,nyr);
amzl = reshape(azlos,ni*nj,nyr);

lme_tp_var1 = NaN*ones(66,1);
lme_tb_var1 = NaN*ones(66,1);
lme_det_var1 = NaN*ones(66,1);
lme_mz_var1 = NaN*ones(66,1);
lme_mzl_var1 = NaN*ones(66,1);
lme_tp_var2 = NaN*ones(66,1);
lme_tb_var2 = NaN*ones(66,1);
lme_det_var2 = NaN*ones(66,1);
lme_mz_var2 = NaN*ones(66,1);
lme_mzl_var2 = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    
    ltp = atp2(lid,:);
    ltb = atb2(lid,:);
    ldet = adet2(lid,:);
    lmz = amz(lid,:);
    lmzl = amzl(lid,:);
    
    lme_tp_var1(L,1) = var(ltp(:),'omitnan');
    lme_tb_var1(L,1) = var(ltb(:),'omitnan');
    lme_det_var1(L,1) = var(ldet(:),'omitnan');
    lme_mz_var1(L,1) = var(lmz(:),'omitnan');
    lme_mzl_var1(L,1) = var(lmzl(:),'omitnan');
    
    lme_tp_var2(L,1) = nanmean(var(ltp,0,2));
    lme_tb_var2(L,1) = nanmean(var(ltb,0,2));
    lme_det_var2(L,1) = nanmean(var(ldet,0,2));
    lme_mz_var2(L,1) = nanmean(var(lmz,0,2));
    lme_mzl_var2(L,1) = nanmean(var(lmzl,0,2));
    
end
% var1 == var2

%% map info
clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

cmYOR=cbrewer('seq','YlOrRd',66,'pchip');
cmYOB=cbrewer('seq','YlOrBr',66,'pchip');
cmOR=cbrewer('seq','OrRd',66,'pchip');

%% map
figure(1)
subplot('Position',[0.01 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,vtp)
colormap(cmYOR)
caxis([0 1])
colorbar%('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'var Tp','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,vtb)
colormap(cmYOR)
caxis([0 0.3])
colorbar%('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'var Tb','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,vzoo)
colormap(cmYOR)
caxis([0 1.5])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'var Zoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,vdet)
colormap(cmYOR)
caxis([0 0.02])
colorbar%('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'var Det','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,vzlos)
colormap(cmYOR)
caxis([0 0.02])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'var Zoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CESM_FOSI_interann_var_forcings.png'])

%% save
save([fpath 'CESM_FOSI_interann_var_forcings.mat'],...
    'vtp','vtb','vdet','vzoo','vzlos',...
    'lme_tp_var1','lme_tb_var1','lme_det_var1','lme_mz_var1','lme_mzl_var1');

