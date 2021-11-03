% CESM FOSI output
% calc interann variability by grid cell and lme

clear all
close all

%% Paths
fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6_noSeas.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');
load([fpath 'LME-mask-POP_gx1v6.mat']);

tlme = double(lme_mask);
tlme(tlme<0) = nan;

%% FEISTY Inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'TLAT','TLONG','TAREA','time','yr');
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo_totloss_allphytoC.mat'],...
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

%% annual means to remove seasonal var
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

save([fpath 'CESM_FOSI_v13_interann_mean_forcings_anom.mat'],...
    'tp','tb','det','zoo','zlos',...
    'atp','atb','adet','azoo','azlos');

%% mean & std by grid cell
tp_mean = nanmean(tp,3);
tb_mean = nanmean(tb,3);
det_mean = nanmean(det,3);
mz_mean = nanmean(zoo,3);
mzl_mean = nanmean(zlos,3);
        
stp = std(tp,0,3,'omitnan');
stb = std(tb,0,3,'omitnan');
sdet = std(det,0,3,'omitnan');
szoo = std(zoo,0,3,'omitnan');
szlos = std(zlos,0,3,'omitnan');

%% std by lme
atp2 = reshape(tp,ni*nj,nyr);
atb2 = reshape(tb,ni*nj,nyr);
adet2 = reshape(det,ni*nj,nyr);
amz = reshape(zoo,ni*nj,nyr);
amzl = reshape(zlos,ni*nj,nyr);

lme_tp_mean = NaN*ones(66,1);
lme_tb_mean = NaN*ones(66,1);
lme_det_mean = NaN*ones(66,1);
lme_mz_mean = NaN*ones(66,1);
lme_mzl_mean = NaN*ones(66,1);
lme_tp_std1 = NaN*ones(66,1);
lme_tb_std1 = NaN*ones(66,1);
lme_det_std1 = NaN*ones(66,1);
lme_mz_std1 = NaN*ones(66,1);
lme_mzl_std1 = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    if ~isempty(lid)
        
        ltp = atp2(lid,:);
        ltb = atb2(lid,:);
        ldet = adet2(lid,:);
        lmz = amz(lid,:);
        lmzl = amzl(lid,:);
        
        lme_tp_mean(L,1) = nanmean(ltp(:));
        lme_tb_mean(L,1) = nanmean(ltb(:));
        lme_det_mean(L,1) = nanmean(ldet(:));
        lme_mz_mean(L,1) = nanmean(lmz(:));
        lme_mzl_mean(L,1) = nanmean(lmzl(:));
        
        lme_tp_std1(L,1) = std(ltp(:),'omitnan');
        lme_tb_std1(L,1) = std(ltb(:),'omitnan');
        lme_det_std1(L,1) = std(ldet(:),'omitnan');
        lme_mz_std1(L,1) = std(lmz(:),'omitnan');
        lme_mzl_std1(L,1) = std(lmzl(:),'omitnan');
        
    end
end

%% Coefficient of variance
cvtp = stp ./ tp_mean;
cvtb = stb ./ tb_mean;
cvdet = sdet ./ det_mean;
cvzoo = szoo ./ mz_mean;
cvzlos = szlos ./ mzl_mean;

lme_tp_cv = lme_tp_std1 ./ lme_tp_mean;
lme_tb_cv = lme_tb_std1 ./ lme_tb_mean;
lme_det_cv = lme_det_std1 ./ lme_det_mean;
lme_mz_cv = lme_mz_std1 ./ lme_mz_mean;
lme_mzl_cv = lme_mzl_std1 ./ lme_mzl_mean;

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
surfm(TLAT,TLONG,cvtp)
colormap(cmYOR)
caxis([0 0.5])
colorbar%('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'CV Tp','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.68 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cvtb)
colormap(cmYOR)
caxis([0 0.5])
colorbar%('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'CV Tb','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cvzoo)
colormap(cmYOR)
caxis([0 0.5])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'CV Zoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.37 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cvdet)
colormap(cmYOR)
caxis([0 0.5])
colorbar%('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'CV Det','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.06 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cvzlos)
colormap(cmYOR)
caxis([0 0.5])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
text(0.2,1.65,'CV Zoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CESM_FOSI_v13_interann_coeffvar_forcings.png'])

%% save
save([fpath 'CESM_FOSI_v13_interann_var_forcings.mat'],...
    'tp_mean','tb_mean','det_mean','mz_mean','mzl_mean',...
    'stp','stb','sdet','szoo','szlos',...
    'cvtp','cvtb','cvdet','cvzoo','cvzlos',...
    'lme_tp_std1','lme_tb_std1','lme_det_std1','lme_mz_std1','lme_mzl_std1',...
    'lme_tp_mean','lme_tb_mean','lme_det_mean','lme_mz_mean','lme_mzl_mean',...
    'lme_tp_cv','lme_tb_cv','lme_det_cv','lme_mz_cv','lme_mzl_cv');

