% CESM FOSI output
% calc interann variability by grid cell and lme

clear all
close all

%% Paths
%fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
fpath='/Volumes/petrik-lab/GCM_Data/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6_noSeas.mat'],'mask');
load([fpath 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');
load([fpath 'LME-mask-POP_gx1v6.mat']);

%% FEISTY Inputs
fpath='/Volumes/petrik-lab/GCM_Data/CESM/FOSI/';
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

%% remove linear trend
%vectorize
tp2 = reshape(tp,ni*nj,nyr);
tb2 = reshape(tb,ni*nj,nyr);
det2 = reshape(det,ni*nj,nyr);
mz = reshape(zoo,ni*nj,nyr);
mzl = reshape(zlos,ni*nj,nyr);

%just ocean cells
nid = length(GRD.ID);
tp2 = tp2(GRD.ID,:);
tb2 = tb2(GRD.ID,:);
det2 = det2(GRD.ID,:);
mz = mz(GRD.ID,:);
mzl = mzl(GRD.ID,:);

xtp = NaN*ones(size(tp2));
xtb = NaN*ones(size(tb2));
xdet = NaN*ones(size(det2));
xz = NaN*ones(size(mz));
xzl = NaN*ones(size(mzl));

for i = 1:nid
    % Drivers
    xi = tp2(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); %0.651784 seconds.
    tH = m*t + b;
    dR = R - tH;
    xtp(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = tb2(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    xtb(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = det2(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    xdet(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = mz(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    xz(i,:) = dR;
    clear R T t b m tH dR data
    
    xi = mzl(i,:);
    inan = ~isnan(xi);
    R = (xi(inan))';
    T = length(R);
    t = (1:T)';
    data = [t R];
    [m,b] = TheilSen(data); 
    tH = m*t + b;
    dR = R - tH;
    xzl(i,:) = dR;
    clear R T t b m tH dR data
    
end

%% anomalies
atp = xtp - nanmean(xtp,2);
atb = xtb - nanmean(xtb,2);
adet = xdet - nanmean(xdet,2);
azoo = xz - nanmean(xz,2);
azlos = xzl - nanmean(xzl,2);

%% std by lme before putting back on grid
tlme = double(lme_mask);
tlme(tlme<0) = nan;
olme = tlme(GRD.ID);

lme_tp_stda = NaN*ones(66,1);
lme_tb_stda = NaN*ones(66,1);
lme_det_stda = NaN*ones(66,1);
lme_mz_stda = NaN*ones(66,1);
lme_mzl_stda = NaN*ones(66,1);

for L=1:66
    lid = find(olme==L);
    if ~isempty(lid)

        ltp = atp(lid,:);
        ltb = atb(lid,:);
        ldet = adet(lid,:);
        lmz = azoo(lid,:);
        lmzl = azlos(lid,:);

        lme_tp_stda(L,1) = std(ltp(:),'omitnan');
        lme_tb_stda(L,1) = std(ltb(:),'omitnan');
        lme_det_stda(L,1) = std(ldet(:),'omitnan');
        lme_mz_stda(L,1) = std(lmz(:),'omitnan');
        lme_mzl_stda(L,1) = std(lmzl(:),'omitnan');

    end
end

%% put anoms on grid
Atp   = nan*ones(ni,nj,nyr);
Atb   = nan*ones(ni,nj,nyr);
Adet  = nan*ones(ni,nj,nyr);
Azoo  = nan*ones(ni,nj,nyr);
Azlos = nan*ones(ni,nj,nyr);

for t=1:nyr
    
    vv = atp(:,t);
    gv = nan*ones(ni,nj);
    gv(GRD.ID) = vv;
    Atp(:,:,t) = gv;
    clear vv gv
    
    vv = atb(:,t);
    gv = nan*ones(ni,nj);
    gv(GRD.ID) = vv;
    Atb(:,:,t) = gv;
    clear vv gv
    
    vv = adet(:,t);
    gv = nan*ones(ni,nj);
    gv(GRD.ID) = vv;
    Adet(:,:,t) = gv;
    clear vv gv
    
    vv = azoo(:,t);
    gv = nan*ones(ni,nj);
    gv(GRD.ID) = vv;
    Azoo(:,:,t) = gv;
    clear vv gv
    
    vv = azlos(:,t);
    gv = nan*ones(ni,nj);
    gv(GRD.ID) = vv;
    Azlos(:,:,t) = gv;
    clear vv gv
    
end

%% mean & std by grid cell
tp_mean = nanmean(tp,3);
tb_mean = nanmean(tb,3);
det_mean = nanmean(det,3);
mz_mean = nanmean(zoo,3);
mzl_mean = nanmean(zlos,3);

tp_std = std(tp,0,3,'omitnan');
tb_std = std(tb,0,3,'omitnan');
det_std = std(det,0,3,'omitnan');
zoo_std = std(zoo,0,3,'omitnan');
zlos_std = std(zlos,0,3,'omitnan');

tp_stda = std(atp,0,3,'omitnan');
tb_stda = std(atb,0,3,'omitnan');
det_stda = std(adet,0,3,'omitnan');
zoo_stda = std(azoo,0,3,'omitnan');
zlos_stda = std(azlos,0,3,'omitnan');

%% Coefficient of variance
cvtp = tp_std ./ tp_mean;
cvtb = tb_std ./ tb_mean;
cvdet = det_std ./ det_mean;
cvzoo = zoo_std ./ mz_mean;
cvzlos = zlos_std ./ mzl_mean;

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

print('-dpng',[pp 'Map_CESM_FOSI_v15_interann_coeffvar_forcings.png'])

%% save anoms
save([fpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat'],...
    'tp','tb','det','zoo','zlos',...
    'Atp','Atb','Adet','Azoo','Azlos');

%% save
save([fpath 'CESM_FOSI_v15_interann_var_forcings.mat'],...
    'tp_mean','tb_mean','det_mean','mz_mean','mzl_mean',...
    'tp_std','tb_std','det_std','zoo_std','zlos_std',...
    'tp_stda','tb_stda','det_stda','zoo_stda','zlos_stda',...
    'cvtp','cvtb','cvdet','cvzoo','cvzlos',...
    'lme_tp_stda','lme_tb_stda','lme_det_stda','lme_mz_stda','lme_mzl_stda');
