% Function called to calculate FOSI climatology
% To add to DPLE drift-corrected anomalies

function [clim_Tp,clim_Tb,clim_POC,clim_zooC,clim_loss] = fosi_clim_for_dple(Cdir)

% Load
load([Cdir 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'TEMP_150m','TEMP_bottom','POC_FLUX_IN_bottom','TLAT','TLONG','TAREA','time');
load([Cdir 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo_totloss_allphytoC.mat'],...
    'LzooC_150m','Lzoo_loss_150m');

% Doubles and nans
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;
LzooC_150m(LzooC_150m >= 9.9e+36) = nan;
Lzoo_loss_150m(Lzoo_loss_150m >= 9.9e+36) = nan;

fosi_Tp = double(TEMP_150m);
fosi_Tb = double(TEMP_bottom);
fosi_POC = double(POC_FLUX_IN_bottom);
fosi_zooC = double(LzooC_150m);
fosi_loss = double(Lzoo_loss_150m);
fosi_time = double(time);

clear TEMP_150m TEMP_bottom POC_FLUX_IN_bottom time LzooC_150m Lzoo_loss_150m

%% Monthly clim
[ni,nj,nmo] = size(fosi_Tp);
nyr = nmo/12;

clim_Tp = nan*ones(ni,nj,12);
clim_Tb = nan*ones(ni,nj,12);
clim_POC = nan*ones(ni,nj,12);
clim_zooC = nan*ones(ni,nj,12);
clim_loss = nan*ones(ni,nj,12);
for m = 1:12
    mo = m:12:nyr;
    clim_Tp(:,:,m) = nanmean(double(fosi_Tp(:,:,mo)),3);
    clim_Tb(:,:,m) = nanmean(double(fosi_Tb(:,:,mo)),3);
    clim_POC(:,:,m) = nanmean(double(fosi_POC(:,:,mo)),3);
    clim_zooC(:,:,m) = nanmean(double(fosi_zooC(:,:,mo)),3);
    clim_loss(:,:,m) = nanmean(double(fosi_loss(:,:,mo)),3);
end

%% Repeat climatol for 10 yrs
clim_Tp = repmat(clim_Tp,1,1,10);
clim_Tb = repmat(clim_Tb,1,1,10);
clim_POC = repmat(clim_POC,1,1,10);
clim_zooC = repmat(clim_zooC,1,1,10);
clim_loss = repmat(clim_loss,1,1,10);

whos clim_POC fosi_POC

%% Add Nov & Dec to beginning
Tp_ND = cat(3,clim_Tp(:,:,11),clim_Tp(:,:,12));
Tb_ND = cat(3,clim_Tb(:,:,11),clim_Tb(:,:,12));
POC_ND = cat(3,clim_POC(:,:,11),clim_POC(:,:,12));
zoo_ND = cat(3,clim_zooC(:,:,11),clim_zooC(:,:,12));
loss_ND = cat(3,clim_loss(:,:,11),clim_loss(:,:,12));

clim_Tp = cat(3,Tp_ND,clim_Tp);
clim_Tb = cat(3,Tb_ND,clim_Tb);
clim_POC = cat(3,POC_ND,clim_POC);
clim_zooC = cat(3,zoo_ND,clim_zooC);
clim_loss = cat(3,loss_ND,clim_loss);

whos clim_POC

%% Check climatol
%figure
%pcolor(squeeze(clim_zooC(:,:,3))); shading flat; colorbar;
%caxis([0 2e4])
%figure
%pcolor(squeeze(clim_zooC(:,:,9))); shading flat; colorbar;
%caxis([0 2e4])

figure
pcolor(squeeze(clim_zooC(:,:,17))); shading flat; colorbar;
caxis([0 2e4])
figure
pcolor(squeeze(clim_zooC(:,:,23))); shading flat; colorbar;
caxis([0 2e4])

clear fosi_Tp fosi_Tb fosi_POC fosi_zooC fosi_loss

end
