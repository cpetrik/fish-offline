% Calc FFT from COBALT & FEISTY ts for diff var
% Estimate rel impt of variability of diff time period lengths
% c.f. LeMezo

clear all
close all

%% COBALT --------------------------------------------------------
spath='/Volumes/MIP/GCM_DATA/ESM2M_PI/';

% anomaly time series
load([spath 'cobalt_pi400_biomes_temp_anom.mat']);
load([spath 'cobalt_pi400_biomes_det_anom_ln.mat']);
load([spath 'cobalt_pi400_biomes_zoo_anom_ln.mat']);
load([spath 'cobalt_pi400_biomes_hploss_anom_ln.mat']);

%% reshape to save only land-free cells
% and transpose
% Exclude 1st 100 yrs (spinup for FEISTY)
xmat(:,:,1) = tp_anom(:,1201:end)';
xmat(:,:,2) = tb_anom(:,1201:end)';
xmat(:,:,3) = det_anom(:,1201:end)';
xmat(:,:,4) = mz_anom(:,1201:end)';
xmat(:,:,5) = lz_anom(:,1201:end)';
xmat(:,:,6) = z_anom(:,1201:end)';
xmat(:,:,7) = hpmz_anom(:,1201:end)';
xmat(:,:,8) = hplz_anom(:,1201:end)';
xmat(:,:,9) = hp_anom(:,1201:end)';

% clear tp_anom tb_anom mz_anom lz_anom hpmz_anom hplz_anom det_anom
% clear z_anom hp_anom

%% FEISTY ---------------------------------------------------------
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

% anomaly time series
load([fpath 'feisty_pi400_biomes_anom300_ln.mat'])
nt = length(yid);

%% already removed 1st 100 yrs (spinup)
%transpose
xmat(:,:,10) = sf_anom';
xmat(:,:,11) = sp_anom';
xmat(:,:,12) = sd_anom';
xmat(:,:,13) = mf_anom';
xmat(:,:,14) = mp_anom'; 
xmat(:,:,15) = md_anom';
xmat(:,:,16) = lp_anom'; 
xmat(:,:,17) = ld_anom'; 

xmat(:,:,18)  = B_anom';
xmat(:,:,19)  = F_anom';
xmat(:,:,20)  = P_anom';
xmat(:,:,21)  = D_anom';
xmat(:,:,22)  = S_anom';
xmat(:,:,23)  = M_anom';
xmat(:,:,24)  = L_anom';
xmat(:,:,25) = all_anom';

%clear sf_anom sp_anom sd_anom mf_anom mp_anom md_anom lp_anom ld_anom

%% FFT all at once
Dmat = NaN*ones(6,4,25);
Fmat = NaN*ones(nt/2,4,25);

for n = 1:25
    xts = squeeze(xmat(:,:,n));
    Y = fft(xts);
    L = length(xts);
    P2 = abs(Y./L);
    P1 = P2(1:L/2+1 , :);
    P1 = P1(2:end , :);
    f1 = (0:(L/2))/L;
    time = 1./f1;
    time = time(2:end);
    yr = time/12;
    
    Fmat(:,:,n) = P1;
    
    %0-2 yrs
    y02 = (yr<2);
    Dmat(1,:,n) = sum(P1(y02,:)) ./ sum(P1);
    
    %2-5 yrs
    y25 = (yr>=2 & yr<5);
    Dmat(2,:,n) = sum(P1(y25)) ./ sum(P1) * (1./(5-2));
    
    %5-12 yrs
    y512 = (yr>=5 & yr<12);
    Dmat(3,:,n)  = sum(P1(y512)) ./ sum(P1) * (1./(12-5));
    
    %12-20 yrs
    y1220 = (yr>=12 & yr<20);
    Dmat(4,:,n)  = sum(P1(y1220)) ./ sum(P1) * (1./(20-12));
    
    %20-30 yrs
    y2030 = (yr>=20 & yr<30);
    Dmat(5,:,n)  = sum(P1(y2030)) ./ sum(P1) * (1./(30-20));
    
    %30-50 yrs
    y3050 = (yr>=30 & yr<50);
    Dmat(6,:,n)  = sum(P1(y3050)) ./ sum(P1) * (1./(50-30));
    
end

%%
tchunk = {'1-2y','2-5y','5-12y','12-20y','20-30y','30-50y'};
vars = {'tp','tb','det','mz','lz','z','hpmz','hplz','hp',...
    'sf','sp','sd','mf','mp','md','lp','ld',...
    'B','F','P','D','S','M','L','All'};

save([fpath 'fft_pi400_cobalt_fesity_biomes_300yr_ln.mat'],...
    'biome4_hist','biome','yr','time','f1','Fmat','Dmat','xmat',...
    'tchunk','vars');

%% Figures
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
ppath = [pp cfile '/'];

%% Bubble chart?
tc = 1:6;
nvars = 1:25;
[vmat,tmat] = meshgrid(nvars,tc);

% bubblechart(tmat,vmat,squeeze(Dmat(:,1,:)))

figure(1)
for t=1:6
    bubblechart(tmat(t,:),vmat(t,:),squeeze(Dmat(t,1,:))); hold on
end
set(gca,'YTick',1:25,'YTickLabel',vars,'XTick',1:6,'XTickLabel',tchunk)
title([biome{1} ' relative importance of variability scales'])
print('-dpng',[ppath 'Pre300_biomes_bubble_biome1.png'])

%%
figure(2)
for t=1:6
    bubblechart(tmat(t,:),vmat(t,:),squeeze(Dmat(t,2,:))); hold on
end
set(gca,'YTick',1:25,'YTickLabel',vars,'XTick',1:6,'XTickLabel',tchunk)
title([biome{2} ' relative importance of variability scales'])
print('-dpng',[ppath 'Pre300_biomes_bubble_biome2.png'])

figure(3)
for t=1:6
    bubblechart(tmat(t,:),vmat(t,:),squeeze(Dmat(t,3,:))); hold on
end
set(gca,'YTick',1:25,'YTickLabel',vars,'XTick',1:6,'XTickLabel',tchunk)
title([biome{3} ' relative importance of variability scales'])
print('-dpng',[ppath 'Pre300_biomes_bubble_biome3.png'])

figure(4)
for t=1:6
    bubblechart(tmat(t,:),vmat(t,:),squeeze(Dmat(t,4,:))); hold on
end
set(gca,'YTick',1:25,'YTickLabel',vars,'XTick',1:6,'XTickLabel',tchunk)
title([biome{4} ' relative importance of variability scales'])
print('-dpng',[ppath 'Pre300_biomes_bubble_biome4.png'])

%% Comp biomes by time
figure(5)
for b=1:4
    bubblechart(repmat(b,25,1),[1:25]',squeeze(Dmat(1,b,:))); hold on
end
set(gca,'YTick',1:25,'YTickLabel',vars,'XTick',1:4,'XTickLabel',biome)
title(['Relative importance of variability at ',tchunk{1} ,' scales'])
print('-dpng',[ppath 'Pre300_biomes_bubble_time1.png'])

figure(6)
for b=1:4
    bubblechart(repmat(b,25,1),[1:25]',squeeze(Dmat(2,b,:))); hold on
end
set(gca,'YTick',1:25,'YTickLabel',vars,'XTick',1:4,'XTickLabel',biome)
title(['Relative importance of variability at ',tchunk{2} ,' scales'])
print('-dpng',[ppath 'Pre300_biomes_bubble_time2.png'])

figure(7)
for b=1:4
    bubblechart(repmat(b,25,1),[1:25]',squeeze(Dmat(3,b,:))); hold on
end
set(gca,'YTick',1:25,'YTickLabel',vars,'XTick',1:4,'XTickLabel',biome)
title(['Relative importance of variability at ',tchunk{3} ,' scales'])
print('-dpng',[ppath 'Pre300_biomes_bubble_time3.png'])

figure(8)
for b=1:4
    bubblechart(repmat(b,25,1),[1:25]',squeeze(Dmat(4,b,:))); hold on
end
set(gca,'YTick',1:25,'YTickLabel',vars,'XTick',1:4,'XTickLabel',biome)
title(['Relative importance of variability at ',tchunk{4} ,' scales'])
print('-dpng',[ppath 'Pre300_biomes_bubble_time4.png'])

figure(9)
for b=1:4
    bubblechart(repmat(b,25,1),[1:25]',squeeze(Dmat(5,b,:))); hold on
end
set(gca,'YTick',1:25,'YTickLabel',vars,'XTick',1:4,'XTickLabel',biome)
title(['Relative importance of variability at ',tchunk{5} ,' scales'])
print('-dpng',[ppath 'Pre300_biomes_bubble_time5.png'])

figure(10)
for b=1:4
    bubblechart(repmat(b,25,1),[1:25]',squeeze(Dmat(6,b,:))); hold on
end
set(gca,'YTick',1:25,'YTickLabel',vars,'XTick',1:4,'XTickLabel',biome)
title(['Relative importance of variability at ',tchunk{6} ,' scales'])
print('-dpng',[ppath 'Pre300_biomes_bubble_time6.png'])
