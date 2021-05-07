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

%% FEISTY ---------------------------------------------------------
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

% anomaly time series
load([fpath 'feisty_pi400_biomes_anom300_ln.mat'])
nt = length(yid);

%% reshape to save only land-free cells
% and transpose

pmat = nan*ones(3600,4,12);
fmat = nan*ones(3600,4,12);
dmat = nan*ones(3600,4,12);

% Exclude 1st 100 yrs (spinup for FEISTY)
pmat(:,:,1) = tp_anom(:,1201:end)';
pmat(:,:,2) = mz_anom(:,1201:end)';
pmat(:,:,3) = lz_anom(:,1201:end)';
pmat(:,:,4) = z_anom(:,1201:end)';
% xmat(:,:,7) = hpmz_anom(:,1201:end)';
% xmat(:,:,8) = hplz_anom(:,1201:end)';
% xmat(:,:,9) = hp_anom(:,1201:end)';

dmat(:,:,1) = tb_anom(:,1201:end)';
dmat(:,:,2) = det_anom(:,1201:end)';

% clear tp_anom tb_anom mz_anom lz_anom hpmz_anom hplz_anom det_anom
% clear z_anom hp_anom

%% already removed 1st 100 yrs (spinup)
%transpose
fmat(:,:,5) = sf_anom';
fmat(:,:,6) = mf_anom';
fmat(:,:,8)  = F_anom';

pmat(:,:,5) = sp_anom';
pmat(:,:,6) = mp_anom';
pmat(:,:,7) = lp_anom';
pmat(:,:,8)  = P_anom';

dmat(:,:,3)  = B_anom';
dmat(:,:,5) = sd_anom';
dmat(:,:,6) = md_anom';
dmat(:,:,7) = ld_anom';
dmat(:,:,8)  = D_anom';

pmat(:,:,9)  = S_anom';
pmat(:,:,10)  = M_anom';
pmat(:,:,11)  = L_anom';
pmat(:,:,12) = all_anom';

%clear sf_anom sp_anom sd_anom mf_anom mp_anom md_anom lp_anom ld_anom

%% FFT all at once
Fmat = NaN*ones(6,4,12);
Pmat = NaN*ones(6,4,12);
Dmat = NaN*ones(6,4,12);
FFTmatF = NaN*ones(nt/2,4,25);
FFTmatP = NaN*ones(nt/2,4,25);
FFTmatD = NaN*ones(nt/2,4,25);

%% P
for n = 1:12
    xts = squeeze(pmat(:,:,n));
    if(~isnan(xts(1,1)))
        Y = fft(xts);
        L = length(xts);
        P2 = abs(Y./L);
        P1 = P2(1:L/2+1 , :);
        P1 = P1(2:end , :);
        f1 = (0:(L/2))/L;
        f1 = f1(2:end);
        time = 1./f1;
        yr = time/12;
        dyr = diff(yr)';
        dyr(1800) = -9.254e-05;
        df = diff(f1)';
        df(1800) = df(1799);
        
        FFTmatP(:,:,n) = P1;
        
        %0-2 yrs
        y02 = find(yr<2);
        Pmat(1,:,n) = trapz(f1(y02),(P1(y02,:).*df(y02))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %2-5 yrs
        y25 = find(yr>=2 & yr<5);
        Pmat(2,:,n) = trapz(f1(y25),(P1(y25,:).*df(y25))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %5-12 yrs
        y512 = find(yr>=5 & yr<12);
        Pmat(3,:,n) = trapz(f1(y512),(P1(y512,:).*df(y512))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %12-20 yrs
        y1220 = find(yr>=12 & yr<20);
        Pmat(4,:,n) = trapz(f1(y1220),(P1(y1220,:).*df(y1220))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %20-30 yrs
        y2030 = find(yr>=20 & yr<30);
        Pmat(5,:,n) = trapz(f1(y2030),(P1(y2030,:).*df(y2030))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %30-50 yrs
        y3050 = find(yr>=30 & yr<50);
        Pmat(6,:,n) = trapz(f1(y3050),(P1(y3050,:).*df(y3050))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
    end
    
end

%% F
for n = 1:12
    xts = squeeze(fmat(:,:,n));
    if(~isnan(xts(1,1)))
        Y = fft(xts);
        L = length(xts);
        P2 = abs(Y./L);
        P1 = P2(1:L/2+1 , :);
        P1 = P1(2:end , :);
        f1 = (0:(L/2))/L;
        f1 = f1(2:end);
        time = 1./f1;
        yr = time/12;
        dyr = diff(yr)';
        dyr(1800) = -9.254e-05;
        df = diff(f1)';
        df(1800) = df(1799);
        
        FFTmatF(:,:,n) = P1;
        
        %0-2 yrs
        y02 = find(yr<2);
        Fmat(1,:,n) = trapz(f1(y02),(P1(y02,:).*df(y02))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %2-5 yrs
        y25 = find(yr>=2 & yr<5);
        Fmat(2,:,n) = trapz(f1(y25),(P1(y25,:).*df(y25))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %5-12 yrs
        y512 = find(yr>=5 & yr<12);
        Fmat(3,:,n) = trapz(f1(y512),(P1(y512,:).*df(y512))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %12-20 yrs
        y1220 = find(yr>=12 & yr<20);
        Fmat(4,:,n) = trapz(f1(y1220),(P1(y1220,:).*df(y1220))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %20-30 yrs
        y2030 = find(yr>=20 & yr<30);
        Fmat(5,:,n) = trapz(f1(y2030),(P1(y2030,:).*df(y2030))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %30-50 yrs
        y3050 = find(yr>=30 & yr<50);
        Fmat(6,:,n) = trapz(f1(y3050),(P1(y3050,:).*df(y3050))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
    end
    
end

%% D
for n = 1:12
    xts = squeeze(dmat(:,:,n));
    if(~isnan(xts(1,1)))
        Y = fft(xts);
        L = length(xts);
        P2 = abs(Y./L);
        P1 = P2(1:L/2+1 , :);
        P1 = P1(2:end , :);
        f1 = (0:(L/2))/L;
        f1 = f1(2:end);
        time = 1./f1;
        yr = time/12;
        dyr = diff(yr)';
        dyr(1800) = -9.254e-05;
        df = diff(f1)';
        df(1800) = df(1799);
        
        FFTmatD(:,:,n) = P1;
        
        %0-2 yrs
        y02 = find(yr<2);
        Dmat(1,:,n) = trapz(f1(y02),(P1(y02,:).*df(y02))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %2-5 yrs
        y25 = find(yr>=2 & yr<5);
        Dmat(2,:,n) = trapz(f1(y25),(P1(y25,:).*df(y25))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %5-12 yrs
        y512 = find(yr>=5 & yr<12);
        Dmat(3,:,n) = trapz(f1(y512),(P1(y512,:).*df(y512))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %12-20 yrs
        y1220 = find(yr>=12 & yr<20);
        Dmat(4,:,n) = trapz(f1(y1220),(P1(y1220,:).*df(y1220))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %20-30 yrs
        y2030 = find(yr>=20 & yr<30);
        Dmat(5,:,n) = trapz(f1(y2030),(P1(y2030,:).*df(y2030))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
        
        %30-50 yrs
        y3050 = find(yr>=30 & yr<50);
        Dmat(6,:,n) = trapz(f1(y3050),(P1(y3050,:).*df(y3050))) ./ trapz(f1,(P1.*df))  * (1./(2-0));
    end
    
end

%%
tchunk = {'1-2y','2-5y','5-12y','12-20y','20-30y','30-50y'};
Pvars = {'tp','mz','lz','z','sp','mp','lp','P','S','M','L','All'};
Fvars = {'',   '',  '', '','sf','mf', '', 'F', '', '','',''};
Dvars = {'tb','det','B','','sd','md','ld','D','', '','',''};
vars = {'t','M2^o','L2^o','2^o','s','m','l','all','S','M','L','All'};

save([fpath 'fft_trapz_pi400_cobalt_fesity_biomes_300yr_ln_v2.mat'],...
    'biome4_hist','biome','yr','time','f1','Fmat','Pmat','Dmat',...
    'fmat','pmat','dmat','FFTmatF','FFTmatP','FFTmatD',...
    'tchunk','Fvars','Pvars','Dvars','vars');

%% Figures
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
ppath = [pp cfile '/'];

%% Bubble chart?
tc = 1:3:18;
nvars = 1:12;
[vmat,tmat] = meshgrid(nvars,tc);

figure(1) %WHY DOESN'T THIS LOOK LIKE 1ST COL OF FIG 5?
for t=1:6
    if (t==1)
        bubblechart(tmat(t,:)-2,vmat(t,:),squeeze(Fmat(t,1,:)),'r'); hold on
        bubblechart(tmat(t,:)-1,vmat(t,:),squeeze(Pmat(t,1,:)),'b'); hold on
        bubblechart(tmat(t,:),vmat(t,:),squeeze(Dmat(t,1,:)),'#149414'); hold on
    elseif(t==2)
        bubblechart(tmat(t,:)-0.75,vmat(t,:),squeeze(Fmat(t,1,:)),'r'); hold on
        bubblechart(tmat(t,:),vmat(t,:),squeeze(Pmat(t,1,:)),'b'); hold on
        bubblechart(tmat(t,:)+0.75,vmat(t,:),squeeze(Dmat(t,1,:)),'#149414'); hold on
    else
        bubblechart(tmat(t,:)-0.5,vmat(t,:),squeeze(Fmat(t,1,:)),'r'); hold on
        bubblechart(tmat(t,:),vmat(t,:),squeeze(Pmat(t,1,:)),'b'); hold on
        bubblechart(tmat(t,:)+0.5,vmat(t,:),squeeze(Dmat(t,1,:)),'#149414'); hold on
    end
end
set(gca,'YTick',1:12,'YTickLabel',vars,'XTick',[0,4:3:18],'XTickLabel',tchunk)
title([biome{1} ' relative importance of variability scales'])
print('-dpng',[ppath 'Pre300_fft_trapz_bubble_biome1_v2.png'])

%%
figure(2)
for t=1:6
    if (t==1)
        bubblechart(tmat(t,:)-2,vmat(t,:),squeeze(Fmat(t,2,:)),'r'); hold on
        bubblechart(tmat(t,:)-1,vmat(t,:),squeeze(Pmat(t,2,:)),'b'); hold on
        bubblechart(tmat(t,:),vmat(t,:),squeeze(Dmat(t,2,:)),'#149414'); hold on
    elseif(t==2)
        bubblechart(tmat(t,:)-0.75,vmat(t,:),squeeze(Fmat(t,2,:)),'r'); hold on
        bubblechart(tmat(t,:),vmat(t,:),squeeze(Pmat(t,2,:)),'b'); hold on
        bubblechart(tmat(t,:)+0.75,vmat(t,:),squeeze(Dmat(t,2,:)),'#149414'); hold on
    else
        bubblechart(tmat(t,:)-0.5,vmat(t,:),squeeze(Fmat(t,2,:)),'r'); hold on
        bubblechart(tmat(t,:),vmat(t,:),squeeze(Pmat(t,2,:)),'b'); hold on
        bubblechart(tmat(t,:)+0.5,vmat(t,:),squeeze(Dmat(t,2,:)),'#149414'); hold on
    end
end
set(gca,'YTick',1:12,'YTickLabel',vars,'XTick',[0,4:3:18],'XTickLabel',tchunk)
title([biome{2} ' relative importance of variability scales'])
print('-dpng',[ppath 'Pre300_biomes_bubble_biome2_v2.png'])

figure(3)
for t=1:6
    if (t==1)
        bubblechart(tmat(t,:)-2,vmat(t,:),squeeze(Fmat(t,3,:)),'r'); hold on
        bubblechart(tmat(t,:)-1,vmat(t,:),squeeze(Pmat(t,3,:)),'b'); hold on
        bubblechart(tmat(t,:),vmat(t,:),squeeze(Dmat(t,3,:)),'#149414'); hold on
    elseif(t==2)
        bubblechart(tmat(t,:)-0.75,vmat(t,:),squeeze(Fmat(t,3,:)),'r'); hold on
        bubblechart(tmat(t,:),vmat(t,:),squeeze(Pmat(t,3,:)),'b'); hold on
        bubblechart(tmat(t,:)+0.75,vmat(t,:),squeeze(Dmat(t,3,:)),'#149414'); hold on
    else
        bubblechart(tmat(t,:)-0.5,vmat(t,:),squeeze(Fmat(t,3,:)),'r'); hold on
        bubblechart(tmat(t,:),vmat(t,:),squeeze(Pmat(t,3,:)),'b'); hold on
        bubblechart(tmat(t,:)+0.5,vmat(t,:),squeeze(Dmat(t,3,:)),'#149414'); hold on
    end
end
set(gca,'YTick',1:12,'YTickLabel',vars,'XTick',[0,4:3:18],'XTickLabel',tchunk)
title([biome{3} ' relative importance of variability scales'])
print('-dpng',[ppath 'Pre300_fft_trapz_bubble_biome3_v2.png'])

figure(4)
for t=1:6
    if (t==1)
        bubblechart(tmat(t,:)-2,vmat(t,:),squeeze(Fmat(t,4,:)),'r'); hold on
        bubblechart(tmat(t,:)-1,vmat(t,:),squeeze(Pmat(t,4,:)),'b'); hold on
        bubblechart(tmat(t,:),vmat(t,:),squeeze(Dmat(t,4,:)),'#149414'); hold on
    elseif(t==2)
        bubblechart(tmat(t,:)-0.75,vmat(t,:),squeeze(Fmat(t,4,:)),'r'); hold on
        bubblechart(tmat(t,:),vmat(t,:),squeeze(Pmat(t,4,:)),'b'); hold on
        bubblechart(tmat(t,:)+0.75,vmat(t,:),squeeze(Dmat(t,4,:)),'#149414'); hold on
    else
        bubblechart(tmat(t,:)-0.5,vmat(t,:),squeeze(Fmat(t,4,:)),'r'); hold on
        bubblechart(tmat(t,:),vmat(t,:),squeeze(Pmat(t,4,:)),'b'); hold on
        bubblechart(tmat(t,:)+0.5,vmat(t,:),squeeze(Dmat(t,4,:)),'#149414'); hold on
    end
end
set(gca,'YTick',1:12,'YTickLabel',vars,'XTick',[0,4:3:18],'XTickLabel',tchunk)
title([biome{4} ' relative importance of variability scales'])
print('-dpng',[ppath 'Pre300_fft_trapz_bubble_biome4_v2.png'])

%% Comp biomes by time
figure(5)
for b=1:4
    bubblechart(repmat(b,12,1)-0.25,[1:12]',squeeze(Fmat(1,b,:)),'r'); hold on
    bubblechart(repmat(b,12,1),[1:12]',squeeze(Pmat(1,b,:)),'b'); hold on
    bubblechart(repmat(b,12,1)+0.25,[1:12]',squeeze(Dmat(1,b,:)),'#149414'); hold on
end
set(gca,'YTick',1:12,'YTickLabel',vars,'XTick',1:4,'XTickLabel',biome)
title(['Relative importance of variability at ',tchunk{1} ,' scales'])
%print('-dpng',[ppath 'Pre300_fft_trapz_biomes_bubble_time1.png'])

figure(6)
for b=1:4
    bubblechart(repmat(b,12,1)-0.25,[1:12]',squeeze(Fmat(2,b,:)),'r'); hold on
    bubblechart(repmat(b,12,1),[1:12]',squeeze(Pmat(2,b,:)),'b'); hold on
    bubblechart(repmat(b,12,1)+0.25,[1:12]',squeeze(Dmat(2,b,:)),'#149414'); hold on
end
set(gca,'YTick',1:12,'YTickLabel',vars,'XTick',1:4,'XTickLabel',biome)
title(['Relative importance of variability at ',tchunk{2} ,' scales'])
%print('-dpng',[ppath 'Pre300_fft_trapz_biomes_bubble_time2.png'])

figure(7)
for b=1:4
    bubblechart(repmat(b,12,1)-0.25,[1:12]',squeeze(Fmat(3,b,:)),'r'); hold on
    bubblechart(repmat(b,12,1),[1:12]',squeeze(Pmat(3,b,:)),'b'); hold on
    bubblechart(repmat(b,12,1)+0.25,[1:12]',squeeze(Dmat(3,b,:)),'#149414'); hold on
end
set(gca,'YTick',1:12,'YTickLabel',vars,'XTick',1:4,'XTickLabel',biome)
title(['Relative importance of variability at ',tchunk{3} ,' scales'])
%print('-dpng',[ppath 'Pre300_fft_trapz_biomes_bubble_time3.png'])

figure(8)
for b=1:4
    bubblechart(repmat(b,12,1)-0.25,[1:12]',squeeze(Fmat(4,b,:)),'r'); hold on
    bubblechart(repmat(b,12,1),[1:12]',squeeze(Pmat(4,b,:)),'b'); hold on
    bubblechart(repmat(b,12,1)+0.25,[1:12]',squeeze(Dmat(4,b,:)),'#149414'); hold on
end
set(gca,'YTick',1:12,'YTickLabel',vars,'XTick',1:4,'XTickLabel',biome)
title(['Relative importance of variability at ',tchunk{4} ,' scales'])
%print('-dpng',[ppath 'Pre300_fft_trapz_biomes_bubble_time4.png'])

figure(9)
for b=1:4
    bubblechart(repmat(b,12,1)-0.25,[1:12]',squeeze(Fmat(5,b,:)),'r'); hold on
    bubblechart(repmat(b,12,1),[1:12]',squeeze(Pmat(5,b,:)),'b'); hold on
    bubblechart(repmat(b,12,1)+0.25,[1:12]',squeeze(Dmat(5,b,:)),'#149414'); hold on
end
set(gca,'YTick',1:12,'YTickLabel',vars,'XTick',1:4,'XTickLabel',biome)
title(['Relative importance of variability at ',tchunk{5} ,' scales'])
%print('-dpng',[ppath 'Pre300_fft_trapz_biomes_bubble_time5.png'])

figure(10)
for b=1:4
    bubblechart(repmat(b,12,1)-0.25,[1:12]',squeeze(Fmat(6,b,:)),'r'); hold on
    bubblechart(repmat(b,12,1),[1:12]',squeeze(Pmat(6,b,:)),'b'); hold on
    bubblechart(repmat(b,12,1)+0.25,[1:12]',squeeze(Dmat(6,b,:)),'#149414'); hold on
end
set(gca,'YTick',1:12,'YTickLabel',vars,'XTick',1:4,'XTickLabel',biome)
title(['Relative importance of variability at ',tchunk{6} ,' scales'])
%print('-dpng',[ppath 'Pre300_fft_trapz_biomes_bubble_time6.png'])
