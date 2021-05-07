% Calc mean and std of powerspectra from COBALT ts for diff var
% By 4 biomes, N&S hemispheres separate

clear all
close all

%% COBALT --------------------------------------------------------
spath='/Volumes/MIP/GCM_DATA/CORE-forced/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
exper = 'Biome_exper_btm_';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/' exper];
ppath = [pp cfile '/Biome_exper/'];

% power spectra
load([fpath 'powerspec_feisty_ln.mat'])

%% Each biome repeated 220 times
% 20 random draws x 11 noise slopes x 8 biomes = 1760

% xsf = reshape(sf,20,[],8);
% xsp = reshape(sp,20,[],8);
% xsd = reshape(sd,20,[],8);
% xmf = reshape(mf,20,[],8);
% xmp = reshape(mp,20,[],8);
% xmd = reshape(md,20,[],8);
% xlp = reshape(lp,20,[],8);
% xld = reshape(ld,20,[],8);
% xB  = reshape(B,20,[],8);
% xF  = reshape(F,20,[],8);
% xP  = reshape(P,20,[],8);
% xD  = reshape(D,20,[],8);
% xS  = reshape(S,20,[],8);
% xM  = reshape(M,20,[],8);
% xL  = reshape(L,20,[],8);
% xall = reshape(All,20,[],8);

%% Means
% msf = squeeze(mean(xsf,1));
% msp = squeeze(mean(xsp,1));
% msd = squeeze(mean(xsd,1));
% mmf = squeeze(mean(xmf,1));
% mmp = squeeze(mean(xmp,1));
% mmd = squeeze(mean(xmd,1));
% mlp = squeeze(mean(xlp,1));
% mld = squeeze(mean(xld,1));
% mB  = squeeze(mean(xB,1));
% mF  = squeeze(mean(xF,1));
% mP  = squeeze(mean(xP,1));
% mD  = squeeze(mean(xD,1));
% mS  = squeeze(mean(xS,1));
% mM  = squeeze(mean(xM,1));
% mL  = squeeze(mean(xL,1));
% mall = squeeze(mean(xall,1));

%% Std dev
% ssf = squeeze(std(xsf,1));
% ssp = squeeze(std(xsp,1));
% ssd = squeeze(std(xsd,1));
% smf = squeeze(std(xmf,1));
% smp = squeeze(std(xmp,1));
% smd = squeeze(std(xmd,1));
% slp = squeeze(std(xlp,1));
% sld = squeeze(std(xld,1));
% sB  = squeeze(std(xB,1));
% sF  = squeeze(std(xF,1));
% sP  = squeeze(std(xP,1));
% sD  = squeeze(std(xD,1));
% sS  = squeeze(std(xS,1));
% sM  = squeeze(std(xM,1));
% sL  = squeeze(std(xL,1));
% sall = squeeze(std(xall,1));

% Std error = std / sqrt(n); Necessary?

%% Try a diff way to make sure I didn't mess up reshaping
xsf = reshape(sf,220,8);
xsp = reshape(sp,220,8);
xsd = reshape(sd,220,8);
xmf = reshape(mf,220,8);
xmp = reshape(mp,220,8);
xmd = reshape(md,220,8);
xlp = reshape(lp,220,8);
xld = reshape(ld,220,8);
xB  = reshape(B,220,8);
xF  = reshape(F,220,8);
xP  = reshape(P,220,8);
xD  = reshape(D,220,8);
xS  = reshape(S,220,8);
xM  = reshape(M,220,8);
xL  = reshape(L,220,8);
xall = reshape(All,220,8);

msf = NaN*ones(11,8);
msp = NaN*ones(11,8);
msd = NaN*ones(11,8);
mmf = NaN*ones(11,8);
mmp = NaN*ones(11,8);
mmd = NaN*ones(11,8);
mlp = NaN*ones(11,8);
mld = NaN*ones(11,8);
mB  = NaN*ones(11,8);
mF  = NaN*ones(11,8);
mP  = NaN*ones(11,8);
mD  = NaN*ones(11,8);
mS  = NaN*ones(11,8);
mM  = NaN*ones(11,8);
mL  = NaN*ones(11,8);
mall = NaN*ones(11,8);

ssf = NaN*ones(11,8);
ssp = NaN*ones(11,8);
ssd = NaN*ones(11,8);
smf = NaN*ones(11,8);
smp = NaN*ones(11,8);
smd = NaN*ones(11,8);
slp = NaN*ones(11,8);
sld = NaN*ones(11,8);
sB  = NaN*ones(11,8);
sF  = NaN*ones(11,8);
sP  = NaN*ones(11,8);
sD  = NaN*ones(11,8);
sS  = NaN*ones(11,8);
sM  = NaN*ones(11,8);
sL  = NaN*ones(11,8);
sall = NaN*ones(11,8);

%%
st=1:20:220;
en=20:20:220;

for n=1:11
    msf(n,:) = mean(xsf(st(n):en(n),:),1);
    msp(n,:) = mean(xsp(st(n):en(n),:),1);
    msd(n,:) = mean(xsd(st(n):en(n),:),1);
    mmf(n,:) = mean(xmf(st(n):en(n),:),1);
    mmp(n,:) = mean(xmp(st(n):en(n),:),1);
    mmd(n,:) = mean(xmd(st(n):en(n),:),1);
    mlp(n,:) = mean(xlp(st(n):en(n),:),1);
    mld(n,:) = mean(xld(st(n):en(n),:),1);
    mB(n,:)  = mean(xB(st(n):en(n),:),1);
    mF(n,:)  = mean(xF(st(n):en(n),:),1);
    mP(n,:)  = mean(xP(st(n):en(n),:),1);
    mD(n,:)  = mean(xD(st(n):en(n),:),1);
    mS(n,:)  = mean(xS(st(n):en(n),:),1);
    mM(n,:)  = mean(xM(st(n):en(n),:),1);
    mL(n,:)  = mean(xL(st(n):en(n),:),1);
    mall(n,:) = mean(xall(st(n):en(n),:),1);
    
    ssf(n,:) = std(xsf(st(n):en(n),:),1);
    ssp(n,:) = std(xsp(st(n):en(n),:),1);
    ssd(n,:) = std(xsd(st(n):en(n),:),1);
    smf(n,:) = std(xmf(st(n):en(n),:),1);
    smp(n,:) = std(xmp(st(n):en(n),:),1);
    smd(n,:) = std(xmd(st(n):en(n),:),1);
    slp(n,:) = std(xlp(st(n):en(n),:),1);
    sld(n,:) = std(xld(st(n):en(n),:),1);
    sB(n,:)  = std(xB(st(n):en(n),:),1);
    sF(n,:)  = std(xF(st(n):en(n),:),1);
    sP(n,:)  = std(xP(st(n):en(n),:),1);
    sD(n,:)  = std(xD(st(n):en(n),:),1);
    sS(n,:)  = std(xS(st(n):en(n),:),1);
    sM(n,:)  = std(xM(st(n):en(n),:),1);
    sL(n,:)  = std(xL(st(n):en(n),:),1);
    sall(n,:) = std(xall(st(n):en(n),:),1);
end

%% Save
save([fpath 'powerspec_feisty_means_stds_ln.mat'],...
    'msf','msp','msd','mmf','mmp','mmd','mlp','mld',...
    'mB','mF','mP','mD','mS','mM','mL','mall',...
    'ssf','ssp','ssd','smf','smp','smd','slp','sld',...
    'sB','sF','sP','sD','sS','sM','sL','sall')

%% Plot mean and std vs. slope for each biome and type
biome = {'SLC','SECCS','SECSS','SCoast','NLC','NECCS','NECSS','NCoast'};
alp = [0:-0.15:-1.5];

cm9=[0.5 0.5 0.5 ;...   %grey
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0 0.7 0]; %...    %green
    %0 0 0];...      %black
    
set(groot,'defaultAxesColorOrder',cm9);

cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

%%
figure(1)
subplot(3,3,1)
plot(alp,msf,'LineWidth',2); hold on;
% plot(alp,msf+ssf); hold on;
% plot(alp,msf-ssf); hold on;
plot(alp,alp,'--k')
title('SF')
ylim([-5 0.5])

subplot(3,3,2)
plot(alp,msp,'LineWidth',2); hold on;
% plot(alp,msp+ssp); hold on;
% plot(alp,msp-ssp); hold on;
plot(alp,alp,'--k')
title('SP')
ylim([-5 0.5])

subplot(3,3,3)
plot(alp,msd,'LineWidth',2); hold on;
% plot(alp,msd+ssd); hold on;
% plot(alp,msd-ssd); hold on;
plot(alp,alp,'--k')
title('SD')
ylim([-5 0.5])

subplot(3,3,4)
plot(alp,mmf,'LineWidth',2); hold on;
% plot(alp,mmf+smf); hold on;
% plot(alp,mmf-smf); hold on;
plot(alp,alp,'--k')
title('MF')
ylim([-5 0.5])

subplot(3,3,5)
plot(alp,mmp,'LineWidth',2); hold on;
% plot(alp,mmp+smp); hold on;
% plot(alp,mmp-smp); hold on;
plot(alp,alp,'--k')
title('MP')
ylim([-5 0.5])

subplot(3,3,6)
plot(alp,mmd,'LineWidth',2); hold on;
% plot(alp,mmd+smd); hold on;
% plot(alp,mmd-smd); hold on;
plot(alp,alp,'--k')
title('MD')
ylim([-6 0.5])

subplot(3,3,7)
plot(alp,mmf,'LineWidth',2); hold on;
ylim([-1 0])
xlim([-1.5 0])
legend(biome)
set(gca,'YTickLabels',[])

subplot(3,3,8)
plot(alp,mlp,'LineWidth',2); hold on;
% plot(alp,mlp+slp); hold on;
% plot(alp,mlp-slp); hold on;
plot(alp,alp,'--k')
title('LP')
ylim([-5 0.5])

subplot(3,3,9)
plot(alp,mld,'LineWidth',2); hold on;
% plot(alp,mld+sld); hold on;
% plot(alp,mld-sld); hold on;
plot(alp,alp,'--k')
title('LD')
ylim([-5 0.5])
print('-dpng',[ppath exper 'ps_slopes_all_stages_ln.png'])

%%
figure(2)
subplot(3,3,1)
plot(alp,mF,'LineWidth',2); hold on;
% plot(alp,msf+ssf); hold on;
% plot(alp,msf-ssf); hold on;
plot(alp,alp,'--k')
title('F')
ylim([-5 0.5])

subplot(3,3,2)
plot(alp,mP,'LineWidth',2); hold on;
% plot(alp,msp+ssp); hold on;
% plot(alp,msp-ssp); hold on;
plot(alp,alp,'--k')
title('P')
ylim([-5 0.5])

subplot(3,3,3)
plot(alp,mD,'LineWidth',2); hold on;
% plot(alp,msd+ssd); hold on;
% plot(alp,msd-ssd); hold on;
plot(alp,alp,'--k')
title('D')
ylim([-5 0.5])

subplot(3,3,4)
plot(alp,mS,'LineWidth',2); hold on;
% plot(alp,mmf+smf); hold on;
% plot(alp,mmf-smf); hold on;
plot(alp,alp,'--k')
title('S')
ylim([-5 0.5])

subplot(3,3,5)
plot(alp,mM,'LineWidth',2); hold on;
% plot(alp,mmp+smp); hold on;
% plot(alp,mmp-smp); hold on;
plot(alp,alp,'--k')
title('M')
ylim([-5 0.5])

subplot(3,3,6)
plot(alp,mL,'LineWidth',2); hold on;
% plot(alp,mmd+smd); hold on;
% plot(alp,mmd-smd); hold on;
plot(alp,alp,'--k')
title('L')
ylim([-5 0.5])

subplot(3,3,7)
plot(alp,mall,'LineWidth',2); hold on;
ylim([-1 0])
xlim([-1.5 0])
legend(biome)
set(gca,'YTickLabels',[])

subplot(3,3,8)
plot(alp,mall,'LineWidth',2); hold on;
% plot(alp,mlp+slp); hold on;
% plot(alp,mlp-slp); hold on;
plot(alp,alp,'--k')
title('All')
ylim([-5 0.5])

subplot(3,3,9)
plot(alp,mB,'LineWidth',2); hold on;
% plot(alp,mld+sld); hold on;
% plot(alp,mld-sld); hold on;
plot(alp,alp,'--k')
title('B')
ylim([-5 0.5])
print('-dpng',[ppath exper 'ps_slopes_all_groups_ln.png'])

%%
figure(3)
for i=1:8
    subplot(3,3,i)
    plot(alp,msf(:,i),'r','LineWidth',2); hold on;
    plot(alp,msf(:,i)+ssf(:,i),'r'); hold on;
    plot(alp,msf(:,i)-ssf(:,i),'r'); hold on;
    
    plot(alp,mmf(:,i),'color',cm9(4,:),'LineWidth',2); hold on;
    plot(alp,mmf(:,i)+smf(:,i),'color',cm9(4,:)); hold on;
    plot(alp,mmf(:,i)-smf(:,i),'color',cm9(4,:)); hold on;
    plot(alp,alp,'--k')
    
    if (i==2)
        title({'F', biome{i}})
    else
    title(biome(i))
    end
end
subplot(3,3,9)
plot(alp,msf(:,i),'r','LineWidth',2); hold on;
plot(alp,mmf(:,i),'color',cm9(4,:),'LineWidth',2); hold on;
legend('SF','MF')
ylim([0 1])
set(gca,'YTickLabel',[])
print('-dpng',[ppath exper 'ps_slopes_allF_ln.png'])

%%
figure(4)
for i=1:8
    subplot(3,3,i)
    plot(alp,msp(:,i),'color',cm9(5,:),'LineWidth',2); hold on;
    plot(alp,msp(:,i)+ssp(:,i),'color',cm9(5,:)); hold on;
    plot(alp,msp(:,i)-ssp(:,i),'color',cm9(5,:)); hold on;
    
    plot(alp,mmp(:,i),'color',cm9(6,:),'LineWidth',2); hold on;
    plot(alp,mmp(:,i)+smp(:,i),'color',cm9(6,:)); hold on;
    plot(alp,mmp(:,i)-smp(:,i),'color',cm9(6,:)); hold on;
    
    plot(alp,mlp(:,i),'color',cm9(7,:),'LineWidth',2); hold on;
    plot(alp,mlp(:,i)+slp(:,i),'color',cm9(7,:)); hold on;
    plot(alp,mlp(:,i)-slp(:,i),'color',cm9(7,:)); hold on;
    plot(alp,alp,'--k')
    if (i==2)
        title({'P', biome{i}})
    else
    title(biome(i))
    end
end
subplot(3,3,9)
plot(alp,msp(:,i),'color',cm9(5,:),'LineWidth',2); hold on;
plot(alp,mmp(:,i),'color',cm9(6,:),'LineWidth',2); hold on;
plot(alp,mlp(:,i),'color',cm9(7,:),'LineWidth',2);
legend('SP','MP','LP')
ylim([0 1])
set(gca,'YTickLabel',[])
print('-dpng',[ppath exper 'ps_slopes_allP_ln.png'])

%%
figure(5)
for i=1:8
    subplot(3,3,i)
    plot(alp,msd(:,i),'color',cm21(15,:),'LineWidth',2); hold on;
    plot(alp,msd(:,i)+ssd(:,i),'color',cm21(15,:)); hold on;
    plot(alp,msd(:,i)-ssd(:,i),'color',cm21(15,:)); hold on;
    
    plot(alp,mmd(:,i),'color',cm9(8,:),'LineWidth',2); hold on;
    plot(alp,mmd(:,i)+smd(:,i),'color',cm9(8,:)); hold on;
    plot(alp,mmd(:,i)-smd(:,i),'color',cm9(8,:)); hold on;
    
    plot(alp,mld(:,i),'color',cm21(16,:),'LineWidth',2); hold on;
    plot(alp,mld(:,i)+sld(:,i),'color',cm21(16,:)); hold on;
    plot(alp,mld(:,i)-sld(:,i),'color',cm21(16,:)); hold on;
    plot(alp,alp,'--k')
    %ylim([-4 0])
    if (i==2)
        title({'D', biome{i}})
    else
    title(biome(i))
    end
end
subplot(3,3,9)
plot(alp,msd(:,i),'color',cm21(15,:),'LineWidth',2); hold on;
plot(alp,mmd(:,i),'color',cm9(8,:),'LineWidth',2); hold on;
plot(alp,mld(:,i),'color',cm21(16,:),'LineWidth',2);
legend('SD','MD','LD')
ylim([0 1])
set(gca,'YTickLabel',[])
print('-dpng',[ppath exper 'ps_slopes_allD_ln.png'])

%%
figure(6)
for i=1:8
    subplot(3,3,i)
    plot(alp,mB(:,i),'k','LineWidth',2); hold on;
    plot(alp,mB(:,i)+sB(:,i),'k'); hold on;
    plot(alp,mB(:,i)-sB(:,i),'k'); hold on;
    plot(alp,alp,'--b')
    ylim([-4 0])
    if (i==2)
        title({'B', biome{i}})
    else
    title(biome(i))
    end
end
print('-dpng',[ppath exper 'ps_slopes_B_ln.png'])

%%
figure(7)
for i=1:8
    subplot(3,3,i)
    plot(alp,mF(:,i),'r','LineWidth',2); hold on;
    plot(alp,mF(:,i)+sF(:,i),'r'); hold on;
    plot(alp,mF(:,i)-sF(:,i),'r'); hold on;
    
    plot(alp,mP(:,i),'b','LineWidth',2); hold on;
    plot(alp,mP(:,i)+sP(:,i),'b'); hold on;
    plot(alp,mP(:,i)-sP(:,i),'b'); hold on;
    
    plot(alp,mD(:,i),'color',cm9(8,:),'LineWidth',2); hold on;
    plot(alp,mD(:,i)+sD(:,i),'color',cm9(8,:)); hold on;
    plot(alp,mD(:,i)-sD(:,i),'color',cm9(8,:)); hold on;
    plot(alp,alp,'--k')
    ylim([-4.5 0])
    if (i==2)
        title({'F,P,D', biome{i}})
    else
    title(biome(i))
    end
end
subplot(3,3,9)
plot(alp,mF(:,i),'r','LineWidth',2); hold on;
plot(alp,mP(:,i),'b','LineWidth',2); hold on;
plot(alp,mD(:,i),'color',cm9(8,:),'LineWidth',2);
legend('F','P','D')
ylim([0 1])
set(gca,'YTickLabel',[])
print('-dpng',[ppath exper 'ps_slopes_all_types_ln.png'])

%%
figure(8)
for i=1:8
    subplot(3,3,i)
    plot(alp,mS(:,i),'color',cm21(20,:),'LineWidth',2); hold on;
    plot(alp,mS(:,i)+sS(:,i),'color',cm21(20,:)); hold on;
    plot(alp,mS(:,i)-sS(:,i),'color',cm21(20,:)); hold on;
    
    plot(alp,mM(:,i),'m','LineWidth',2); hold on;
    plot(alp,mM(:,i)+sM(:,i),'m'); hold on;
    plot(alp,mM(:,i)-sM(:,i),'m'); hold on;
    
    plot(alp,mL(:,i),'color',cm21(6,:),'LineWidth',2); hold on;
    plot(alp,mL(:,i)+sL(:,i),'color',cm21(6,:)); hold on;
    plot(alp,mL(:,i)-sL(:,i),'color',cm21(6,:)); hold on;
    plot(alp,alp,'--k')
    ylim([-5 0])
    if (i==2)
        title({'S,M,L', biome{i}})
    else
    title(biome(i))
    end
end
subplot(3,3,9)
plot(alp,mS(:,i),'color',cm21(20,:),'LineWidth',2); hold on;
plot(alp,mM(:,i),'m','LineWidth',2); hold on;
plot(alp,mL(:,i),'color',cm21(6,:),'LineWidth',2);
legend('S','M','L')
ylim([0 1])
set(gca,'YTickLabel',[])
print('-dpng',[ppath exper 'ps_slopes_all_sizes_ln.png'])

%%
figure(9)
for i=1:8
    subplot(3,3,i)
    plot(alp,mall(:,i),'k','LineWidth',2); hold on;
    plot(alp,mall(:,i)+sall(:,i),'k'); hold on;
    plot(alp,mall(:,i)-sall(:,i),'k'); hold on;
    plot(alp,alp,'--b')
    ylim([-5 0])
    if (i==2)
        title({'All', biome{i}})
    else
    title(biome(i))
    end
end
print('-dpng',[ppath exper 'ps_slopes_All_fish_ln.png'])

