% Visualize output of FEISTY
% CESM DPLE
% Time series plots and maps

clear
close all

%% Map data
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

[ni,nj]=size(TLONG);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/DPLE/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/DPLE/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%pick year
StartYr = 1954;
%loop over members
submem = 1:40;

%%
Ftmean = nan*ones(length(submem),122);
Ptmean = Ftmean;
Dtmean = Ftmean;
Btmean = Ftmean;
Stmean = Ftmean;
Mtmean = Ftmean;
Ltmean = Ftmean;

Fsmean = nan*ones(ni,nj,length(submem));
Psmean = Fsmean;
Dsmean = Fsmean;
Bsmean = Fsmean;
Ssmean = Fsmean;
Msmean = Fsmean;
Lsmean = Fsmean;

for mem=1:length(submem) %will loop over

    Member = submem(mem);
    exper = ['v15_Y' num2str(StartYr) '_M' num2str(Member) '_All_fish03_' ];

    load([fpath 'Time_Means_DPLE_' exper cfile '.mat']);
    load([fpath 'Space_Means_DPLE_' exper cfile '.mat']);

    %% Types together
    F = sf_tmean+mf_tmean;
    P = sp_tmean+mp_tmean+lp_tmean;
    D = sd_tmean+md_tmean+ld_tmean;
    B = b_tmean;
    S = sf_tmean+sp_tmean+sd_tmean;
    M = mf_tmean+mp_tmean+md_tmean;
    L = lp_tmean+ld_tmean;

    %% Plots in space

    Zsf=NaN*ones(ni,nj);
    Zsp=NaN*ones(ni,nj);
    Zsd=NaN*ones(ni,nj);
    Zmf=NaN*ones(ni,nj);
    Zmp=NaN*ones(ni,nj);
    Zmd=NaN*ones(ni,nj);
    Zlp=NaN*ones(ni,nj);
    Zld=NaN*ones(ni,nj);
    Zb=NaN*ones(ni,nj);

    Zsf(GRD.ID)=sf_sbio;
    Zsp(GRD.ID)=sp_sbio;
    Zsd(GRD.ID)=sd_sbio;
    Zmf(GRD.ID)=mf_sbio;
    Zmp(GRD.ID)=mp_sbio;
    Zmd(GRD.ID)=md_sbio;
    Zlp(GRD.ID)=lp_sbio;
    Zld(GRD.ID)=ld_sbio;
    Zb(GRD.ID)=b_sbio;

    All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
    AllF = Zsf+Zmf;
    AllP = Zsp+Zmp+Zlp;
    AllD = Zsd+Zmd+Zld;
    AllS = Zsp+Zsf+Zsd;
    AllM = Zmp+Zmf+Zmd;
    AllL = Zlp+Zld;
    FracPD = AllP ./ (AllP+AllD);
    FracPF = AllP ./ (AllP+AllF);
    FracLM = AllL ./ (AllL+AllM);

    %% save outputs for comparison
    save([fpath 'Time_Means_DPLE_' exper cfile '.mat'],'F','P','D','B',...
        'S','M','L','-append');

    save([fpath 'Space_Means_DPLE_' exper cfile '.mat'],...
        'AllF','AllP','AllD','AllS','AllM','AllL','-append');

    Ftmean(mem,:) = F;
    Ptmean(mem,:) = P;
    Dtmean(mem,:) = D;
    Btmean(mem,:) = B;
    Stmean(mem,:) = S;
    Mtmean(mem,:) = M;
    Ltmean(mem,:) = L;

    Fsmean(:,:,mem) = AllF;
    Psmean(:,:,mem) = AllP;
    Dsmean(:,:,mem) = AllD;
    Bsmean(:,:,mem) = Zb;
    Ssmean(:,:,mem) = AllS;
    Msmean(:,:,mem) = AllM;
    Lsmean(:,:,mem) = AllL;


end
exper2 = ['v15_Y' num2str(StartYr) '_All_fish03_' ];
save([fpath 'Plot_Means_DPLE_' exper2 cfile '.mat'],...
    'Ftmean','Ptmean','Dtmean','Btmean','Stmean','Mtmean','Ltmean',...
    'Fsmean','Psmean','Dsmean','Bsmean','Ssmean','Msmean','Lsmean');

%%
load([fpath 'Plot_Means_DPLE_' exper2 cfile '.mat'])

%% Plots in time
t = 1:122; %time;
y = (t/12) + (10/12) + StartYr;

% All types
figure(1)
subplot(2,2,1)
plot(y,log10(Ftmean),'r','Linewidth',1); hold on;
plot(y,log10(mean(Ftmean)),'color',[0.75 0 0],'Linewidth',3); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('F')

subplot(2,2,2)
plot(y,log10(Ptmean),'b','Linewidth',1); hold on;
plot(y,log10(mean(Ptmean)),'color',[0 0 0.75],'Linewidth',3); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('P')

subplot(2,2,3)
plot(y,log10(Dtmean),'color',[0 0.75 0.5],'Linewidth',1); hold on;
plot(y,log10(mean(Dtmean)),'color',[0 0.5 0.25],'Linewidth',3); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('D')

subplot(2,2,4)
plot(y,log10(Btmean),'color',[0.5 0.5 0.5],'Linewidth',3); hold on;
plot(y,log10(mean(Btmean)),'k','Linewidth',2); hold on;
xlim([y(1) y(end)])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('B')
stamp(exper2)
print('-dpng',[ppath 'DPLE_',exper2,'all_types.png'])


