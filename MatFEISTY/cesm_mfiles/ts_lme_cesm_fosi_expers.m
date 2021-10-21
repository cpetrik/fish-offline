% Compare FEISTY runs w/ CESM FOSI experiments
% Time series plots and maps

clear all
close all

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI//';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_noCC_RE00100';
mod = 'v13_All_fish03';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% Full
load([fpath 'LME_fosi_fished_v13_All_fish03_',cfile '.mat']);

FullmF = lme_msfb + lme_mmfb;
FullmP = lme_mspb + lme_mmpb + lme_mlpb;
FullmD = lme_msdb + lme_mmdb + lme_mldb;
FullmB = lme_mbb;
FullmS = lme_msfb + lme_mspb + lme_msdb;
FullmM = lme_mmfb + lme_mmpb + lme_mmdb;
FullmL = lme_mlpb + lme_mldb;
FullmA = FullmF+FullmP+FullmD;

FulltF = lme_ssfb + lme_smfb;
FulltP = lme_sspb + lme_smpb + lme_slpb;
FulltD = lme_ssdb + lme_smdb + lme_sldb;
FulltB = lme_sbb;
FulltS = lme_ssfb + lme_sspb + lme_ssdb;
FulltM = lme_smfb + lme_smpb + lme_smdb;
FulltL = lme_slpb + lme_sldb;
FulltA = FulltF+FulltP+FulltD;

clear lme_msfb lme_mspb lme_msdb lme_mmfb lme_mmpb lme_mmdb
clear lme_mlpb lme_mldb lme_mbb
clear lme_ssfb lme_sspb lme_ssdb lme_smfb lme_smpb lme_smdb
clear lme_slpb lme_sldb lme_sbb

%%
%exper = varFood, varTemp, climatol
load([fpath 'LME_fosi_fished_v13_climatol_',cfile '.mat']);

ClimmF = lme_msfb + lme_mmfb;
ClimmP = lme_mspb + lme_mmpb + lme_mlpb;
ClimmD = lme_msdb + lme_mmdb + lme_mldb;
ClimmB = lme_mbb;
ClimmS = lme_msfb + lme_mspb + lme_msdb;
ClimmM = lme_mmfb + lme_mmpb + lme_mmdb;
ClimmL = lme_mlpb + lme_mldb;
ClimmA = ClimmF+ClimmP+ClimmD;

ClimtF = lme_ssfb + lme_smfb;
ClimtP = lme_sspb + lme_smpb + lme_slpb;
ClimtD = lme_ssdb + lme_smdb + lme_sldb;
ClimtB = lme_sbb;
ClimtS = lme_ssfb + lme_sspb + lme_ssdb;
ClimtM = lme_smfb + lme_smpb + lme_smdb;
ClimtL = lme_slpb + lme_sldb;
ClimtA = ClimtF+ClimtP+ClimtD;

clear lme_msfb lme_mspb lme_msdb lme_mmfb lme_mmpb lme_mmdb
clear lme_mlpb lme_mldb lme_mbb
clear lme_ssfb lme_sspb lme_ssdb lme_smfb lme_smpb lme_smdb
clear lme_slpb lme_sldb lme_sbb

%%
load([fpath 'LME_fosi_fished_v13_varTemp_',cfile '.mat']);

TempmF = lme_msfb + lme_mmfb;
TempmP = lme_mspb + lme_mmpb + lme_mlpb;
TempmD = lme_msdb + lme_mmdb + lme_mldb;
TempmB = lme_mbb;
TempmS = lme_msfb + lme_mspb + lme_msdb;
TempmM = lme_mmfb + lme_mmpb + lme_mmdb;
TempmL = lme_mlpb + lme_mldb;
TempmA = TempmF+TempmP+TempmD;

TemptF = lme_ssfb + lme_smfb;
TemptP = lme_sspb + lme_smpb + lme_slpb;
TemptD = lme_ssdb + lme_smdb + lme_sldb;
TemptB = lme_sbb;
TemptS = lme_ssfb + lme_sspb + lme_ssdb;
TemptM = lme_smfb + lme_smpb + lme_smdb;
TemptL = lme_slpb + lme_sldb;
TemptA = TemptF+TemptP+TemptD;

clear lme_msfb lme_mspb lme_msdb lme_mmfb lme_mmpb lme_mmdb
clear lme_mlpb lme_mldb lme_mbb
clear lme_ssfb lme_sspb lme_ssdb lme_smfb lme_smpb lme_smdb
clear lme_slpb lme_sldb lme_sbb

%%
load([fpath 'LME_fosi_fished_v13_varFood_',cfile '.mat']);

FoodmF = lme_msfb + lme_mmfb;
FoodmP = lme_mspb + lme_mmpb + lme_mlpb;
FoodmD = lme_msdb + lme_mmdb + lme_mldb;
FoodmB = lme_mbb;
FoodmS = lme_msfb + lme_mspb + lme_msdb;
FoodmM = lme_mmfb + lme_mmpb + lme_mmdb;
FoodmL = lme_mlpb + lme_mldb;
FoodmA = FoodmF+FoodmP+FoodmD;

FoodtF = lme_ssfb + lme_smfb;
FoodtP = lme_sspb + lme_smpb + lme_slpb;
FoodtD = lme_ssdb + lme_smdb + lme_sldb;
FoodtB = lme_sbb;
FoodtS = lme_ssfb + lme_sspb + lme_ssdb;
FoodtM = lme_smfb + lme_smpb + lme_smdb;
FoodtL = lme_slpb + lme_sldb;
FoodtA = FoodtF+FoodtP+FoodtD;

clear lme_msfb lme_mspb lme_msdb lme_mmfb lme_mmpb lme_mmdb
clear lme_mlpb lme_mldb lme_mbb
clear lme_ssfb lme_sspb lme_ssdb lme_smfb lme_smpb lme_smdb
clear lme_slpb lme_sldb lme_sbb

%% Diffs from climatol
dFoodtF = FoodtF - ClimtF;
dFoodtP = FoodtP - ClimtP;
dFoodtD = FoodtD - ClimtD;
dFoodtB = FoodtB - ClimtB;
dFoodtS = FoodtS - ClimtS;
dFoodtM = FoodtM - ClimtM;
dFoodtL = FoodtL - ClimtL;
dFoodtA = FoodtA - ClimtA;

dFoodmF = FoodmF - ClimmF;
dFoodmP = FoodmP - ClimmP;
dFoodmD = FoodmD - ClimmD;
dFoodmB = FoodmB - ClimmB;
dFoodmS = FoodmS - ClimmS;
dFoodmM = FoodmM - ClimmM;
dFoodmL = FoodmL - ClimmL;
dFoodmA = FoodmA - ClimmA;

dTemptF = TemptF - ClimtF;
dTemptP = TemptP - ClimtP;
dTemptD = TemptD - ClimtD;
dTemptB = TemptB - ClimtB;
dTemptS = TemptS - ClimtS;
dTemptM = TemptM - ClimtM;
dTemptL = TemptL - ClimtL;
dTemptA = TemptA - ClimtA;

dTempmF = TempmF - ClimmF;
dTempmP = TempmP - ClimmP;
dTempmD = TempmD - ClimmD;
dTempmB = TempmB - ClimmB;
dTempmS = TempmS - ClimmS;
dTempmM = TempmM - ClimmM;
dTempmL = TempmL - ClimmL;
dTempmA = TempmA - ClimmA;

dFulltF = FulltF - ClimtF;
dFulltP = FulltP - ClimtP;
dFulltD = FulltD - ClimtD;
dFulltB = FulltB - ClimtB;
dFulltS = FulltS - ClimtS;
dFulltM = FulltM - ClimtM;
dFulltL = FulltL - ClimtL;
dFulltA = FulltA - ClimtA;

dFullmF = FullmF - ClimmF;
dFullmP = FullmP - ClimmP;
dFullmD = FullmD - ClimmD;
dFullmB = FullmB - ClimmB;
dFullmS = FullmS - ClimmS;
dFullmM = FullmM - ClimmM;
dFullmL = FullmL - ClimmL;
dFullmA = FullmA - ClimmA;

%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black
    
set(groot,'defaultAxesColorOrder',cm10);

%% US LMEs only
lid = [54,1:2,10,3,5:7];
lname = {'CHK','EBS','GAK','HI','CCE','GMX','SE','NE'};

% Plots in time
t = 1:length(ClimtF); %time;
y = t/12;

%% Benthos
figure(1)
for z = 1:length(lid)
    lme = lid(z);
    subplot(3,3,z)
    plot(y,log10(FullmB(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmB(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmB(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmB(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    %ylim([-5 2])
    if z==8
        xlabel('Time (y)')
    end
    if z==4
        ylabel('log_1_0 Biomass (g m^-^2)')
    end
    if z==2
        title(['Benthos ' lname(z)])
    else
        title(lname(z))
    end
end
lg  = legend({'full','clim','temp','prey'}); 
lg.Position(1:2) = [.69 .21];
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_USLME_exper_ts_B.png'])

%% Forage
figure(2)
for z = 1:length(lid)
    lme = lid(z);
    subplot(3,3,z)
    plot(y,log10(FullmF(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmF(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmF(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmF(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    %ylim([-5 2])
    if z==8
        xlabel('Time (y)')
    end
    if z==4
        ylabel('log_1_0 Biomass (g m^-^2)')
    end
    if z==2
        title(['Forage ' lname(z)])
    else
        title(lname(z))
    end
end
lg  = legend({'full','clim','temp','prey'}); 
lg.Position(1:2) = [.69 .21];
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_USLME_exper_ts_F.png'])

%% Lg Pel
figure(3)
for z = 1:length(lid)
    lme = lid(z);
    subplot(3,3,z)
    plot(y,log10(FullmP(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmP(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmP(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmP(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    %ylim([-5 2])
    if z==8
        xlabel('Time (y)')
    end
    if z==4
        ylabel('log_1_0 Biomass (g m^-^2)')
    end
    if z==2
        title(['Lg Pelagic ' lname(z)])
    else
        title(lname(z))
    end
end
lg  = legend({'full','clim','temp','prey'}); 
lg.Position(1:2) = [.69 .21];
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_USLME_exper_ts_P.png'])

%% Dem
figure(4)
for z = 1:length(lid)
    lme = lid(z);
    subplot(3,3,z)
    plot(y,log10(FullmD(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmD(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmD(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmD(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    %ylim([-5 2])
    if z==8
        xlabel('Time (y)')
    end
    if z==4
        ylabel('log_1_0 Biomass (g m^-^2)')
    end
    if z==2
        title(['Demersal ' lname(z)])
    else
        title(lname(z))
    end
end
lg  = legend({'full','clim','temp','prey'}); 
lg.Position(1:2) = [.69 .21];
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_USLME_exper_ts_D.png'])

%% All
figure(5)
for z = 1:length(lid)
    lme = lid(z);
    subplot(3,3,z)
    plot(y,log10(FullmA(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmA(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmA(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmA(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    %ylim([-5 2])
    if z==8
        xlabel('Time (y)')
    end
    if z==4
        ylabel('log_1_0 Biomass (g m^-^2)')
    end
    if z==2
        title(['All ' lname(z)])
    else
        title(lname(z))
    end
end
lg  = legend({'full','clim','temp','prey'}); 
lg.Position(1:2) = [.69 .21];
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_USLME_exper_ts_All.png'])

%% Small
figure(6)
for z = 1:length(lid)
    lme = lid(z);
    subplot(3,3,z)
    plot(y,log10(FullmS(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmS(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmS(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmS(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    %ylim([-5 2])
    if z==8
        xlabel('Time (y)')
    end
    if z==4
        ylabel('log_1_0 Biomass (g m^-^2)')
    end
    if z==2
        title(['Small ' lname(z)])
    else
        title(lname(z))
    end
end
lg  = legend({'full','clim','temp','prey'}); 
lg.Position(1:2) = [.69 .21];
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_USLME_exper_ts_S.png'])

%% Med
figure(7)
for z = 1:length(lid)
    lme = lid(z);
    subplot(3,3,z)
    plot(y,log10(FullmM(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmM(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmM(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmM(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    %ylim([-5 2])
    if z==8
        xlabel('Time (y)')
    end
    if z==4
        ylabel('log_1_0 Biomass (g m^-^2)')
    end
    if z==2
        title(['Medium ' lname(z)])
    else
        title(lname(z))
    end
end
lg  = legend({'full','clim','temp','prey'}); 
lg.Position(1:2) = [.69 .21];
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_USLME_exper_ts_M.png'])

%% Lrg
figure(8)
for z = 1:length(lid)
    lme = lid(z);
    subplot(3,3,z)
    plot(y,log10(FullmL(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmL(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmL(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmL(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    %ylim([-5 2])
    if z==8
        xlabel('Time (y)')
    end
    if z==4
        ylabel('log_1_0 Biomass (g m^-^2)')
    end
    if z==2
        title(['Large ' lname(z)])
    else
        title(lname(z))
    end
end
lg  = legend({'full','clim','temp','prey'}); 
lg.Position(1:2) = [.69 .21];
stamp('')
print('-dpng',[ppath 'FOSI_',mod,'_USLME_exper_ts_L.png'])

%% Subplots of each fn type per LME
for z = 1:length(lid)
    lme = lid(z);
    f9 = figure(9);
    clf
    subplot(3,3,1)
    plot(y,log10(FullmS(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmS(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmS(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmS(lme,:)),'-.b','Linewidth',2); hold on;
    title('Small')
    
    subplot(3,3,2)
    plot(y,log10(FullmM(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmM(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmM(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmM(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title([lname(z) ' Medium'])
    
    subplot(3,3,3)
    plot(y,log10(FullmL(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmL(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmL(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmL(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    %ylim([-5 2])
    title('Large')
    
    subplot(3,3,4)
    plot(y,log10(FullmF(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmF(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmF(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmF(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    ylabel('log_1_0 Biomass (g m^-^2)')
    title('Forage')
    
    subplot(3,3,5)
    plot(y,log10(FullmP(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmP(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmP(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmP(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Lg Pelagic')
    
    subplot(3,3,6)
    plot(y,log10(FullmD(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmD(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmD(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmD(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Demersal')
    
    subplot(3,3,7)
    plot(y,log10(FullmA(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmA(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmA(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmA(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    xlabel('Time (y)')
    title('All Fish')
    
    subplot(3,3,9)
    plot(y,log10(FullmB(lme,:)),'k','Linewidth',2); hold on;
    plot(y,log10(ClimmB(lme,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(TempmB(lme,:)),'-.r','Linewidth',2); hold on;
    plot(y,log10(FoodmB(lme,:)),'-.b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Benthos')
    lg  = legend({'full','clim','temp','prey'});
    lg.Position(1:2) = [.45 .21];
    stamp('')
    print(f9,'-dpng',[ppath 'FOSI_',mod,'_USLME_exper_ts_',lname{z},'.png'])
end




