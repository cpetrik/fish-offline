% Visualize output of FEISTY testcase

clear all
close all

%%
datap = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/';
figp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/testcase/';

dp = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ100_mMZ045_nmort1_BE08_noCC_RE00100';
sname = 'testcase_locs3_All_fish03';
harv = 'All_fish03';
dpath = [datap char(dp) '/'];
fpath = [figp char(dp) '/'];
if (~isfolder(fpath))
    mkdir(fpath)
end
cfile = char(dp);

%%
load([dpath sname '.mat'])
[nt,nx,nid] = size(biomass);
stages = {'SF','SP','SD','MF','MP','MD','LP','LD','B'};
locs = {'EBS', 'HOT', 'PUP'}';

%%
load('cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);

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

%%
all_mean=NaN*ones(3,4,nx);
z = NaN*ones(nx,3);

t=1:nt;

%% Biomass
F = biomass(:,:,1) + biomass(:,:,4);
P = biomass(:,:,2) + biomass(:,:,5) + biomass(:,:,7);
D = biomass(:,:,3) + biomass(:,:,6) + biomass(:,:,8);
B = biomass(:,:,9);

% Growth rate = nu - energy for biomass production)
% Consump per biomass (I) by type
% Consump per biomass (I)
% Feeding level = con/cmax
% Production (= nu * biom)
% Gross growth efficiency (= nu/consump)
gge = energy_avail_rate ./ ingestion_rate;
% Reproduction
% Recruitment
% Metabolism
% Predation (per biomass)
% Natural mortality (nmort/biom)
% Fishing mort = fish
% Fishing rate = frate
% Fishing catch = totcatch
% Total mortality w/o fishing = pred + nat
tot1 = mortality_rate + predation_rate;
% Total mortality w/ fishing (all rates per biomass) = pred + nat + frate
tot2 = mortality_rate + predation_rate + fish_catch_rate;
    
%% Loop over locations

for s=1:nx
    
    %% TIME SERIES ----------------------------------------------------
    loc = locs{s};
    lname = [loc];
    y=t;
    
    %% Biomass
    figure(1)
    clf
    plot(y,log10(squeeze(biomass(:,s,:))),'Linewidth',1); hold on;
    legend(stages)
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-5.5 2])
    xlabel('Day')
    ylabel('log10 Biomass (g m^-^2)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_logmean_biomass_' lname '.png'])
    
    figure(2)
    clf
    plot(y,log10(B(:,s,:)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
    plot(y,log10(F(:,s,:)),'r','Linewidth',2); hold on;
    plot(y,log10(P(:,s,:)),'b','Linewidth',2); hold on;
    plot(y,log10(D(:,s,:)),'k','Linewidth',2); hold on;
    legend('B','F','P','D')
    legend('location','northwest')
    xlim([y(1) y(end)])
    ylim([-5 2])
    xlabel('Day')
    ylabel('log10 Biomass (g m^-^2)')
    title(['testcase v1 ' lname])
    stamp(cfile)
    print('-dpng',[fpath sname harv '_timeseries_Logmean_biomass_types_' loc '.png'])
     
    %% Temp habitat
    figure(3)
    clf
    plot(y,(squeeze(T_habitat(:,s,:))),'Linewidth',1); hold on;
    legend(stages(1:8))
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('Thab (^oC)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_Thab_' lname '.png'])
    
    %% nu - energy for biomass production
    figure(4)
    clf
    plot(y,(squeeze(energy_avail_rate(:,s,:))),'Linewidth',1); hold on;
    legend(stages(1:8))
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('energy avail rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_nu_' lname '.png'])
    
    %% Production (= nu * biom)
    figure(5)
    clf
    plot(y,(squeeze(energy_avail_rate(:,s,:).*biomass(:,s,:))),'Linewidth',1); hold on;
    legend(stages(1:8))
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('production (g)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_prod_' lname '.png'])
    
    %% Reproduction
    figure(6)
    clf
    subplot(2,1,1)
    plot(y,(squeeze(reproduction_rate(:,s,4))),'Linewidth',1,'color',cm10(4,:)); hold on;
    legend('MF')
    legend('location','west')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    ylabel('repro rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    
    subplot(2,1,2)
    plot(y,(squeeze(reproduction_rate(:,s,7))),'Linewidth',1,'color',cm10(7,:)); hold on;
    plot(y,(squeeze(reproduction_rate(:,s,8))),'Linewidth',1,'color',cm10(8,:)); hold on;
    legend(stages(7:8))
    legend('location','northwest')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('repro rate (g g^-^1 m^-^2)')
    stamp(sname)
    print('-dpng',[fpath sname '_ts_repro_' lname '.png'])
    
    %% Metabolism
    figure(7)
    clf
    plot(y,(squeeze(metabolism_rate(:,s,:))),'Linewidth',1); hold on;
    legend(stages(1:8))
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('metabolism rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_metab_' lname '.png'])
    
    %% Ingestion
    figure(8)
    clf
    plot(y,(squeeze(ingestion_rate(:,s,:))),'Linewidth',1); hold on;
    legend(stages(1:8))
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-0.05 0.25])
    xlabel('Day')
    ylabel('ingestion rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_ingest_' lname '.png'])
    
    %% Growth = gamma
    figure(9)
    clf
    plot(y,(squeeze(growth_rate(:,s,:))),'Linewidth',1); hold on;
    legend(stages(1:8))
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-0.01 0.15])
    xlabel('Day')
    ylabel('growth rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_gamma_' lname '.png'])
    
    %% Recruitment
    figure(10)
    clf
    subplot(3,1,1)
    plot(y,(squeeze(recruitment_flux(:,s,1))),'Linewidth',1,'color',cm10(1,:)); hold on;
    legend(stages(1))
    legend('location','west')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    ylabel('recruitment flux (g m^-^2)')
    title(['testcase v1 ' lname])
    
    subplot(3,1,2)
    for n=2:3
    plot(y,(squeeze(recruitment_flux(:,s,n))),'Linewidth',1,'color',cm10(n,:)); hold on;
    end
    legend(stages(2:3))
    legend('location','west')
    xlim([y(1) y(end)])
    %ylim([-0.2e-9 1.2e-9])
    ylabel('recruitment flux (g m^-^2)')
    
    subplot(3,1,3)
    for n=4:8
    plot(y,(squeeze(recruitment_flux(:,s,n))),'Linewidth',1,'color',cm10(n,:)); hold on;
    end
    legend(stages(4:8))
    legend('location','northwest')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('recruitment flux (g m^-^2)')
    stamp(sname)
    print('-dpng',[fpath sname '_ts_recruit_' lname '.png'])
    
    %% Predation flux & rate
    figure(11)
    clf
    plot(y,(squeeze(predation_flux(:,s,:))),'Linewidth',1); hold on;
    legend(stages(1:8))
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('predation flux (g m^-^2)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_predflux_' lname '.png'])
    
    figure(12)
    clf
    plot(y,(squeeze(predation_rate(:,s,:))),'Linewidth',1); hold on;
    legend(stages(1:8))
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('predation rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_predrate_' lname '.png'])
    
    %% Nat mort
    figure(13)
    clf
    plot(y,(squeeze(mortality_rate(:,s,:))),'Linewidth',1); hold on;
    legend(stages(1:8))
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([2.7e-04 2.8e-04])
    xlabel('Day')
    ylabel('natural mortality rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_nmort_' lname '.png'])
    
    %% Fish catch
    figure(14)
    clf
    plot(y,(squeeze(fish_catch_rate(:,s,:))),'Linewidth',1); hold on;
    legend(stages(1:8))
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('fishing rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_fishing_' lname '.png'])
    
    %% Total mortality w/o fishing
    figure(15)
    clf
    plot(y,(squeeze(tot1(:,s,:))),'Linewidth',1); hold on;
    legend(stages(1:8))
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('total non-fishing mortality rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_prod_' lname '.png'])
    
    %% Total mortality w/ fishing
    figure(16)
    clf
    plot(y,(squeeze(tot2(:,s,:))),'Linewidth',1); hold on;
    legend(stages(1:8))
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('total mortality rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_prod_' lname '.png'])
    
    %% Gross growth efficiency (= nu/consump)
    figure(17)
    clf
    subplot(3,1,1)
    plot(y,(squeeze(gge(:,s,1:5))),'Linewidth',1); hold on;
    legend(stages(1:5))
    legend('location','east')
    xlim([y(1) y(end)])
    ylim([-0.2 1])
    ylabel('gross growth effic')
    title(['testcase v1 ' lname])
    
    subplot(3,1,2)
    plot(y,(squeeze(gge(:,s,6))),'Linewidth',1,'color',cm10(6,:)); hold on;
    plot(y,(squeeze(gge(:,s,8))),'Linewidth',1,'color',cm10(8,:)); hold on;
    legend([stages(6) stages(8)])
    legend('location','east')
    xlim([y(1) y(end)])
    ylim([-0.2 1])
    ylabel('gross growth effic')
    
    subplot(3,1,3)
    plot(y,(squeeze(gge(:,s,7))),'Linewidth',1,'color',cm10(7,:)); hold on;
    legend(stages(7))
    legend('location','east')
    xlim([y(1) y(end)])
    ylim([-0.2 1])
    xlabel('Day')
    ylabel('gross growth effic')
    stamp(sname)
    print('-dpng',[fpath sname '_ts_gge_' lname '.png'])
    
end

%% Loop over types

for s=1:nid
    
    %% TIME SERIES ----------------------------------------------------
    loc = stages{s};
    lname = stages{s};
    y=t;
    
    %% Biomass
    figure(1)
    clf
    subplot(3,3,1)
    plot(y,log10(squeeze(biomass(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('log10 Biomass (g m^-^2)')
    title(['testcase v1 ' lname])
    
    subplot(3,3,2)
    plot(y,log10(squeeze(biomass(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    ylim([1 3])
    legend(locs)
    
%     subplot(3,3,2)
%     plot(y,log10(B(:,:,s)),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
%     plot(y,log10(F(:,:,s)),'r','Linewidth',2); hold on;
%     plot(y,log10(P(:,:,s)),'b','Linewidth',2); hold on;
%     plot(y,log10(D(:,:,s)),'k','Linewidth',2); hold on;
%     legend('B','F','P','D')
%     legend('location','northwest')
%     xlim([y(1) y(end)])
%     ylim([-5 2])
%     xlabel('Day')
%     ylabel('log10 Biomass (g m^-^2)')
%     title(['testcase v1 ' lname])
     
    %% Temp habitat
    subplot(3,3,3)
    plot(y,(squeeze(T_habitat(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('Thab (^oC)')
    title(['testcase v1 ' lname])
    
    %% nu - energy for biomass production
    subplot(3,3,4)
    plot(y,(squeeze(energy_avail_rate(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('energy avail rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    
    %% Production (= nu * biom)
    subplot(3,3,5)
    plot(y,(squeeze(energy_avail_rate(:,:,s).*biomass(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('production (g)')
    title(['testcase v1 ' lname])
    
    %% Reproduction
    subplot(3,3,6)
    plot(y,(squeeze(reproduction_rate(:,:,s))),'Linewidth',1,'color',cm10(4,:)); hold on;
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    ylabel('repro rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    
    %% Metabolism
    subplot(3,3,7)
    plot(y,(squeeze(metabolism_rate(:,:,s))),'Linewidth',1); hold on;
    legend(locs)
    %legend('location','eastoutside')
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('metabolism rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    
    %% Ingestion
    subplot(3,3,8)
    plot(y,(squeeze(ingestion_rate(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    ylim([-0.05 0.25])
    xlabel('Day')
    ylabel('ingestion rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    
    %% Growth = gamma
    subplot(3,3,9)
    plot(y,(squeeze(growth_rate(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    ylim([-0.01 0.15])
    xlabel('Day')
    ylabel('growth rate (g g^-^1 m^-^2)')
    title(['testcase v1 ' lname])
    stamp(sname)
    print('-dpng',[fpath sname '_ts_gamma_' lname '.png'])
    
    %% Recruitment
    figure(10)
    clf
    subplot(3,3,1)
    plot(y,(squeeze(recruitment_flux(:,:,s))),'Linewidth',1,'color',cm10(1,:)); hold on;
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    ylabel('recruitment flux (g m^-^2)')
    
    subplot(3,3,2)
    plot(y,(squeeze(predation_flux(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('predation flux (g m^-^2)')
    title(['testcase v1 ' lname])
    
    subplot(3,3,3)
    plot(y,(squeeze(predation_rate(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('predation rate (g g^-^1 m^-^2)')
    
    %% Nat mort
    subplot(3,3,4)
    plot(y,(squeeze(mortality_rate(:,:,s))),'Linewidth',1); hold on;
    legend(locs)
    xlim([y(1) y(end)])
    ylim([2.7e-04 2.8e-04])
    xlabel('Day')
    ylabel('natural mortality rate (g g^-^1 m^-^2)')
    
    %% Fish catch
    subplot(3,3,5)
    plot(y,(squeeze(fish_catch_rate(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('fishing rate (g g^-^1 m^-^2)')
    
    %% Total mortality w/o fishing
    tot1 = mortality_rate + predation_rate;
    
    subplot(3,3,6)
    plot(y,(squeeze(tot1(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('total non-fishing mortality rate (g g^-^1 m^-^2)')
    
    %% Total mortality w/ fishing
    
    subplot(3,3,7)
    plot(y,(squeeze(tot2(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    xlabel('Day')
    ylabel('total mortality rate (g g^-^1 m^-^2)')
    
    %% Gross growth efficiency (= nu/consump)
    subplot(3,3,8)
    plot(y,(squeeze(gge(:,:,s))),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    %ylim([-5.5 0.5])
    ylabel('gross growth effic')
    stamp(sname)
    print('-dpng',[fpath sname '_ts_gge_' lname '.png'])
    
end

%%  SPACE & TIME ----------------------------------------------------

% Biomass
f21 = figure(21);
pcolor(T_pelagic')
shading flat
ylabel('time')
xlabel('loc')
colorbar