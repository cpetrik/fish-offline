% Plot ts of LME means of FEISTY inputs & outputs w/ climate anoms
% CESM FOSI

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/' cfile '/corrs/'];
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

mod = 'v15_All_fish03_';

%% FOSI input forcing
% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],'adety',...
    'atb','atp','azlosy','lme_dety_stda','lme_tb_stda','lme_tp_stda',...
    'lme_mzly_stda');

%% Fish data
% Anoms of area-weighted means with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'lme_ann_mean_anoms.mat'],'aa','va')
Ball = aa;
Bstd = va;
clear aa va

load([fpath 'FEISTY_FOSI_',mod,'lme_nu_ann_mean_anoms.mat'],'aa','va')
Pall = aa;
Pstd = va;
clear aa va

load([ypath 'FishMIP_Phase3a_LME_CPUE_1961-2010_ann_mean_anoms.mat'],'aall',...
    'vall')
Call = aall;
Cstd = vall;
clear aall vall

load([ypath 'FishMIP_Phase3a_LME_catch_1961-2010_ann_mean_anoms.mat'],...
    'aall','vall')
Mall = aall;
Mstd = vall;
clear aall vall

% Ecosystem type
%load([fpath 'LME_fosi_fished_',mod,cfile '.mat'],'etype','etex','Ebiom');

% Clusters
%load([spath,'LME_biom_nu_cpue_cme_A_mlr_coeffs_reduc_alllags3_R2_cluster.mat'])
load([spath,'LME_biom_nu_cpue_catch_A_mlr_cluster_ind_drivers.mat'])

%%
fyr = 1948:2015;
cyr = 1961:2010;
[yr,fid] = intersect(fyr,cyr);

%% subset effort years
adety  = adety(:,fid);
atb    = atb(:,fid);
atp    = atp(:,fid);
azlosy = azlosy(:,fid);
Ball   = Ball(:,fid);
Pall   = Pall(:,fid);

%% divide by 2 std
adety  = adety ./ repmat(2*lme_dety_stda,1,50);
atb    = atb ./ repmat(2*lme_tb_stda,1,50);
atp    = atp ./ repmat(2*lme_tp_stda,1,50);
azlosy = azlosy ./ repmat(2*lme_mzly_stda,1,50);
Ball   = Ball ./ repmat(2*Bstd,1,50);
Pall   = Pall ./ repmat(2*Pstd,1,50);
Call   = Call ./ repmat(2*Cstd,1,50);
Mall   = Mall ./ repmat(2*Mstd,1,50);

%%
% LMEs
lid = [54,1:2,10,3,5:7,65]; %ADD 65 = Aleutian Islands
lname = {'CHK','EBS','GAK','HI','CCE','GMX','SE','NE','AI'};

% pid = [54,1:2,3];
% pname = {'CHK','EBS','GAK','CCE'};
%
% aid = [7,6,5];
% aname = {'NE','SE','GMX'};


%% colors
% cm10=[0.5 0.5 0;... %tan/army
%     0 0.7 0;...   %g
%     1 0 1;...     %m
%     1 0 0;...     %r
%     0.5 0 0;...   %maroon
%     0/255 206/255 209/255;... %turq
%     0 0.5 0.75;...   %med blue
%     0 0 0.75;...    %b
%     0.5 0.5 0.5; ...    %med grey
%     0 0 0];...      %black

% cm10=[0.5 0.5 0.5; ... %med grey
%     0 0 0; ...      %black
%     0 0 0.75;...    %b
%     1 0 0;...       %r
%     0 0.5 0.75;...   %med blue
%     0.5 0 0;...   %maroon
%     0/255 206/255 209/255;... %turq
%     1 0 1;...     %m
%     0 0.7 0;...   %g
%     ];

cm10 = [
    34/255 136/255 51/255;...   %green
    170/255 51/255 119/255;...  %purple
    238/255 102/255 119/255;... %red
    0/255 68/255 136/255;...    %blue
    51/255 187/255 238/255;...  %cyan
    153/255 153/255 51/255;...  %olive
    0 0 0;...                   %black
    0.50 0.50 0.50;...          % grey
    ];

set(groot,'defaultAxesColorOrder',cm10);

load('paul_tol_cmaps.mat')

%try muted and add greys
mcol = muted ./ 255;
mcol = mcol(1:8,:);

%% plot time series  =================================
%

for i=1:length(lid)
    lme = lid(i);

    figure(1)
    clf
    tiledlayout(4,4, 'TileSpacing', 'compact')

    nexttile % TP-biom
    %driver
    yyaxis left
    plot(yr,atp(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    ylabel('TP')
    %fish
    yyaxis right
    plot(yr,Ball(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    title('Biom')

    nexttile % TP-nu
    %driver
    yyaxis left
    plot(yr,atp(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('TP')
    %fish
    yyaxis right
    plot(yr,Pall(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    title([lname(i); 'Prod'])

    nexttile % TP-CPUE
    %driver
    yyaxis left
    plot(yr,atp(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('TP')
    %fish
    yyaxis right
    plot(yr,Call(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    title('CPUE')

    nexttile % TP-CME
    %driver
    yyaxis left
    plot(yr,atp(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('TP')
    %fish
    yyaxis right
    plot(yr,Mall(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    title('CME')

    nexttile % TB-biom
    %driver
    yyaxis left
    plot(yr,atb(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    ylabel('TB')
    %fish
    yyaxis right
    plot(yr,Ball(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('Biom')

    nexttile % TB-nu
    %driver
    yyaxis left
    plot(yr,atb(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('TB')
    %fish
    yyaxis right
    plot(yr,Pall(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('Prod')

    nexttile % TB-CPUE
    %driver
    yyaxis left
    plot(yr,atb(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('TB')
    %fish
    yyaxis right
    plot(yr,Call(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('CPUE')

    nexttile % TB-CME
    %driver
    yyaxis left
    plot(yr,atb(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('TB')
    %fish
    yyaxis right
    plot(yr,Mall(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('CME')

    nexttile % ZL-biom
    %driver
    yyaxis left
    plot(yr,azlosy(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    ylabel('ZL')
    %fish
    yyaxis right
    plot(yr,Ball(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('Biom')

    nexttile % ZL-nu
    %driver
    yyaxis left
    plot(yr,azlosy(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('ZL')
    %fish
    yyaxis right
    plot(yr,Pall(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('Prod')

    nexttile % ZL-CPUE
    %driver
    yyaxis left
    plot(yr,azlosy(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('ZL')
    %fish
    yyaxis right
    plot(yr,Call(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('CPUE')

    nexttile % ZL-CME
    %driver
    yyaxis left
    plot(yr,azlosy(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('ZL')
    %fish
    yyaxis right
    plot(yr,Mall(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('CME')

    nexttile % Det-biom
    %driver
    yyaxis left
    plot(yr,adety(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    ylabel('Det')
    %fish
    yyaxis right
    plot(yr,Ball(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('Biom')

    nexttile % Det-nu
    %driver
    yyaxis left
    plot(yr,adety(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('Det')
    %fish
    yyaxis right
    plot(yr,Pall(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('Prod')

    nexttile % Det-CPUE
    %driver
    yyaxis left
    plot(yr,adety(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('Det')
    %fish
    yyaxis right
    plot(yr,Call(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('CPUE')

    nexttile % Det-CME
    %driver
    yyaxis left
    plot(yr,adety(lme,:),'color',[0 0.25 0.75]); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('Det')
    %fish
    yyaxis right
    plot(yr,Mall(lme,:),'-k'); hold on;
    xlim([yr(1) yr(end)])
    %ylabel('CME')

    print('-dpng',[ppath 'ts_drivers_FOSI_',mod,'biom_nu_cpue_cme_',lname{i},'.png'])

end
% lg = legend(nexttile(4),{'PDO','Sm','Md','Lg'});
% lg.Location = 'southoutside';
% lg.Orientation = 'horizontal';

%% Plot all LMEs by cluster against main driver(s)

% nrows=11;
% ncols=6;
% pos = subfigrid(nrows,ncols,[0.05 0.025 0.025 0.05],[0.72 0.75]);
nrows=7;
ncols=9;
pos = subfigrid(nrows,ncols,[0.05 0.025 0.025 0.05],[0.72 0.75]);

%% Biomass
%figure('Units','inches','Position',[1 3 6.5 8.5]);
figure('Units','inches','Position',[1 3 9 6.5]);
i=0;
%while(i<64)
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        iid = LIDs(i);
        subplot('position',pos(m,:,n))

        %driver
        yyaxis left
        if (ClusterB(i,5)==1)
            plot(yr,adety(iid,:),'color',mcol(1,:)); hold on;
        end
        if (ClusterB(i,5)==2)
            plot(yr,adety(iid,:),'color',mcol(2,:)); hold on;
        end
        if (ClusterB(i,4)==3)
            plot(yr,azlosy(iid,:),'color',mcol(3,:)); hold on;
        end
        if (ClusterB(i,4)==4)
            plot(yr,azlosy(iid,:),'color',mcol(4,:)); hold on;
        end
        if (ClusterB(i,3)==5)
            plot(yr,atb(iid,:),'color',mcol(5,:)); hold on;
        end
        if (ClusterB(i,3)==6)
            plot(yr,atb(iid,:),'color',mcol(6,:)); hold on;
        end
        if (ClusterB(i,2)==7)
            plot(yr,atp(iid,:),'color',mcol(7,:)); hold on;
        end
        if (ClusterB(i,2)==8)
            plot(yr,atp(iid,:),'color',mcol(8,:)); hold on;
        end
        xlim([yr(1) yr(end)])
        %fish
        yyaxis right
        plot(yr,Ball(iid,:),'-k'); hold on;
        xlim([yr(1) yr(end)])
        %title('Biom')
        title(num2str(iid))

        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
%end
%stamp(mod)
print('-dpng',[ppath 'ts_FOSI_',mod,'biom_mlr_cluster_drivers.png'])

%% Prod
figure('Units','inches','Position',[1 3 9 6.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        iid = LIDs(i);
        subplot('position',pos(m,:,n))

        %driver
        yyaxis left
        if (ClusterP(i,5)==1)
            plot(yr,adety(iid,:),'color',mcol(1,:)); hold on;
        end
        if (ClusterP(i,5)==2)
            plot(yr,adety(iid,:),'color',mcol(2,:)); hold on;
        end
        if (ClusterP(i,4)==3)
            plot(yr,azlosy(iid,:),'color',mcol(3,:)); hold on;
        end
        if (ClusterP(i,4)==4)
            plot(yr,azlosy(iid,:),'color',mcol(4,:)); hold on;
        end
        if (ClusterP(i,3)==5)
            plot(yr,atb(iid,:),'color',mcol(5,:)); hold on;
        end
        if (ClusterP(i,3)==6)
            plot(yr,atb(iid,:),'color',mcol(6,:)); hold on;
        end
        if (ClusterP(i,2)==7)
            plot(yr,atp(iid,:),'color',mcol(7,:)); hold on;
        end
        if (ClusterP(i,2)==8)
            plot(yr,atp(iid,:),'color',mcol(8,:)); hold on;
        end
        xlim([yr(1) yr(end)])
        %fish
        yyaxis right
        plot(yr,Pall(iid,:),'-k'); hold on;
        xlim([yr(1) yr(end)])
        %title('Prod')
        title(num2str(iid))

        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_FOSI_',mod,'nu_mlr_cluster_drivers.png'])

%% CPUE
figure('Units','inches','Position',[1 3 9 6.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        iid = LIDs(i);
        subplot('position',pos(m,:,n))

        %driver
        yyaxis left
        if (ClusterC(i,5)==1)
            plot(yr,adety(iid,:),'color',mcol(1,:)); hold on;
        end
        if (ClusterC(i,5)==2)
            plot(yr,adety(iid,:),'color',mcol(2,:)); hold on;
        end
        if (ClusterC(i,4)==3)
            plot(yr,azlosy(iid,:),'color',mcol(3,:)); hold on;
        end
        if (ClusterC(i,4)==4)
            plot(yr,azlosy(iid,:),'color',mcol(4,:)); hold on;
        end
        if (ClusterC(i,3)==5)
            plot(yr,atb(iid,:),'color',mcol(5,:)); hold on;
        end
        if (ClusterC(i,3)==6)
            plot(yr,atb(iid,:),'color',mcol(6,:)); hold on;
        end
        if (ClusterC(i,2)==7)
            plot(yr,atp(iid,:),'color',mcol(7,:)); hold on;
        end
        if (ClusterC(i,2)==8)
            plot(yr,atp(iid,:),'color',mcol(8,:)); hold on;
        end
        xlim([yr(1) yr(end)])
        %fish
        yyaxis right
        plot(yr,Call(iid,:),'-k'); hold on;
        xlim([yr(1) yr(end)])
        %title('CPUE')
        title(num2str(iid))

        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_FOSI_',mod,'cpue_mlr_cluster_drivers.png'])

%% CME
figure('Units','inches','Position',[1 3 9 6.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        iid = LIDs(i);
        subplot('position',pos(m,:,n))

        %driver
        yyaxis left
        if (ClusterM(i,5)==1)
            plot(yr,adety(iid,:),'color',mcol(1,:)); hold on;
        end
        if (ClusterM(i,5)==2)
            plot(yr,adety(iid,:),'color',mcol(2,:)); hold on;
        end
        if (ClusterM(i,4)==3)
            plot(yr,azlosy(iid,:),'color',mcol(3,:)); hold on;
        end
        if (ClusterM(i,4)==4)
            plot(yr,azlosy(iid,:),'color',mcol(4,:)); hold on;
        end
        if (ClusterM(i,3)==5)
            plot(yr,atb(iid,:),'color',mcol(5,:)); hold on;
        end
        if (ClusterM(i,3)==6)
            plot(yr,atb(iid,:),'color',mcol(6,:)); hold on;
        end
        if (ClusterM(i,2)==7)
            plot(yr,atp(iid,:),'color',mcol(7,:)); hold on;
        end
        if (ClusterM(i,2)==8)
            plot(yr,atp(iid,:),'color',mcol(8,:)); hold on;
        end
        xlim([yr(1) yr(end)])
        %fish
        yyaxis right
        plot(yr,Mall(iid,:),'-k'); hold on;
        xlim([yr(1) yr(end)])
        %title('CME')
        title(num2str(iid))

        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_FOSI_',mod,'cme_mlr_cluster_drivers.png'])

%% All fish ts together by LME
ecol = drainbow(2:3:end,:) ./ 255;

%1st half
figure('Units','inches','Position',[1 3 9 6.5]);
i=0;
for m = 1:2:6
    for n = 1:ncols
        i=i+1;
        iid = LIDs(i);
        
        subplot('position',pos(m,:,n))
        %biom
        yyaxis left
        plot(yr,Ball(iid,:),'-','color',cm10(1,:)); hold on;
        xlim([yr(1) yr(end)])
        %nu
        yyaxis right
        plot(yr,Pall(iid,:),'-','color',cm10(2,:)); hold on;
        xlim([yr(1) yr(end)])
        title(num2str(iid))

        subplot('position',pos(m+1,:,n))
        %cpue
        yyaxis left
        plot(yr,Call(iid,:),'-','color',cm10(4,:)); hold on;
        xlim([yr(1) yr(end)])
        %cme
        yyaxis right
        plot(yr,Mall(iid,:),'-','color',cm10(7,:)); hold on;
        xlim([yr(1) yr(end)])
        title(num2str(iid))

        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_FOSI_',mod,'biom_nu_cpue_cme_LMEs1-27.png'])

%2nd half
figure('Units','inches','Position',[1 3 9 6.5]);
i=27;
for m = 1:2:6
    for n = 1:ncols
        i=i+1;
        iid = LIDs(i);
        
        subplot('position',pos(m,:,n))
        %biom
        yyaxis left
        plot(yr,Ball(iid,:),'-','color',cm10(1,:)); hold on;
        xlim([yr(1) yr(end)])
        %nu
        yyaxis right
        plot(yr,Pall(iid,:),'-','color',cm10(2,:)); hold on;
        xlim([yr(1) yr(end)])
        title(num2str(iid))

        subplot('position',pos(m+1,:,n))
        %cpue
        yyaxis left
        plot(yr,Call(iid,:),'-','color',cm10(4,:)); hold on;
        xlim([yr(1) yr(end)])
        %cme
        yyaxis right
        plot(yr,Mall(iid,:),'-','color',cm10(7,:)); hold on;
        xlim([yr(1) yr(end)])
        title(num2str(iid))

        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_FOSI_',mod,'biom_nu_cpue_cme_LMEs28-54.png'])

%% last
figure('Units','inches','Position',[1 3 9 6.5]);
i=54;
for m = 1
    for n = 1:ncols
        i=i+1;
        iid = LIDs(i);
        
        subplot('position',pos(m,:,n))
        %biom
        yyaxis left
        plot(yr,Ball(iid,:),'-','color',cm10(1,:)); hold on;
        xlim([yr(1) yr(end)])
        %nu
        yyaxis right
        plot(yr,Pall(iid,:),'-','color',cm10(2,:)); hold on;
        xlim([yr(1) yr(end)])
        title(num2str(iid))

        subplot('position',pos(m+1,:,n))
        %cpue
        yyaxis left
        plot(yr,Call(iid,:),'-','color',cm10(4,:)); hold on;
        xlim([yr(1) yr(end)])
        %cme
        yyaxis right
        plot(yr,Mall(iid,:),'-','color',cm10(7,:)); hold on;
        xlim([yr(1) yr(end)])
        title(num2str(iid))

        if i==63
            subplot('position',pos(m+2,:,n))
            plot(yr,4*ones(size(yr)),'-','color',cm10(1,:)); hold on;
            plot(yr,3*ones(size(yr)),'-','color',cm10(2,:)); hold on;
            plot(yr,2*ones(size(yr)),'-','color',cm10(4,:)); hold on;
            plot(yr,ones(size(yr)),'-','color',cm10(7,:)); hold on;
            legend('Biom','Prod','CPUE','CME')
        end
    end
end
print('-dpng',[ppath 'ts_FOSI_',mod,'biom_nu_cpue_cme_LMEs55-66.png'])


