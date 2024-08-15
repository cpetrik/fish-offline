% Plot corr of driver-fish corrs
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/' cfile '/corrs/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

load([spath,'LMEs_corr_biom_drivers_0_5_lag.mat'])

%%  ---------------------------------------------------------
ctex = {'TP','TB','Det','Zmeso','ZmLoss'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

mcol = [238/255 102/255 119/255;... %red
    170/255 51/255 119/255;...  %purple
    34/255 136/255 51/255;...   %green
    0/255 68/255 136/255;...    %blue
    51/255 187/255 238/255;...  %cyan
    ];

%colorblind friendly
% cb=[34/255 136/255 51/255;...   %green
%     153/255 153/255 51/255;...  %olive
%     51/255 187/255 238/255;...  %cyan
%     0/255 68/255 136/255;...    %blue
%     238/255 102/255 119/255;... %red
%     170/255 51/255 119/255;...  %purple
%     0 0 0;...                   %black
%     0.25 0.25 0.25;...             %dk grey
%     0.50 0.50 0.50;...             % grey
%     0.75 0.75 0.75];               %lt grey

%% Map
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

tlme = double(lme_mask);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
clatlim=[plotminlat plotmaxlat];
clonlim=[plotminlon plotmaxlon];

load coastlines;

%%
% %significance
% Fsig  = ones(ni,nj);
% Psig  = ones(ni,nj);
% Dsig  = ones(ni,nj);
% Vsig  = ones(ni,nj);
% Bsig  = ones(ni,nj);
% Ssig  = ones(ni,nj);
% Msig  = ones(ni,nj);
% Lsig  = ones(ni,nj);  

% Loop over drivers
for j = 1:length(tanom)
    % on grid
    Fcorr  = nan(ni,nj);
    Pcorr  = nan(ni,nj);
    Dcorr  = nan(ni,nj);
    Acorr  = nan(ni,nj);
    Bcorr  = nan(ni,nj);
    Scorr  = nan(ni,nj);
    Mcorr  = nan(ni,nj);
    Lcorr  = nan(ni,nj);

    % Use 1yr lag for F, D, A; 2yr for P, 0yr for B

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        Fcorr(id) = squeeze(FtabC(j,2,i));
        Pcorr(id) = squeeze(PtabC(j,3,i));
        Dcorr(id) = squeeze(DtabC(j,2,i));
        Acorr(id) = squeeze(AtabC(j,2,i));
        Bcorr(id) = squeeze(BtabC(j,1,i));
        Scorr(id) = squeeze(StabC(j,1,i));
        Mcorr(id) = squeeze(MtabC(j,2,i));
        Lcorr(id) = squeeze(LtabC(j,3,i));
    end

    %%
    figure(1)
    clf
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    hold on
    surfm(TLAT,TLONG,Acorr)
    cmocean('balance')
    colorbar
    caxis([-1 1])
    title('All fishes')
    c=colorbar;
    c.Label.String = ['Corr with ' tanom{j}];
    print('-dpng',[ppath 'Map_LMEs_',tanom{j},'_corr_allfish.png'])

    %% same driver, diff fn types
    figure(2)
    clf
    subplot('Position',[0.01 0.575 0.32 0.4]) %F
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    hold on
    surfm(TLAT,TLONG,Fcorr)
    cmocean('balance')
    caxis([-1 1])
    title('Forage fishes')

    subplot('Position',[0.33 0.575 0.32 0.4]) %P
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    hold on
    surfm(TLAT,TLONG,Pcorr)
    cmocean('balance')
    caxis([-1 1])
    title('Large pelagics')

    subplot('Position',[0.65 0.575 0.32 0.4]) %A
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    hold on
    surfm(TLAT,TLONG,Acorr)
    cmocean('balance')
    caxis([-1 1])
    title('All fishes')

    subplot('Position',[0.01 0.10 0.32 0.4]) %D
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    hold on
    surfm(TLAT,TLONG,Dcorr)
    cmocean('balance')
    caxis([-1 1])
    title('Demersals')

    subplot('Position',[0.33 0.10 0.32 0.4]) %B
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    hold on
    surfm(TLAT,TLONG,Bcorr)
    cmocean('balance')
    caxis([-1 1])
    colorbar('Position',[0.7 0.2 0.03 0.25])
    title('Benthos')
    c=colorbar('Position',[0.7 0.2 0.03 0.25]);
    c.Label.String = ['Corr with ' tanom{j}];
    print('-dpng',[ppath 'Map_LMEs_',tanom{j},'_corr_fntypes.png'])

   
    %% all drivers, F
    f3=figure(3);
    if (j==1)
        subplot('Position',[0.01 0.575 0.32 0.4]); %Tp
    elseif (j==2)
        subplot('Position',[0.33 0.575 0.32 0.4]); %Zm
    elseif (j==3)
        subplot('Position',[0.65 0.575 0.32 0.4]); %ZmLoss
    elseif (j==4)
        subplot('Position',[0.01 0.10 0.32 0.4]); %Tb
    else
        subplot('Position',[0.33 0.10 0.32 0.4]); %Det
    end

    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    hold on
    surfm(TLAT,TLONG,Fcorr)
    cmocean('balance')
    caxis([-1 1])
    title(tanom{j})
    if j==5
        c=colorbar('Position',[0.7 0.2 0.03 0.25]);
        c.Label.String = 'Corr with F';
    end

    %% all drivers, P
    f4=figure(4);
    if (j==1)
        subplot('Position',[0.01 0.575 0.32 0.4]); %Tp
    elseif (j==2)
        subplot('Position',[0.33 0.575 0.32 0.4]); %Zm
    elseif (j==3)
        subplot('Position',[0.65 0.575 0.32 0.4]); %ZmLoss
    elseif (j==4)
        subplot('Position',[0.01 0.10 0.32 0.4]); %Tb
    else
        subplot('Position',[0.33 0.10 0.32 0.4]); %Det
    end

    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    hold on
    surfm(TLAT,TLONG,Pcorr)
    cmocean('balance')
    caxis([-1 1])
    title(tanom{j})
    if j==5
        c=colorbar('Position',[0.7 0.2 0.03 0.25]);
        c.Label.String = 'Corr with P';
    end

    %% all drivers, D
    f5=figure(5);
    if (j==1)
        subplot('Position',[0.01 0.575 0.32 0.4]); %Tp
    elseif (j==2)
        subplot('Position',[0.33 0.575 0.32 0.4]); %Zm
    elseif (j==3)
        subplot('Position',[0.65 0.575 0.32 0.4]); %ZmLoss
    elseif (j==4)
        subplot('Position',[0.01 0.10 0.32 0.4]); %Tb
    else
        subplot('Position',[0.33 0.10 0.32 0.4]); %Det
    end

    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    hold on
    surfm(TLAT,TLONG,Dcorr)
    cmocean('balance')
    caxis([-1 1])
    title(tanom{j})
    if j==5
        c=colorbar('Position',[0.7 0.2 0.03 0.25]);
        c.Label.String = 'Corr with D';
    end

    %% all drivers, A
    f6=figure(6);
    if (j==1)
        subplot('Position',[0.01 0.575 0.32 0.4]); %Tp
    elseif (j==2)
        subplot('Position',[0.33 0.575 0.32 0.4]); %Zm
    elseif (j==3)
        subplot('Position',[0.65 0.575 0.32 0.4]); %ZmLoss
    elseif (j==4)
        subplot('Position',[0.01 0.10 0.32 0.4]); %Tb
    else
        subplot('Position',[0.33 0.10 0.32 0.4]); %Det
    end

    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    hold on
    surfm(TLAT,TLONG,Acorr)
    cmocean('balance')
    caxis([-1 1])
    title(tanom{j})
    if j==5
        c=colorbar('Position',[0.7 0.2 0.03 0.25]);
        c.Label.String = 'Corr with All';
    end

    %% all drivers, B
    f7=figure(7);
    if (j==1)
        subplot('Position',[0.01 0.575 0.32 0.4]); %Tp
    elseif (j==2)
        subplot('Position',[0.33 0.575 0.32 0.4]); %Zm
    elseif (j==3)
        subplot('Position',[0.65 0.575 0.32 0.4]); %ZmLoss
    elseif (j==4)
        subplot('Position',[0.01 0.10 0.32 0.4]); %Tb
    else
        subplot('Position',[0.33 0.10 0.32 0.4]); %Det
    end

    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    hold on
    surfm(TLAT,TLONG,Bcorr)
    cmocean('balance')
    caxis([-1 1])
    title(tanom{j})
    if j==5
        c=colorbar('Position',[0.7 0.2 0.03 0.25]);
        c.Label.String = 'Corr with B';
    end

end
print(f3,'-dpng',[ppath 'Map_LMEs_drivers_corr_F.png'])
print(f4,'-dpng',[ppath 'Map_LMEs_drivers_corr_P.png'])
print(f5,'-dpng',[ppath 'Map_LMEs_drivers_corr_D.png'])
print(f6,'-dpng',[ppath 'Map_LMEs_drivers_corr_A.png'])
print(f7,'-dpng',[ppath 'Map_LMEs_drivers_corr_B.png'])

