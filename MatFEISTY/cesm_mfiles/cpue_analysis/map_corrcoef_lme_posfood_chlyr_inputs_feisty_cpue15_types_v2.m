% Mapp corr coefs of ind driver combos with cpue
% Subplots together for comparison
% Restricted analysis to chl yrs 1997-2015
% Const effort
% Add biom & prod from obsfish

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regress_cpue/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',...
    cfile,'/corrs_cpue/'];

mod = 'v15_All_fish03';
mod2 = 'v15_obsfish2015';

%%
load([spath 'LMEs_corr_cpue_chlyrs15_inputs_feisty_mostsiglag_posfood.mat'])
%dim: LME x driver x (corr, p-val, lag)

LAmat2 = LAmat;
LFmat2 = LFmat;
LPmat2 = LPmat;
LDmat2 = LDmat;

clear LAmat LFmat LPmat LDmat

%%
load([spath 'LMEs_corr_cpue_chlyrs_inputs_obsfish2015_mostsiglag_posfood.mat'])

LAmat2(:,9:10,:) = LAtab(:,7:8,:);
LFmat2(:,9:10,:) = LFtab(:,7:8,:);
LPmat2(:,9:10,:) = LPtab(:,7:8,:);
LDmat2(:,9:10,:) = LDtab(:,7:8,:);

tanom2 = {'TP','TB','Det','ZmLoss','SST','Chl','CBiom','CProd','OBiom','OProd'};

%%
clear LAtab LFtab LPtab LDtab

%% Map data
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
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

cmRB = cbrewer('div','RdBu',21,'PCHIP');
cmRB = flipud(cmRB);

subpos = [0.015 0.6 0.43 0.2;... 
    0.015 0.4 0.43 0.2;... 
    0.48 0.4 0.43 0.2;... 
    0.48 0.6 0.43 0.2;... 
    0.015 0.8 0.43 0.2;... 
    0.48 0.8 0.43 0.2;... 
    0.015 0.2 0.43 0.2;... 
    0.015 0.0 0.43 0.2;... 
    0.48 0.2 0.43 0.2;... 
    0.48 0.0 0.43 0.2];   

%% All fish 
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
for j=1:10 
    Cmat = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (LAmat2(i,j,2) <= 0.05)
            Cmat(id) = LAmat2(i,j,1);
        end
    end

    subplot('Position',subpos(j,:))
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(TLAT,TLONG,Cmat)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 1])
    set(gcf,'renderer','painters')
    text(0,1.75,tanom2{j},'HorizontalAlignment','center')

end
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
print('-dpng',[ppath 'Map_LMEs_chlyr_cpue15_driver_feisty_obsfish2015_corrcoef_All.png'])

%% Forage
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
for j=1:10 
    Cmat = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (LFmat2(i,j,2) <= 0.05)
            Cmat(id) = LFmat2(i,j,1);
        end
    end

    subplot('Position',subpos(j,:))
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(TLAT,TLONG,Cmat)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 1])
    set(gcf,'renderer','painters')
    text(0,1.75,tanom2{j},'HorizontalAlignment','center')

end
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
print('-dpng',[ppath 'Map_LMEs_chlyr_cpue15_driver_feisty_obsfish2015_corrcoef_F.png'])

%% Lg Pel
f3 = figure('Units','inches','Position',[1 3 6.5 8]);
for j=1:10 
    Cmat = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (LPmat2(i,j,2) <= 0.05)
            Cmat(id) = LPmat2(i,j,1);
        end
    end

    subplot('Position',subpos(j,:))
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(TLAT,TLONG,Cmat)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 1])
    set(gcf,'renderer','painters')
    text(0,1.75,tanom2{j},'HorizontalAlignment','center')

end
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
print('-dpng',[ppath 'Map_LMEs_chlyr_cpue15_driver_feisty_obsfish2015_corrcoef_P.png'])

%% Demersal
f4 = figure('Units','inches','Position',[1 3 6.5 8]);
for j=1:10 
    Cmat = nan(ni,nj);

    for i=1:length(lid)
        L=lid(i);
        id = find(tlme==L);

        if (LDmat2(i,j,2) <= 0.05)
            Cmat(id) = LDmat2(i,j,1);
        end
    end

    subplot('Position',subpos(j,:))
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(TLAT,TLONG,Cmat)
    colormap(cmRB)
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    clim([-1 1])
    set(gcf,'renderer','painters')
    text(0,1.75,tanom2{j},'HorizontalAlignment','center')

end
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
print('-dpng',[ppath 'Map_LMEs_chlyr_cpue15_driver_feisty_obsfish2015_corrcoef_D.png'])


%% 8plot by driver
% f1 = figure('Units','inches','Position',[1 3 6.5 8]);
% %f1.Units = 'inches';
% 
% %A - F biom
% subplot('Position',[0.015 0.75 0.43 0.25])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(TLAT2,TLON2,log10(AllF2))
% colormap(cmBP50)
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-1.5 2.0])
% set(gcf,'renderer','painters')
% text(0,1.75,'F biom','HorizontalAlignment','center')
% 
% %B - P biom
% subplot('Position',[0.015 0.5 0.43 0.25])
% 
% %C - D biom
% subplot('Position',[0.015 0.25 0.43 0.25])
% 
% %D - A biom
% subplot('Position',[0.015 0.0 0.43 0.25])
% 
% %E - F prod
% subplot('Position',[0.51 0.75 0.43 0.25])
% 
% %F - P prod
% subplot('Position',[0.51 0.5 0.43 0.25])
% 
% %G - D prod
% subplot('Position',[0.51 0.25 0.43 0.25])
% 
% %H - All prod
% subplot('Position',[0.51 0.0 0.43 0.25])
