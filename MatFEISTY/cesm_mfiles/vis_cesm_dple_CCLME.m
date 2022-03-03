% FEISTY DPLE outputs saved as NetCDF
% Initialized 2015
% CCLME only

clear all
close all

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

ID = GRD.ID;

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
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%pick year
StartYr = 2015;
%loop over members
submem = 1:40;
%%
for mem=1%:length(submem) %will loop over
    %%
    Member = submem(mem);
    exper = ['v14_Y' num2str(StartYr) '_M' num2str(Member) '_All_fish03_' ];
    
    load([fpath 'CESM_DPLE_CCLME_outputs_monthly_' exper cfile '.mat'])
    
    AllFish = AllF + AllP + AllD;
    [ni,nj,nt] = size(AllF);
    
    %% Create animation of distribution over time
    figure(1)
    h = imagescn(lat2,lon2,AllF(:,:,1));
    cb = colorbar;
    ylabel(cb,'F abund (g m^-^2)')
    cmocean matter
    title(datestr(datenum(0,0,1),'mmmm'))
    caxis([0 5])
    hold on
    gif(['CESM_DPLE_',exper,'CCLME_F_distr.gif'],'frame',gcf,'delaytime',1/12,'nodither')
    for k=2:nt
        h.CData = AllF(:,:,k);
        title(datestr(datenum(0,k,1),'mmmm'))
        gif
    end
    
    figure(2)
    h = imagescn(lat2,lon2,AllP(:,:,1));
    cb = colorbar;
    ylabel(cb,'P abund (g m^-^2)')
    cmocean matter
    title(datestr(datenum(0,0,1),'mmmm'))
    caxis([0 4])
    hold on
    gif(['CESM_DPLE_',exper,'CCLME_P_distr.gif'],'frame',gcf,'delaytime',1/12,'nodither')
    for k=2:nt
        h.CData = AllP(:,:,k);
        title(datestr(datenum(0,k,1),'mmmm'))
        gif
    end
    
    figure(3)
    h = imagescn(lat2,lon2,AllD(:,:,1));
    cb = colorbar;
    ylabel(cb,'D abund (g m^-^2)')
    cmocean matter
    title(datestr(datenum(0,0,1),'mmmm'))
    caxis([0 2])
    hold on
    gif(['CESM_DPLE_',exper,'CCLME_D_distr.gif'],'frame',gcf,'delaytime',1/12,'nodither')
    for k=2:nt
        h.CData = AllD(:,:,k);
        title(datestr(datenum(0,k,1),'mmmm'))
        gif
    end
    
    figure(4)
    h = imagescn(lat2,lon2,log10(AllB(:,:,1)));
    cb = colorbar;
    ylabel(cb,'B log_1_0 abund (g m^-^2)')
    cmocean matter
    title(datestr(datenum(0,0,1),'mmmm'))
    caxis([-1 1])
    hold on
    gif(['CESM_DPLE_',exper,'CCLME_B_distr.gif'],'frame',gcf,'delaytime',1/12,'nodither')
    for k=2:nt
        h.CData = AllB(:,:,k);
        title(datestr(datenum(0,k,1),'mmmm'))
        gif
    end


%     axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%         'Grid','off','FLineWidth',1)
%     surfm(TLAT,TLONG,log10(Zb))
%     colormap(cmBP50)                %decent looking coastlines
%     h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%     caxis([-1 2]);
%     hcb = colorbar('h');
%     set(gcf,'renderer','painters')
%     title('DPLE log10 mean benthic biomass (g m^-^2)')
%     stamp(exper)
end