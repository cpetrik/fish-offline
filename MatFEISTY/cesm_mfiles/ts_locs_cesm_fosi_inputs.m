% Visualize FOSI forcing of FEISTY 
% from monthly netcdfs

clear 
close all

%% Paths

fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([fpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([fpath 'Data_grid_POP_gx1v6_noSeas.mat']);

cfile = 'FOSI';

%% FEISTY Inputs
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat'],...
    'FillValue','missing_value','TEMP_150m','TEMP_150m_units','TEMP_bottom',...
    'TEMP_bottom_units','POC_FLUX_IN_bottom','POC_FLUX_IN_bottom_units',...
    'TLAT','TLONG','TAREA','time','yr');
load([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo.mat'],...
    'LzooC_150m','Lzoo_loss_150m');

%% nans & zeros
TEMP_150m = double(TEMP_150m);
TEMP_bottom = double(TEMP_bottom);
POC_FLUX_IN_bottom = double(POC_FLUX_IN_bottom);

TEMP_bottom(TEMP_bottom >= 9.9e+36) = nan;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;

LzooC_150m(LzooC_150m<0) = 0.0;
Lzoo_loss_150m(Lzoo_loss_150m<0) = 0.0;
Lzoo_quad_150m(Lzoo_quad_150m<0) = 0.0;
POC_FLUX_IN_bottom(POC_FLUX_IN_bottom<0) = 0.0;

%% reshape
[ni,nj,nt] = size(TEMP_150m);
TEMP_150m = reshape(TEMP_150m,ni*nj,nt);
TEMP_bottom = reshape(TEMP_bottom,ni*nj,nt);
POC_FLUX_IN_bottom = reshape(POC_FLUX_IN_bottom,ni*nj,nt);
LzooC_150m = reshape(LzooC_150m,ni*nj,nt);
Lzoo_loss_150m = reshape(Lzoo_loss_150m,ni*nj,nt);
Lzoo_quad_150m = reshape(Lzoo_quad_150m,ni*nj,nt);

y = (1:nt)/12;

%% Locs
load([fpath 'cesm_grid_id_locs_area_dep.mat'],'ids','abbrev');
spots = abbrev;
ID = GRD.ID(ids);
NX = length(ID);
spots=spots';

%% convert units
%quad already in correct units
LzooC_150m = LzooC_150m * 1e-9 * 1e4 * 12.01 * 9.0;
Lzoo_loss_150m = Lzoo_loss_150m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
POC_FLUX_IN_bottom = POC_FLUX_IN_bottom * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

%%
for s=1%:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    lid = ID(s);
    
    %% TIME SERIES ----------------------------------------------------
    
    figure(1)
    clf
    subplot(3,2,1)
    plot(y,TEMP_150m(lid,:),'k','Linewidth',1.5); 
    xlim([y(1) y(end)])
    %ylim([-5 2])
    xlabel('Year')
    %ylabel('log10 Biomass (g m^-^2)')
    title([loc ' Tpel'])
    
    subplot(3,2,2)
    plot(y,TEMP_bottom(lid,:),'k','Linewidth',1.5); 
    xlim([y(1) y(end)])
    %ylim([-5 2])
    xlabel('Year')
    %ylabel('log10 Biomass (g m^-^2)')
    title([loc ' Tbtm'])
    
    subplot(3,2,3)
    plot(y,POC_FLUX_IN_bottom(lid,:),'k','Linewidth',1.5); 
    xlim([y(1) y(end)])
    %ylim([-5 2])
    xlabel('Year')
    ylabel('(g m^-^2 d^-1)')
    title('Btm POC')
    
    subplot(3,2,4)
    plot(y,LzooC_150m(lid,:),'k','Linewidth',1.5); 
    xlim([y(1) y(end)])
    %ylim([-5 2])
    xlabel('Year')
    ylabel('(g m^-^2)')
    title('LZ biomass')
    
    subplot(3,2,5)
    plot(y,Lzoo_loss_150m(lid,:),'k','Linewidth',1.5); 
    xlim([y(1) y(end)])
    %ylim([-5 2])
    xlabel('Year')
    ylabel('(g m^-^2)')
    title('LZ loss total')
    
    subplot(3,2,6)
    plot(y,Lzoo_quad_150m(lid,:),'k','Linewidth',1.5); 
    xlim([y(1) y(end)])
    %ylim([-5 2])
    xlabel('Year')
    ylabel('(g m^-^2 d^-1)')
    title('LZ loss quad')
    stamp(cfile)
    
    print('-dpng',[pp cfile '_ts_inputs_' loc '.png'])
    
    %  TIME SERIES ----------------------------------------------------
end

%% map of locations
locll = TEMP_150m(:,1);
locll(~isnan(locll)) = zeros;
locll(ID) = ones;
locll = reshape(locll,ni,nj);


clatlim=[-90 90];
clonlim=[-280 80];
load coastlines

%% locations not correct
figure(2)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,locll)
colormap('cool')
%caxis([0 35])
colorbar%('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
%text(0.2,1.65,'Tp','HorizontalAlignment','center','FontWeight','bold')
%h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
for i =1:16
    %textm(GRD.LAT(ids(i)),GRD.LON(ids(i)),abbrev{i},'Color','black','HorizontalAlignment','center');
    textm(TLAT(ID(i)),TLONG(ID(i)),abbrev{i},'Color','black','HorizontalAlignment','center');
end

