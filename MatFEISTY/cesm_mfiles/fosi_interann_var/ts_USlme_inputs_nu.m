% Time series of forcing and fish
% nu instead of biomass

clear
close all


% LMEs
lid = [54,1:2,65,10,3,5:7]; %ADD 65 = Aleutian Islands
lname = {'CHK','EBS','GAK','AI','HI','CCE','GMX','SE','NE'};

%% FOSI input forcing

%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adet','adety','atb','atp','azlos','azlosy','azoo');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

% put into a matrix & use annual nu
manom(:,:,1) = atp;
manom(:,:,2) = atb;
manom(:,:,3) = azoo;
manom(:,:,4) = azlosy;
manom(:,:,5) = adety;

tanom = {'TP','TB','Zmeso','ZmLossY','DetY'};

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/'];
cpath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/' cfile '/corrs/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

% Anoms with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms.mat'],...
    'aa','ad','af','ap')

%%
cm10=[0.5 0.5 0.5; ... %med grey
    0 0 0; ...      %black
    0 0 0.75;...    %b
    1 0 0;...       %r
    0 0.5 0.75;...   %med blue
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    1 0 1;...     %m
    0 0.7 0;...   %g
    ];

set(groot,'defaultAxesColorOrder',cm10);

%%
figure(1)
tiledlayout(3,3, 'TileSpacing', 'compact')
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = lname{L};

    nexttile

    for j = 1:length(tanom)

        %input forcing
        driver = tanom{j};
        %              LME  time    driver                LME   time   driver
        sclim = ((manom(i,:,j))') ./ (2*std((manom(i,:,j))));

        plot(1:68,sclim); hold on;
        ylim([-1.5 1.5])
        title(ilme)

    end % driver
    %%

end %LME
lg = legend(nexttile(6),tanom);
lg.Location = 'eastoutside';
%lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_USLME_forcing_anoms_div2SD.png'])

%% not divided
figure(2)
tiledlayout(3,3, 'TileSpacing', 'compact')
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = lname{L};

    nexttile

    for j = 1:2

        %input forcing
        driver = tanom{j};
        %              LME  time    driver
        sclim = ((manom(i,:,j))') ;

        plot(1:68,sclim); hold on;
        %ylim([-1.5 1.5])
        title(ilme)

    end % driver
    %%

end %LME
lg = legend(nexttile(6),tanom{1:2});
lg.Location = 'eastoutside';
%lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_USLME_forcing_temp_anoms.png'])

figure(3)
tiledlayout(3,3, 'TileSpacing', 'compact')
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = lname{L};

    nexttile

    for j = 3

        %input forcing
        driver = tanom{j};
        %              LME  time    driver
        sclim = ((manom(i,:,j))') ;

        plot(1:68,sclim); hold on;
        %ylim([-1.5 1.5])
        title(ilme)

    end % driver
    %%

end %LME
print('-dpng',[ppath 'ts_USLME_forcing_zmeso_anom.png'])

figure(4)
tiledlayout(3,3, 'TileSpacing', 'compact')
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = lname{L};

    nexttile

    for j = 4:length(tanom)

        %input forcing
        driver = tanom{j};
        %              LME  time    driver
        sclim = ((manom(i,:,j))') ;

        plot(1:68,sclim); hold on;
        %ylim([-1.5 1.5])
        title(ilme)

    end % driver
    %%

end %LME
lg = legend(nexttile(6),tanom{4:5});
lg.Location = 'eastoutside';
%lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_USLME_forcing_ZmLossDet_anoms.png'])

%% Fish
figure(5)
tiledlayout(3,3, 'TileSpacing', 'compact')
for L = 1:length(lid)
    %LME
    i = lid(L);
    ilme = lname{L};

    nexttile

    plot(1:68,af(i,:)); hold on;
    plot(1:68,ap(i,:)); hold on;
    plot(1:68,ad(i,:)); hold on;
    plot(1:68,aa(i,:)); hold on;
    
    %ylim([-1.5 1.5])
    title(ilme)
end %LME
lg = legend(nexttile(6),{'F','P','D','All'});
lg.Location = 'eastoutside';
%lg.Orientation = 'horizontal';
print('-dpng',[cpath 'ts_USLME_fntypes_nu_anoms.png'])

%% Fish div by 2SD
figure(6)
tiledlayout(3,3, 'TileSpacing', 'compact')
for L = 1:length(lid)
    %LME
    i = lid(L);
    ilme = lname{L};

    nexttile
    
    plot(1:68,(af(i,:)./(2*std(af(i,:))))); hold on;
    plot(1:68,(ap(i,:)./(2*std(ap(i,:))))); hold on;
    plot(1:68,(ad(i,:)./(2*std(ad(i,:))))); hold on;
    plot(1:68,(aa(i,:)./(2*std(aa(i,:))))); hold on;
    ylim([-1.5 1.5])
    title(ilme)
end %LME
lg = legend(nexttile(6),{'F','P','D','All'});
lg.Location = 'eastoutside';
%lg.Orientation = 'horizontal';
print('-dpng',[cpath 'ts_USLME_fntypes_nu_anoms_div2SD.png'])


