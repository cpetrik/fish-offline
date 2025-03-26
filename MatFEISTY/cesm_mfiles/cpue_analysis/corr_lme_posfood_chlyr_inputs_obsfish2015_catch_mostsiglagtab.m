% Use calc corr of catch with forcing, biomass, nu
% find most sig lag of each driver in each LME
% min yrs as sat chl
% obs fishing effort

clear
close all

%% FOSI input forcing
%cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
cpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/cpue2015/';

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

% FEISTY outputs
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

% fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/cpue2015/';
% spath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/cpue2015/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regress_cpue/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',...
    cfile,'/corrs_cpue'];


mod = 'v15_obsfish2015';

% Fishing data
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

%% All corrs
CFtab = nan*ones(63,8,4);
PFtab = nan*ones(63,8,4);
CPtab = CFtab;
PPtab = CFtab;
CDtab = CFtab;
PDtab = CFtab;
CAtab = CFtab;
PAtab = CFtab;

%% sat & inputs
load([spath 'LMEs_corr_catch_chlyrs15_driver_lags.mat'])
stex = tanom;

load([spath 'LMEs_corr_catch_chlyrs15_feisty_lags.mat'],'lid')

%inputss & sat
CAtab(:,1:6,:) = AtabC(lid,:,:);
CFtab(:,1:6,:) = FtabC(lid,:,:);
CPtab(:,1:6,:) = PtabC(lid,:,:);
CDtab(:,1:6,:) = DtabC(lid,:,:);

PAtab(:,1:6,:) = AtabP(lid,:,:);
PFtab(:,1:6,:) = FtabP(lid,:,:);
PPtab(:,1:6,:) = PtabP(lid,:,:);
PDtab(:,1:6,:) = DtabP(lid,:,:);

clear AtabC AtabP FtabC FtabP PtabC PtabP DtabC DtabP tanom

%%
load([spath 'LMEs_corr_catch_chlyrs15_obsfish2015_lags.mat'])
ftex = tanom;

%sat
CAtab(:,7:9,1:3) = AtabC(:,1:3,:);
CFtab(:,7:9,1:3) = FtabC(:,1:3,:);
CPtab(:,7:9,1:3) = PtabC(:,1:3,:);
CDtab(:,7:9,1:3) = DtabC(:,1:3,:);

PAtab(:,7:9,1:3) = AtabP(:,1:3,:);
PFtab(:,7:9,1:3) = FtabP(:,1:3,:);
PPtab(:,7:9,1:3) = PtabP(:,1:3,:);
PDtab(:,7:9,1:3) = DtabP(:,1:3,:);

clear AtabC AtabP FtabC FtabP PtabC PtabP DtabC DtabP tanom

%%
tanom = {'TP','TB','Det','ZmLoss','SST','Chl','Biom','Prod','Yield'};
cnam = {'corr','p','lag'};

%Lags
yr = 0:3;  %reduce lags 0:4

LFtab = nan*ones(length(lid),length(tanom),3);
LPtab = nan*ones(length(lid),length(tanom),3);
LDtab = nan*ones(length(lid),length(tanom),3);
LAtab = nan*ones(length(lid),length(tanom),3);

%%
for j=1:length(tanom)
    driver = tanom{j};

    for L = 1:length(lid)

        %LME
        i = lid(L);
        ilme = num2str(i);

        AtabC = squeeze(CAtab(L,j,:));
        FtabC = squeeze(CFtab(L,j,:));
        PtabC = squeeze(CPtab(L,j,:));
        DtabC = squeeze(CDtab(L,j,:));

        AtabP = squeeze(PAtab(L,j,:));
        FtabP = squeeze(PFtab(L,j,:));
        PtabP = squeeze(PPtab(L,j,:));
        DtabP = squeeze(PDtab(L,j,:));

        %% force prey & fish corrs to be pos or zero (3,4,6,7,8)
        if j~=1
            if j~=2
                if j~=5
                    AtabP(AtabC<0) = 1;
                    FtabP(FtabC<0) = 1;
                    PtabP(PtabC<0) = 1;
                    DtabP(DtabC<0) = 1;

                    AtabC(AtabC<0) = 0;
                    FtabC(FtabC<0) = 0;
                    PtabC(PtabC<0) = 0;
                    DtabC(DtabC<0) = 0;
                end
            end
        end

        %%
        maxC = max(abs(AtabC(:)));
        pid = find(abs(AtabC(:))==maxC);
        if(length(pid)>1)
            pid = pid(1);
        end
        LAtab(L,j,1) = AtabC(pid);
        LAtab(L,j,2) = AtabP(pid);
        LAtab(L,j,3) = yr(pid);
        clear pid maxC

        maxC = max(abs(FtabC(:)));
        if(~isnan(maxC))
            pid = find(abs(FtabC(:))==maxC);
            if(length(pid)>1)
                pid = pid(1);
            end
            LFtab(L,j,1) = FtabC(pid);
            LFtab(L,j,2) = FtabP(pid);
            LFtab(L,j,3) = yr(pid);
        end
        clear pid maxC

        maxC = max(abs(PtabC(:)));
        if(~isnan(maxC))
            pid = find(abs(PtabC(:))==maxC);
            if(length(pid)>1)
                pid = pid(1);
            end
            LPtab(L,j,1) = PtabC(pid);
            LPtab(L,j,2) = PtabP(pid);
            LPtab(L,j,3) = yr(pid);
        end
        clear pid maxC

        maxC = max(abs(DtabC(:)));
        pid = find(abs(DtabC(:))==maxC);
        if(length(pid)>1)
            pid = pid(1);
        end
        LDtab(L,j,1) = DtabC(pid);
        LDtab(L,j,2) = DtabP(pid);
        LDtab(L,j,3) = yr(pid);
        clear pid maxC

        clear AtabC AtabP FtabC FtabP PtabC PtabP DtabC DtabP

    end %LME

end %driver

save([spath,'LMEs_corr_catch_chlyrs15_inputs_obsfish2015_mostsiglag_posfood.mat'],...
    'LFtab','LPtab','LDtab','LAtab','lid','tanom','cnam','yr');

