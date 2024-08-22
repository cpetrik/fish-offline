% Use calc corr of cpue with forcing, biomass, nu
% find most sig lag of each driver in each LME
% min yrs as sat sst
% const fishing effort

clear
close all

%% FOSI input forcing
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

% FEISTY outputs
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regress_cpue/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',...
    cfile,'/corrs_cpue'];


mod = 'v15_All_fish03';

% Fishing data
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

%% All corrs
CFtab = nan*ones(63,7,4);
PFtab = nan*ones(63,7,4);
CPtab = CFtab;
PPtab = CFtab;
CDtab = CFtab;
PDtab = CFtab;
CAtab = CFtab;
PAtab = CFtab;

%% sat & inputs
load([spath 'LMEs_corr_cpue_sstyrs15_driver_lags.mat'])
stex = tanom;

load([spath 'LMEs_corr_cpue_sstyrs15_feisty_lags.mat'],'lid')

%inputss & sat
CAtab(:,1:5,:) = AtabC(lid,:,:);
CFtab(:,1:5,:) = FtabC(lid,:,:);
CPtab(:,1:5,:) = PtabC(lid,:,:);
CDtab(:,1:5,:) = DtabC(lid,:,:);

PAtab(:,1:5,:) = AtabP(lid,:,:);
PFtab(:,1:5,:) = FtabP(lid,:,:);
PPtab(:,1:5,:) = PtabP(lid,:,:);
PDtab(:,1:5,:) = DtabP(lid,:,:);

clear AtabC AtabP FtabC FtabP PtabC PtabP DtabC DtabP tanom

%%
load([spath 'LMEs_corr_cpue_sstyrs15_feisty_lags.mat'])
ftex = tanom;

%sat
CAtab(:,6:7,1:3) = AtabC(:,1:2,:);
CFtab(:,6:7,1:3) = FtabC(:,1:2,:);
CPtab(:,6:7,1:3) = PtabC(:,1:2,:);
CDtab(:,6:7,1:3) = DtabC(:,1:2,:);

PAtab(:,6:7,1:3) = AtabP(:,1:2,:);
PFtab(:,6:7,1:3) = FtabP(:,1:2,:);
PPtab(:,6:7,1:3) = PtabP(:,1:2,:);
PDtab(:,6:7,1:3) = DtabP(:,1:2,:);

clear AtabC AtabP FtabC FtabP PtabC PtabP DtabC DtabP tanom

%%
tanom = {'TP','TB','Det','ZmLoss','SST','Biom','Prod'};
cnam = {'corr','p','lag'};

%Lags
yr = 0:3;  %reduce lags 0:4

LFmat = nan*ones(length(lid),length(tanom),3);
LPmat = nan*ones(length(lid),length(tanom),3);
LDmat = nan*ones(length(lid),length(tanom),3);
LAmat = nan*ones(length(lid),length(tanom),3);

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
        LAmat(L,j,1) = AtabC(pid);
        LAmat(L,j,2) = AtabP(pid);
        LAmat(L,j,3) = yr(pid);
        clear pid maxC

        maxC = max(abs(FtabC(:)));
        if(~isnan(maxC))
            pid = find(abs(FtabC(:))==maxC);
            if(length(pid)>1)
                pid = pid(1);
            end
            LFmat(L,j,1) = FtabC(pid);
            LFmat(L,j,2) = FtabP(pid);
            LFmat(L,j,3) = yr(pid);
        end
        clear pid maxC

        maxC = max(abs(PtabC(:)));
        if(~isnan(maxC))
            pid = find(abs(PtabC(:))==maxC);
            if(length(pid)>1)
                pid = pid(1);
            end
            LPmat(L,j,1) = PtabC(pid);
            LPmat(L,j,2) = PtabP(pid);
            LPmat(L,j,3) = yr(pid);
        end
        clear pid maxC

        maxC = max(abs(DtabC(:)));
        pid = find(abs(DtabC(:))==maxC);
        if(length(pid)>1)
            pid = pid(1);
        end
        LDmat(L,j,1) = DtabC(pid);
        LDmat(L,j,2) = DtabP(pid);
        LDmat(L,j,3) = yr(pid);
        clear pid maxC

        clear AtabC AtabP FtabC FtabP PtabC PtabP DtabC DtabP

    end %LME

end %driver

save([spath,'LMEs_corr_cpue_sstyrs15_inputs_feisty_mostsiglag_posfood.mat'],...
    'LFmat','LPmat','LDmat','LAmat','lid','tanom','cnam','yr');

