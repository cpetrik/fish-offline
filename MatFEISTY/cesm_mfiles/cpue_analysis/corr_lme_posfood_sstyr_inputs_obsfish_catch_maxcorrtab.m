% Use calc corr of catch with forcing, biomass, nu
% find most sig driver and lag
% min yrs as sat chl

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


mod = 'v15_obsfish';

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
load([spath 'LMEs_corr_catch_sstyrs_driver_lags.mat'])
stex = tanom;

load([spath 'LMEs_corr_catch_sstyrs_feisty_lags.mat'],'lid')

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
load([spath 'LMEs_corr_catch_sstyrs_obsfish_lags.mat'])
ftex = tanom;

%fish
CAtab(:,6:8,1:3) = AtabC;
CFtab(:,6:8,1:3) = FtabC;
CPtab(:,6:8,1:3) = PtabC;
CDtab(:,6:8,1:3) = DtabC;

PAtab(:,6:8,1:3) = AtabP;
PFtab(:,6:8,1:3) = FtabP;
PPtab(:,6:8,1:3) = PtabP;
PDtab(:,6:8,1:3) = DtabP;

clear AtabC AtabP FtabC FtabP PtabC PtabP DtabC DtabP tanom

%%
tanom = {'TP','TB','Det','ZmLoss','SST','Biom','Prod','Yield'};
cnam = {'corr','p','lag','idriver','driver'};

%Lags
yr = 0:3;  %reduce lags 0:4

% Drivers
tanom2=tanom';
tanom2(:,2)=tanom2(:,1);
tanom2(:,3)=tanom2(:,1);
tanom2(:,4)=tanom2(:,1);
tanom2(:,5)=tanom2(:,1);
tanom2(:,6)=tanom2(:,1);
tanom2(:,7)=tanom2(:,1);
tanom2(:,8)=tanom2(:,1);
tanom2(:,9)=tanom2(:,1);
tanom2(:,10)=tanom2(:,1);

[Ymat,Jmat] = meshgrid(yr,1:length(tanom));

LFtab = nan*ones(length(lid),4);
LPtab = nan*ones(length(lid),4);
LDtab = nan*ones(length(lid),4);
LAtab = nan*ones(length(lid),4);

LFt = cell(length(lid),1);
LPt = cell(length(lid),1);
LDt = cell(length(lid),1);
LAt = cell(length(lid),1);

%%
for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = num2str(i);

    AtabC = squeeze(CAtab(L,:,:));
    FtabC = squeeze(CFtab(L,:,:));
    PtabC = squeeze(CPtab(L,:,:));
    DtabC = squeeze(CDtab(L,:,:));

    AtabP = squeeze(PAtab(L,:,:));
    FtabP = squeeze(PFtab(L,:,:));
    PtabP = squeeze(PPtab(L,:,:));
    DtabP = squeeze(PDtab(L,:,:));

    %% force prey & fish corrs to be pos or zero (3,4,6,7,8)
    AtabC(3,AtabC(3,:)<0) = 0;
    AtabC(4,AtabC(4,:)<0) = 0;
    AtabC(6,AtabC(6,:)<0) = 0;
    AtabC(7,AtabC(7,:)<0) = 0;
    AtabC(8,AtabC(8,:)<0) = 0;
    
    FtabC(3,FtabC(3,:)<0) = 0;
    FtabC(4,FtabC(4,:)<0) = 0;
    FtabC(6,FtabC(6,:)<0) = 0;
    FtabC(7,FtabC(7,:)<0) = 0;
    FtabC(8,FtabC(8,:)<0) = 0;
    
    PtabC(3,PtabC(3,:)<0) = 0;
    PtabC(4,PtabC(4,:)<0) = 0;
    PtabC(6,PtabC(6,:)<0) = 0;
    PtabC(7,PtabC(7,:)<0) = 0;
    PtabC(8,PtabC(8,:)<0) = 0;
    
    DtabC(3,DtabC(3,:)<0) = 0;
    DtabC(4,DtabC(4,:)<0) = 0;
    DtabC(6,DtabC(6,:)<0) = 0;
    DtabC(7,DtabC(7,:)<0) = 0;
    DtabC(8,DtabC(8,:)<0) = 0;
    
    %%
    maxC = max(abs(AtabC(:)));
    pid = find(abs(AtabC(:))==maxC);
    LAtab(L,1) = AtabC(pid);
    LAtab(L,2) = AtabP(pid);
    LAtab(L,3) = Ymat(pid);
    LAtab(L,4) = Jmat(pid);
    LAt(L) = tanom2(pid);
    clear pid maxC

    maxC = max(abs(FtabC(:)));
    if(~isnan(maxC))
        pid = find(abs(FtabC(:))==maxC);
        LFtab(L,1) = FtabC(pid);
        LFtab(L,2) = FtabP(pid);
        LFtab(L,3) = Ymat(pid);
        LFtab(L,4) = Jmat(pid);
        LFt(L) = tanom2(pid);
    end
    clear pid maxC

    maxC = max(abs(PtabC(:)));
    if(~isnan(maxC))
        pid = find(abs(PtabC(:))==maxC);
        LPtab(L,1) = PtabC(pid);
        LPtab(L,2) = PtabP(pid);
        LPtab(L,3) = Ymat(pid);
        LPtab(L,4) = Jmat(pid);
        LPt(L) = tanom2(pid);
    end
    clear pid maxC

    maxC = max(abs(DtabC(:)));
    pid = find(abs(DtabC(:))==maxC);
    LDtab(L,1) = DtabC(pid);
    LDtab(L,2) = DtabP(pid);
    LDtab(L,3) = Ymat(pid);
    LDtab(L,4) = Jmat(pid);
    LDt(L) = tanom2(pid);
    clear pid maxC

    clear AtabC AtabP FtabC FtabP PtabC PtabP DtabC DtabP

end %LME

%%
lname = cellstr(num2str(lid));
% cnam, lname, tanom2
Atab1 = array2table(LAtab,"RowNames",lname);
Atab1(:,5) = LAt;
Atab1.Properties.VariableNames = cnam;

Ftab1 = array2table(LFtab,"RowNames",lname);
Ftab1(:,5) = LFt;
Ftab1.Properties.VariableNames = cnam;

Ptab1 = array2table(LPtab,"RowNames",lname);
Ptab1(:,5) = LPt;
Ptab1.Properties.VariableNames = cnam;

Dtab1 = array2table(LDtab,"RowNames",lname);
Dtab1(:,5) = LDt;
Dtab1.Properties.VariableNames = cnam;


%%
writetable(Atab1,[spath,'LMEs_corr_catch_sstyrs_inputs_obsfish_maxcorr_posfood_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ftab1,[spath,'LMEs_corr_catch_sstyrs_inputs_obsfish_maxcorr_posfood_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'LMEs_corr_catch_sstyrs_inputs_obsfish_maxcorr_posfood_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'LMEs_corr_catch_sstyrs_inputs_obsfish_maxcorr_posfood_D.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'LMEs_corr_catch_sstyrs_inputs_obsfish_maxcorr_posfood.mat'],...
    'LFtab','LPtab','LDtab','LAtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','lid');
