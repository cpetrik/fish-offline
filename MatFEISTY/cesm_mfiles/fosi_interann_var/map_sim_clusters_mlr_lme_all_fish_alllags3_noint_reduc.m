% Map clusters of driver-fish mult linear regressions
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/' cfile '/corrs/'];

mod = 'v15_All_fish03_';

% Ecosystem type
load([fpath 'LME_fosi_fished_',mod,cfile '.mat'],'etype','etex','Ebiom');

%spath ='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';
load([spath,'LME_biom_nu_cpue_cme_A_mlr_coeffs_reduc_alllags3_R2_cluster.mat'])

%% vector of clusters to use old code
ClusterB = Abiom(:,7);
ClusterP = Anu(:,7);
ClusterC = Acpue(:,7);
ClusterM = Acme(:,7);

%% Cluster descrips
% alltex = {'-Det',...    % = 1
    % '+Det',...          % = 2
    % '+Det, +ZL',...     % = 3
    % '+Det, -ZL',...     % = 4
    % '-Det, +ZL',...     % = 5
    % '+ZL',...           % = 6
    % '+TB',...           % = 7
    % '+TB, -TP',...      % = 8
    % '+TB, -TP, +ZL',... % = 9
    % '-TB',...           % = 10
    % '-TB, +TP',...      % = 11
    % '+TP'};             % = 12

btexD = {'+Det',... % = 1
    '',...      % = 2
    '',...      % = 3
    '+Det',...% = 4
    '-Det',...     % = 5
    '+Det',...     % = 6
    '-Det',...     % = 7
    '+Det',...     % = 8
    '+Det'};            % = 9

btexZ = {'',... % = 1
    '',...      % = 2
    '+ZL',...      % = 3
    '+ZL',...% = 4
    '+ZL',...     % = 5
    '',...     % = 6
    '+ZL',...     % = 7
    '',...     % = 8
    ''};            % = 9

btexP = {'',... % = 1
    '-TP',...      % = 2
    '',...      % = 3
    '',...% = 4
    '',...     % = 5
    '',...     % = 6
    '',...     % = 7
    '',...     % = 8
    ''};            % = 9

btexB = {'',... % = 1
    '+TB',...      % = 2
    '',...      % = 3
    '',...% = 4
    '',...     % = 5
    '',...     % = 6
    '',...     % = 7
    '-TB',...     % = 8
    ''};            % = 9

ptexD = {'',...       % = 1
    '-Det',...     % = 2
    '',...           % = 3
    '+Det',...          % = 4
    '',...      % = 5
    ''};   % = 6

ptexZ = {'',...       % = 1
    '+ZL',...     % = 2
    '+ZL',...           % = 3
    '',...          % = 4
    '',...      % = 5
    '+ZL'};   % = 6

ptexP = {'+TP',...       % = 1
    '',...     % = 2
    '',...           % = 3
    '',...          % = 4
    '+TP',...      % = 5
    '-TP'};   % = 6

ptexB = {'',...       % = 1
    '',...     % = 2
    '',...           % = 3
    '',...          % = 4
    '-TB',...      % = 5
    '+TB'};   % = 6
               
%cpue
ctexD = {'-Det',...      % = 1
    '-Det',...     % = 2
    '',...           % = 3
    '+Det',...     % = 4
    '',...           % = 5
    ''};        % = 6

ctexZ = {'',...      % = 1
    '+ZL',...     % = 2
    '',...           % = 3
    '-ZL',...     % = 4
    '',...           % = 5
    ''};        % = 6

ctexP = {'',...      % = 1
    '',...     % = 2
    '',...           % = 3
    '',...     % = 4
    '',...           % = 5
    '+TP'};        % = 6

ctexB = {'',...      % = 1
    '',...     % = 2
    '-TB',...           % = 3
    '',...     % = 4
    '+TB',...           % = 5
    '-TB'};        % = 6
       
%cme
mtexD = {'-Det',...       % = 1
    '',...       % = 2
    '+Det',...      % = 3
    '',...       % = 4
    '',...       % = 5
    '+Det',...      % = 6
    '+Det',...      % = 7
    ''};         % = 8

mtexZ = {'',...       % = 1
    '',...       % = 2
    '-ZL',...      % = 3
    '',...       % = 4
    '',...       % = 5
    '+ZL',...      % = 6
    '-ZL',...      % = 7
    ''};         % = 8

mtexP = {'',...       % = 1
    '+TP',...       % = 2
    '',...      % = 3
    '-TP',...       % = 4
    '+TP',...       % = 5
    '',...      % = 6
    '',...      % = 7
    '-TP'};         % = 8

mtexB = {'',...       % = 1
    '-TB',...       % = 2
    '',...      % = 3
    '+TB',...       % = 4
    '-TB',...       % = 5
    '',...      % = 6
    '',...      % = 7
    '+TB'};         % = 8
   
%% Add a text descript vector that matches to clusters
% 
%biomass
CatBP = cell(length(ClusterB),1);
CatBP(ClusterB==1,1) = btexP(1);
CatBP(ClusterB==2,1) = btexP(2);
CatBP(ClusterB==3,1) = btexP(3);
CatBP(ClusterB==4,1) = btexP(4);
CatBP(ClusterB==5,1) = btexP(5);
CatBP(ClusterB==6,1) = btexP(6);
CatBP(ClusterB==7,1) = btexP(7);
CatBP(ClusterB==8,1) = btexP(8);
CatBP(ClusterB==9,1) = btexP(9);
CatBP(isnan(ClusterB)) = {''};

CatBB = cell(length(ClusterB),1);
CatBB(ClusterB==1,1) = btexB(1);
CatBB(ClusterB==2,1) = btexB(2);
CatBB(ClusterB==3,1) = btexB(3);
CatBB(ClusterB==4,1) = btexB(4);
CatBB(ClusterB==5,1) = btexB(5);
CatBB(ClusterB==6,1) = btexB(6);
CatBB(ClusterB==7,1) = btexB(7);
CatBB(ClusterB==8,1) = btexB(8);
CatBB(ClusterB==9,1) = btexB(9);
CatBB(isnan(ClusterB)) = {''};

CatBZ = cell(length(ClusterB),1);
CatBZ(ClusterB==1,1) = btexZ(1);
CatBZ(ClusterB==2,1) = btexZ(2);
CatBZ(ClusterB==3,1) = btexZ(3);
CatBZ(ClusterB==4,1) = btexZ(4);
CatBZ(ClusterB==5,1) = btexZ(5);
CatBZ(ClusterB==6,1) = btexZ(6);
CatBZ(ClusterB==7,1) = btexZ(7);
CatBZ(ClusterB==8,1) = btexZ(8);
CatBZ(ClusterB==9,1) = btexZ(9);
CatBZ(isnan(ClusterB)) = {''};

CatBD = cell(length(ClusterB),1);
CatBD(ClusterB==1,1) = btexD(1);
CatBD(ClusterB==2,1) = btexD(2);
CatBD(ClusterB==3,1) = btexD(3);
CatBD(ClusterB==4,1) = btexD(4);
CatBD(ClusterB==5,1) = btexD(5);
CatBD(ClusterB==6,1) = btexD(6);
CatBD(ClusterB==7,1) = btexD(7);
CatBD(ClusterB==8,1) = btexD(8);
CatBD(ClusterB==9,1) = btexD(9);
CatBD(isnan(ClusterB)) = {''};

%nu, prod
CatPP = cell(length(ClusterP),1);
CatPP(ClusterP==1,1) = ptexP(1);
CatPP(ClusterP==2,1) = ptexP(2);
CatPP(ClusterP==3,1) = ptexP(3);
CatPP(ClusterP==4,1) = ptexP(4);
CatPP(ClusterP==5,1) = ptexP(5);
CatPP(ClusterP==6,1) = ptexP(6);
CatPP(isnan(ClusterP)) = {''};

CatPB = cell(length(ClusterP),1);
CatPB(ClusterP==1,1) = ptexB(1);
CatPB(ClusterP==2,1) = ptexB(2);
CatPB(ClusterP==3,1) = ptexB(3);
CatPB(ClusterP==4,1) = ptexB(4);
CatPB(ClusterP==5,1) = ptexB(5);
CatPB(ClusterP==6,1) = ptexB(6);
CatPB(isnan(ClusterP)) = {''};

CatPZ = cell(length(ClusterP),1);
CatPZ(ClusterP==1,1) = ptexZ(1);
CatPZ(ClusterP==2,1) = ptexZ(2);
CatPZ(ClusterP==3,1) = ptexZ(3);
CatPZ(ClusterP==4,1) = ptexZ(4);
CatPZ(ClusterP==5,1) = ptexZ(5);
CatPZ(ClusterP==6,1) = ptexZ(6);
CatPZ(isnan(ClusterP)) = {''};

CatPD = cell(length(ClusterP),1);
CatPD(ClusterP==1,1) = ptexD(1);
CatPD(ClusterP==2,1) = ptexD(2);
CatPD(ClusterP==3,1) = ptexD(3);
CatPD(ClusterP==4,1) = ptexD(4);
CatPD(ClusterP==5,1) = ptexD(5);
CatPD(ClusterP==6,1) = ptexD(6);
CatPD(isnan(ClusterP)) = {''};

%cpue
CatCP = cell(length(ClusterC),1);
CatCP(ClusterC==1,1) = ctexP(1);
CatCP(ClusterC==2,1) = ctexP(2);
CatCP(ClusterC==3,1) = ctexP(3);
CatCP(ClusterC==4,1) = ctexP(4);
CatCP(ClusterC==5,1) = ctexP(5);
CatCP(ClusterC==6,1) = ctexP(6);
CatCP(isnan(ClusterC)) = {''};

CatCB = cell(length(ClusterC),1);
CatCB(ClusterC==1,1) = ctexB(1);
CatCB(ClusterC==2,1) = ctexB(2);
CatCB(ClusterC==3,1) = ctexB(3);
CatCB(ClusterC==4,1) = ctexB(4);
CatCB(ClusterC==5,1) = ctexB(5);
CatCB(ClusterC==6,1) = ctexB(6);
CatCB(isnan(ClusterC)) = {''};

CatCZ = cell(length(ClusterC),1);
CatCZ(ClusterC==1,1) = ctexZ(1);
CatCZ(ClusterC==2,1) = ctexZ(2);
CatCZ(ClusterC==3,1) = ctexZ(3);
CatCZ(ClusterC==4,1) = ctexZ(4);
CatCZ(ClusterC==5,1) = ctexZ(5);
CatCZ(ClusterC==6,1) = ctexZ(6);
CatCZ(isnan(ClusterC)) = {''};

CatCD = cell(length(ClusterC),1);
CatCD(ClusterC==1,1) = ctexD(1);
CatCD(ClusterC==2,1) = ctexD(2);
CatCD(ClusterC==3,1) = ctexD(3);
CatCD(ClusterC==4,1) = ctexD(4);
CatCD(ClusterC==5,1) = ctexD(5);
CatCD(ClusterC==6,1) = ctexD(6);
CatCD(isnan(ClusterC)) = {''};

%cme
CatMP = cell(length(ClusterM),1);
CatMP(ClusterM==1,1) = mtexP(1);
CatMP(ClusterM==2,1) = mtexP(2);
CatMP(ClusterM==3,1) = mtexP(3);
CatMP(ClusterM==4,1) = mtexP(4);
CatMP(ClusterM==5,1) = mtexP(5);
CatMP(ClusterM==6,1) = mtexP(6);
CatMP(ClusterM==7,1) = mtexP(7);
CatMP(ClusterM==8,1) = mtexP(8);
CatMP(isnan(ClusterM)) = {''};

CatMB = cell(length(ClusterM),1);
CatMB(ClusterM==1,1) = mtexB(1);
CatMB(ClusterM==2,1) = mtexB(2);
CatMB(ClusterM==3,1) = mtexB(3);
CatMB(ClusterM==4,1) = mtexB(4);
CatMB(ClusterM==5,1) = mtexB(5);
CatMB(ClusterM==6,1) = mtexB(6);
CatMB(ClusterM==7,1) = mtexB(7);
CatMB(ClusterM==8,1) = mtexB(8);
CatMB(isnan(ClusterM)) = {''};

CatMZ = cell(length(ClusterM),1);
CatMZ(ClusterM==1,1) = mtexZ(1);
CatMZ(ClusterM==2,1) = mtexZ(2);
CatMZ(ClusterM==3,1) = mtexZ(3);
CatMZ(ClusterM==4,1) = mtexZ(4);
CatMZ(ClusterM==5,1) = mtexZ(5);
CatMZ(ClusterM==6,1) = mtexZ(6);
CatMZ(ClusterM==7,1) = mtexZ(7);
CatMZ(ClusterM==8,1) = mtexZ(8);
CatMZ(isnan(ClusterM)) = {''};

CatMD = cell(length(ClusterM),1);
CatMD(ClusterM==1,1) = mtexD(1);
CatMD(ClusterM==2,1) = mtexD(2);
CatMD(ClusterM==3,1) = mtexD(3);
CatMD(ClusterM==4,1) = mtexD(4);
CatMD(ClusterM==5,1) = mtexD(5);
CatMD(ClusterM==6,1) = mtexD(6);
CatMD(ClusterM==7,1) = mtexD(7);
CatMD(ClusterM==8,1) = mtexD(8);
CatMD(isnan(ClusterM)) = {''};

%%
CatBP = string(CatBP);
CatPP = string(CatPP);
CatCP = string(CatCP);
CatMP = string(CatMP);

CatBB = string(CatBB);
CatPB = string(CatPB);
CatCB = string(CatCB);
CatMB = string(CatMB);

CatBZ = string(CatBZ);
CatPZ = string(CatPZ);
CatCZ = string(CatCZ);
CatMZ = string(CatMZ);

CatBD = string(CatBD);
CatPD = string(CatPD);
CatCD = string(CatCD);
CatMD = string(CatMD);

% Put text in table
% cCatM = char(CatM);
% cCatP = char(CatP);
% cCatC = char(CatC);
% cCatB = char(CatB);
% Atab = table(cCatB,cCatP,cCatC,cCatM,'VariableNames',...
%     {'Biomass','Production','CPUE','Catch-Effort'});

%writetable(Atab,[spath 'LMEs_driver_mlr_AllFish_cluster_v3_fntypes.csv'])


%% Create one master colormap for all categories
alltex = {'+Det',...    % = 1
    '-Det',...          % = 2
    '-ZL',...           % = 3
    '+ZL',...           % = 4
    '+TB',...           % = 5
    '-TB',...           % = 6
    '-TP',...           % = 7
    '+TP'};             % = 8

%Biomass
ClusterB(:,2) = nan; %TP
ClusterB((CatBP==alltex(1)),2) = 1;
ClusterB((CatBP==alltex(2)),2) = 2;
ClusterB((CatBP==alltex(3)),2) = 3;
ClusterB((CatBP==alltex(4)),2) = 4;
ClusterB((CatBP==alltex(5)),2) = 5;
ClusterB((CatBP==alltex(6)),2) = 6;
ClusterB((CatBP==alltex(7)),2) = 7;
ClusterB((CatBP==alltex(8)),2) = 8;

ClusterB(:,3) = nan; %TB
ClusterB((CatBB==alltex(1)),3) = 1;
ClusterB((CatBB==alltex(2)),3) = 2;
ClusterB((CatBB==alltex(3)),3) = 3;
ClusterB((CatBB==alltex(4)),3) = 4;
ClusterB((CatBB==alltex(5)),3) = 5;
ClusterB((CatBB==alltex(6)),3) = 6;
ClusterB((CatBB==alltex(7)),3) = 7;
ClusterB((CatBB==alltex(8)),3) = 8;

ClusterB(:,4) = nan; %ZmL
ClusterB((CatBZ==alltex(1)),4) = 1;
ClusterB((CatBZ==alltex(2)),4) = 2;
ClusterB((CatBZ==alltex(3)),4) = 3;
ClusterB((CatBZ==alltex(4)),4) = 4;
ClusterB((CatBZ==alltex(5)),4) = 5;
ClusterB((CatBZ==alltex(6)),4) = 6;
ClusterB((CatBZ==alltex(7)),4) = 7;
ClusterB((CatBZ==alltex(8)),4) = 8;

ClusterB(:,5) = nan; %Det
ClusterB((CatBD==alltex(1)),5) = 1;
ClusterB((CatBD==alltex(2)),5) = 2;
ClusterB((CatBD==alltex(3)),5) = 3;
ClusterB((CatBD==alltex(4)),5) = 4;
ClusterB((CatBD==alltex(5)),5) = 5;
ClusterB((CatBD==alltex(6)),5) = 6;
ClusterB((CatBD==alltex(7)),5) = 7;
ClusterB((CatBD==alltex(8)),5) = 8;

% Prod
ClusterP(:,2) = nan; %TP
ClusterP((CatPP==alltex(1)),2) = 1;
ClusterP((CatPP==alltex(2)),2) = 2;
ClusterP((CatPP==alltex(3)),2) = 3;
ClusterP((CatPP==alltex(4)),2) = 4;
ClusterP((CatPP==alltex(5)),2) = 5;
ClusterP((CatPP==alltex(6)),2) = 6;
ClusterP((CatPP==alltex(7)),2) = 7;
ClusterP((CatPP==alltex(8)),2) = 8;

ClusterP(:,3) = nan; %TB
ClusterP((CatPB==alltex(1)),3) = 1;
ClusterP((CatPB==alltex(2)),3) = 2;
ClusterP((CatPB==alltex(3)),3) = 3;
ClusterP((CatPB==alltex(4)),3) = 4;
ClusterP((CatPB==alltex(5)),3) = 5;
ClusterP((CatPB==alltex(6)),3) = 6;
ClusterP((CatPB==alltex(7)),3) = 7;
ClusterP((CatPB==alltex(8)),3) = 8;

ClusterP(:,4) = nan; %ZmL
ClusterP((CatPZ==alltex(1)),4) = 1;
ClusterP((CatPZ==alltex(2)),4) = 2;
ClusterP((CatPZ==alltex(3)),4) = 3;
ClusterP((CatPZ==alltex(4)),4) = 4;
ClusterP((CatPZ==alltex(5)),4) = 5;
ClusterP((CatPZ==alltex(6)),4) = 6;
ClusterP((CatPZ==alltex(7)),4) = 7;
ClusterP((CatPZ==alltex(8)),4) = 8;

ClusterP(:,5) = nan; %Det
ClusterP((CatPD==alltex(1)),5) = 1;
ClusterP((CatPD==alltex(2)),5) = 2;
ClusterP((CatPD==alltex(3)),5) = 3;
ClusterP((CatPD==alltex(4)),5) = 4;
ClusterP((CatPD==alltex(5)),5) = 5;
ClusterP((CatPD==alltex(6)),5) = 6;
ClusterP((CatPD==alltex(7)),5) = 7;
ClusterP((CatPD==alltex(8)),5) = 8;

%CPUE
ClusterC(:,2) = nan; %TP
ClusterC((CatCP==alltex(1)),2) = 1;
ClusterC((CatCP==alltex(2)),2) = 2;
ClusterC((CatCP==alltex(3)),2) = 3;
ClusterC((CatCP==alltex(4)),2) = 4;
ClusterC((CatCP==alltex(5)),2) = 5;
ClusterC((CatCP==alltex(6)),2) = 6;
ClusterC((CatCP==alltex(7)),2) = 7;
ClusterC((CatCP==alltex(8)),2) = 8;

ClusterC(:,3) = nan; %TB
ClusterC((CatCB==alltex(1)),3) = 1;
ClusterC((CatCB==alltex(2)),3) = 2;
ClusterC((CatCB==alltex(3)),3) = 3;
ClusterC((CatCB==alltex(4)),3) = 4;
ClusterC((CatCB==alltex(5)),3) = 5;
ClusterC((CatCB==alltex(6)),3) = 6;
ClusterC((CatCB==alltex(7)),3) = 7;
ClusterC((CatCB==alltex(8)),3) = 8;

ClusterC(:,4) = nan; %ZmL
ClusterC((CatCZ==alltex(1)),4) = 1;
ClusterC((CatCZ==alltex(2)),4) = 2;
ClusterC((CatCZ==alltex(3)),4) = 3;
ClusterC((CatCZ==alltex(4)),4) = 4;
ClusterC((CatCZ==alltex(5)),4) = 5;
ClusterC((CatCZ==alltex(6)),4) = 6;
ClusterC((CatCZ==alltex(7)),4) = 7;
ClusterC((CatCZ==alltex(8)),4) = 8;

ClusterC(:,5) = nan; %Det
ClusterC((CatCD==alltex(1)),5) = 1;
ClusterC((CatCD==alltex(2)),5) = 2;
ClusterC((CatCD==alltex(3)),5) = 3;
ClusterC((CatCD==alltex(4)),5) = 4;
ClusterC((CatCD==alltex(5)),5) = 5;
ClusterC((CatCD==alltex(6)),5) = 6;
ClusterC((CatCD==alltex(7)),5) = 7;
ClusterC((CatCD==alltex(8)),5) = 8;

%CME
ClusterM(:,2) = nan; %TP
ClusterM((CatMP==alltex(1)),2) = 1;
ClusterM((CatMP==alltex(2)),2) = 2;
ClusterM((CatMP==alltex(3)),2) = 3;
ClusterM((CatMP==alltex(4)),2) = 4;
ClusterM((CatMP==alltex(5)),2) = 5;
ClusterM((CatMP==alltex(6)),2) = 6;
ClusterM((CatMP==alltex(7)),2) = 7;
ClusterM((CatMP==alltex(8)),2) = 8;

ClusterM(:,3) = nan; %TB
ClusterM((CatMB==alltex(1)),3) = 1;
ClusterM((CatMB==alltex(2)),3) = 2;
ClusterM((CatMB==alltex(3)),3) = 3;
ClusterM((CatMB==alltex(4)),3) = 4;
ClusterM((CatMB==alltex(5)),3) = 5;
ClusterM((CatMB==alltex(6)),3) = 6;
ClusterM((CatMB==alltex(7)),3) = 7;
ClusterM((CatMB==alltex(8)),3) = 8;

ClusterM(:,4) = nan; %ZmL
ClusterM((CatMZ==alltex(1)),4) = 1;
ClusterM((CatMZ==alltex(2)),4) = 2;
ClusterM((CatMZ==alltex(3)),4) = 3;
ClusterM((CatMZ==alltex(4)),4) = 4;
ClusterM((CatMZ==alltex(5)),4) = 5;
ClusterM((CatMZ==alltex(6)),4) = 6;
ClusterM((CatMZ==alltex(7)),4) = 7;
ClusterM((CatMZ==alltex(8)),4) = 8;

ClusterM(:,5) = nan; %Det
ClusterM((CatMD==alltex(1)),5) = 1;
ClusterM((CatMD==alltex(2)),5) = 2;
ClusterM((CatMD==alltex(3)),5) = 3;
ClusterM((CatMD==alltex(4)),5) = 4;
ClusterM((CatMD==alltex(5)),5) = 5;
ClusterM((CatMD==alltex(6)),5) = 6;
ClusterM((CatMD==alltex(7)),5) = 7;
ClusterM((CatMD==alltex(8)),5) = 8;

%% Save
LIDs = Abiom(:,1);
save([spath,'LME_biom_nu_cpue_cme_A_mlr_cluster_ind_drivers.mat'],...
    'ClusterB','ClusterP','ClusterC','ClusterM','alltex','LIDs');

%%  Colormap

% mcol = [
%     34/255 136/255 51/255;...   %green
%     170/255 51/255 119/255;...  %purple
%     238/255 102/255 119/255;... %red
%     0/255 68/255 136/255;...    %blue
%     51/255 187/255 238/255;...  %cyan
%     153/255 153/255 51/255;...  %olive
%     0 0 0;...                   %black
%     0.50 0.50 0.50;...          % grey
%     ];


% colorblind friendly
load('paul_tol_cmaps.mat')

%ecosystem type
ecol = drainbow(2:3:end,:) ./ 255;

%try muted and add greys
mcol = muted ./ 255;
%add greys 
mcol(10,:) = zmeso(10,:);
mcol(11,:) = zmeso(9,:);
mcol(12,:) = zmeso(8,:);

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
% on grid
BclusP  = nan(ni,nj);
PclusP  = nan(ni,nj);
CclusP  = nan(ni,nj);
MclusP  = nan(ni,nj);
BclusB  = nan(ni,nj);
PclusB  = nan(ni,nj);
CclusB  = nan(ni,nj);
MclusB  = nan(ni,nj);
BclusZ  = nan(ni,nj);
PclusZ  = nan(ni,nj);
CclusZ  = nan(ni,nj);
MclusZ  = nan(ni,nj);
BclusD  = nan(ni,nj);
PclusD  = nan(ni,nj);
CclusD  = nan(ni,nj);
MclusD  = nan(ni,nj);

lid = Abiom(:,1);
%use col = 2 after finding matching # to text
for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
    MclusP(id) = ClusterM(i,2);
    PclusP(id) = ClusterP(i,2);
    CclusP(id) = ClusterC(i,2);
    BclusP(id) = ClusterB(i,2);

    MclusB(id) = ClusterM(i,3);
    PclusB(id) = ClusterP(i,3);
    CclusB(id) = ClusterC(i,3);
    BclusB(id) = ClusterB(i,3);

    MclusZ(id) = ClusterM(i,4);
    PclusZ(id) = ClusterP(i,4);
    CclusZ(id) = ClusterC(i,4);
    BclusZ(id) = ClusterB(i,4);

    MclusD(id) = ClusterM(i,5);
    PclusD(id) = ClusterP(i,5);
    CclusD(id) = ClusterC(i,5);
    BclusD(id) = ClusterB(i,5);
end


%% All same 12 colors and text ---------------------------
% TP
f1 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,BclusP)
colormap(mcol(1:8,:))
clim([1 8])
title('Biomass')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,PclusP)
colormap(mcol(1:8,:))
clim([1 8])
title('Prod')

subplot('Position',[0.01 0.10 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,CclusP)
colormap(mcol(1:8,:))
clim([1 8])
title('CPUE')

subplot('Position',[0.33 0.10 0.32 0.4]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,MclusP)
colormap(mcol(1:8,:))
clim([1 8])
title('Catch-Effort regress')
colorbar('Position',[0.675 0.125 0.03 0.35],'Ticks',1:8,'TickLabels',alltex,...
    'Direction','reverse')

%Food web
subplot('Position',[0.65 0.45 0.32 0.4]) 
ax1=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Ebiom)
colormap(ax1,ecol)
caxis([1 5]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('Food web structure')
colorbar('Ticks',1.5:0.75:5,'TickLabels',etex)

print('-dpng',[ppath 'Map_LMEs_driver_mlr_AllFish_TPcluster_fntypes.png'])

%% TB
f2 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,BclusB)
colormap(mcol(1:8,:))
clim([1 8])
title('Biomass')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,PclusB)
colormap(mcol(1:8,:))
clim([1 8])
title('Prod')

subplot('Position',[0.01 0.10 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,CclusB)
colormap(mcol(1:8,:))
clim([1 8])
title('CPUE')

%subplot('Position',[0.65 0.575 0.32 0.4]) 

subplot('Position',[0.33 0.10 0.32 0.4]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,MclusB)
colormap(mcol(1:8,:))
clim([1 8])
title('Catch-Effort regress')
colorbar('Position',[0.675 0.125 0.03 0.35],'Ticks',1:8,'TickLabels',alltex,...
    'Direction','reverse')

%Food web
subplot('Position',[0.65 0.45 0.32 0.4]) 
ax1=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Ebiom)
colormap(ax1,ecol)
caxis([1 5]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('Food web structure')
colorbar('Ticks',1.5:0.75:5,'TickLabels',etex)
print('-dpng',[ppath 'Map_LMEs_driver_mlr_AllFish_TBcluster_fntypes.png'])

%% ZmL ----------

f3 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,BclusZ)
colormap(mcol(1:8,:))
clim([1 8])
title('Biomass')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,PclusZ)
colormap(mcol(1:8,:))
clim([1 8])
title('Prod')

subplot('Position',[0.01 0.10 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,CclusZ)
colormap(mcol(1:8,:))
clim([1 8])
title('CPUE')

%subplot('Position',[0.65 0.575 0.32 0.4]) 

subplot('Position',[0.33 0.10 0.32 0.4]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,MclusZ)
colormap(mcol(1:8,:))
clim([1 8])
title('Catch-Effort regress')
colorbar('Position',[0.675 0.125 0.03 0.35],'Ticks',1:8,'TickLabels',alltex,...
    'Direction','reverse')

%Food web
subplot('Position',[0.65 0.45 0.32 0.4]) 
ax1=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Ebiom)
colormap(ax1,ecol)
caxis([1 5]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('Food web structure')
colorbar('Ticks',1.5:0.75:5,'TickLabels',etex)
print('-dpng',[ppath 'Map_LMEs_driver_mlr_AllFish_ZmLcluster_fntypes.png'])

%% Det --------------------------

f4 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.01 0.575 0.32 0.4]) %F
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,BclusD)
colormap(mcol(1:8,:))
clim([1 8])
title('Biomass')

subplot('Position',[0.33 0.575 0.32 0.4]) %P
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,PclusD)
colormap(mcol(1:8,:))
clim([1 8])
title('Prod')

subplot('Position',[0.01 0.10 0.32 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,CclusD)
colormap(mcol(1:8,:))
clim([1 8])
title('CPUE')

%subplot('Position',[0.65 0.575 0.32 0.4]) 

subplot('Position',[0.33 0.10 0.32 0.4]) %B
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.95 0.95 0.95]);
hold on
surfm(TLAT,TLONG,MclusD)
colormap(mcol(1:8,:))
clim([1 8])
title('Catch-Effort regress')
colorbar('Position',[0.675 0.125 0.03 0.35],'Ticks',1:8,'TickLabels',alltex,...
    'Direction','reverse')

%Food web
subplot('Position',[0.65 0.45 0.32 0.4]) 
ax1=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Ebiom)
colormap(ax1,ecol)
caxis([1 5]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('Food web structure')
colorbar('Ticks',1.5:0.75:5,'TickLabels',etex)
print('-dpng',[ppath 'Map_LMEs_driver_mlr_AllFish_Detcluster_fntypes.png'])






