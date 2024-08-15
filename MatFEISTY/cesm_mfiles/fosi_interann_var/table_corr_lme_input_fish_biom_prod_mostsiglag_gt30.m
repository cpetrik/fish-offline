% Table of most sig var in driver-fish corrs
% Biomass and Prod together
% Lag with max R2
% For all 63 LMEs
% Const effort only

clear 
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',cfile,'/corrs/'];

mod = 'v15_All_fish03_';

%spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';

%% Biomass
load([spath,'LMEs_corr_driver_maxcorrs.mat'],...
    'LAtab','LFtab','LPtab','LDtab','lid')

sigAB = (LAtab(:,2) <= 0.05);
sigFB = (LFtab(:,2) <= 0.05);
sigPB = (LPtab(:,2) <= 0.05);
sigDB = (LDtab(:,2) <= 0.05);

BAtab = LAtab(sigAB,:);
BFtab = LFtab(sigFB,:);
BPtab = LPtab(sigPB,:);
BDtab = LDtab(sigDB,:);

clear LAtab LFtab LPtab LDtab

%% Prod
load([spath,'LMEs_nu_corr_driver_maxcorrs.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

sigAP = (LAtab(:,2) <= 0.05);
sigFP = (LFtab(:,2) <= 0.05);
sigPP = (LPtab(:,2) <= 0.05);
sigDP = (LDtab(:,2) <= 0.05);

PAtab = LAtab(sigAP,:);
PFtab = LFtab(sigFP,:);
PPtab = LPtab(sigPP,:);
PDtab = LDtab(sigDP,:);

clear LAtab LFtab LPtab LDtab

%% LMEs with >30% biomass
load([fpath 'LME_fosi_fished_',mod,cfile '.mat'],'lme_area','lme_mtype')

lme_type = lme_mtype(:,1:3);
lme_btot = sum(lme_type,2);

lbio = lme_type ./ repmat(lme_btot,1,3);

%remove NaNs are 23, 33, 62 (inland seas)
iis = setdiff(1:66,[23, 33, 62]);
% Relative biomass
rbio= lbio(iis,:);

domF = (rbio(:,1)>0.3);
domP = (rbio(sigPB,2)>0.3);
domD = (rbio(:,3)>0.3);

domPP = (rbio(sigPP,2)>0.3);

nF = sum(domF) %37
nP = sum(domP) %12
nD = sum(domD) %34

%  drivers 
dtex = {'TP','TB','Det','ZmLoss'};

%% All results
%All fishes
matA(1,1) = sum(BAtab(:,4)==1);
matA(2,1) = sum(BAtab(:,4)==2);
matA(3,1) = sum(BAtab(:,4)==3);
matA(4,1) = sum(BAtab(:,4)==4);

matA(1,2) = sum(PAtab(:,4)==1);
matA(2,2) = sum(PAtab(:,4)==2);
matA(3,2) = sum(PAtab(:,4)==3);
matA(4,2) = sum(PAtab(:,4)==4);

%Forage
matF(1,1) = sum(BFtab(:,4)==1);
matF(2,1) = sum(BFtab(:,4)==2);
matF(3,1) = sum(BFtab(:,4)==3);
matF(4,1) = sum(BFtab(:,4)==4);

matF(1,2) = sum(PFtab(:,4)==1);
matF(2,2) = sum(PFtab(:,4)==2);
matF(3,2) = sum(PFtab(:,4)==3);
matF(4,2) = sum(PFtab(:,4)==4);

%Lg Pel
matP(1,1) = sum(BPtab(:,4)==1);
matP(2,1) = sum(BPtab(:,4)==2);
matP(3,1) = sum(BPtab(:,4)==3);
matP(4,1) = sum(BPtab(:,4)==4);

matP(1,2) = sum(PPtab(:,4)==1);
matP(2,2) = sum(PPtab(:,4)==2);
matP(3,2) = sum(PPtab(:,4)==3);
matP(4,2) = sum(PPtab(:,4)==4);

%Dem
matD(1,1) = sum(BDtab(:,4)==1);
matD(2,1) = sum(BDtab(:,4)==2);
matD(3,1) = sum(BDtab(:,4)==3);
matD(4,1) = sum(BDtab(:,4)==4);

matD(1,2) = sum(PDtab(:,4)==1);
matD(2,2) = sum(PDtab(:,4)==2);
matD(3,2) = sum(PDtab(:,4)==3);
matD(4,2) = sum(PDtab(:,4)==4);

%% Only fn types if >30% biomass
%Forage
matF3(1,1) = sum(BFtab(domF,4)==1);
matF3(2,1) = sum(BFtab(domF,4)==2);
matF3(3,1) = sum(BFtab(domF,4)==3);
matF3(4,1) = sum(BFtab(domF,4)==4);

matF3(1,2) = sum(PFtab(domF,4)==1);
matF3(2,2) = sum(PFtab(domF,4)==2);
matF3(3,2) = sum(PFtab(domF,4)==3);
matF3(4,2) = sum(PFtab(domF,4)==4);

%Lg Pel
matP3(1,1) = sum(BPtab(domP,4)==1);
matP3(2,1) = sum(BPtab(domP,4)==2);
matP3(3,1) = sum(BPtab(domP,4)==3);
matP3(4,1) = sum(BPtab(domP,4)==4);

matP3(1,2) = sum(PPtab(domPP,4)==1);
matP3(2,2) = sum(PPtab(domPP,4)==2);
matP3(3,2) = sum(PPtab(domPP,4)==3);
matP3(4,2) = sum(PPtab(domPP,4)==4);

%Dem
matD3(1,1) = sum(BDtab(domD,4)==1);
matD3(2,1) = sum(BDtab(domD,4)==2);
matD3(3,1) = sum(BDtab(domD,4)==3);
matD3(4,1) = sum(BDtab(domD,4)==4);

matD3(1,2) = sum(PDtab(domD,4)==1);
matD3(2,2) = sum(PDtab(domD,4)==2);
matD3(3,2) = sum(PDtab(domD,4)==3);
matD3(4,2) = sum(PDtab(domD,4)==4);

%% R2 > 0.5 = |r| > |0.7071|
matR2(1,1) = sum((abs(BAtab(:,1)))>0.707);
matR2(2,1) = sum((abs(BFtab(:,1)))>0.707);
matR2(3,1) = sum((abs(BPtab(:,1)))>0.707);
matR2(4,1) = sum((abs(BDtab(:,1)))>0.707);

matR2(1,2) = sum((abs(PAtab(:,1)))>0.707);
matR2(2,2) = sum((abs(PFtab(:,1)))>0.707);
matR2(3,2) = sum((abs(PPtab(:,1)))>0.707);
matR2(4,2) = sum((abs(PDtab(:,1)))>0.707);

%%
lname = cellstr(num2str(lid));
cname = {'A Biom','A Prod','F Biom','F Prod',...
    'P Biom','P Prod','D Biom','D Prod'};
cname2 = {'A Biom','F Biom','P Biom','D Biom',...
    'A Prod','F Prod','P Prod','D Prod'};
rname = {'TP','TB','Det','ZmLoss'};

% cnam, lname, tanom2
matW = [matA,matF,matP,matD];
TabW = array2table(matW,"RowNames",rname,"VariableNames",cname);

%Best one - grouped by fn type
matL = [matA';matF';matP';matD'];
TabL = array2table(matL,"RowNames",cname,"VariableNames",rname);

%Grouped by biomass or prod
matL2 = matL(1:2:end,:);
matL2(5:8,:) = matL(2:2:end,:);
TabL2 = array2table(matL2,"RowNames",cname2,"VariableNames",rname);

matL30 = [matA';matF3';matP3';matD3'];
TabL30 = array2table(matL30,"RowNames",cname,"VariableNames",rname);

ftext = {'All','Forage','LgPel','Dem'};
vtext = {'Biomass','Production'};
TabR2 = array2table(matR2,"RowNames",ftext,"VariableNames",vtext);

%%
writetable(TabL,[spath,'Num_LMEs_corr_driver_feisty_fntype.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(TabL30,[spath,'Num_LMEs_corr_driver_feisty_fntype_gt30.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(TabR2,[spath,'Num_LMEs_corr_driver_feisty_fntype_R2gt50.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'Num_LMEs_corr_driver_feisty_biom_nu.mat'],...
    'matF','matP','matD','matA',...
    'matP3','matD3','matF3',...
    'matL','matL30','TabL','TabL30','lid');


