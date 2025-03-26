% Table of most sig var in driver-catch corrs
% Lag with max R2
% For all 63 LMEs
% Const & Obs fishing effort

clear
close all

%% % ------------------------------------------------------------
%cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
cpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/cpue2015/';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

% fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/cpue2015/';
% spath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/cpue2015/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regress_cpue/'];

mod = 'v15_All_fish03';
mod2 = 'v15_obsfish2015';

cnam = {'coef','p','lag','idriver','driver'};

%%  ---------------- sat --------------------------
load([spath,'LMEs_corr_catch_chlyrs15_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')
load([spath 'LMEs_corr_cpue_chlyrs15_feisty_lags.mat'],'lid')

stex = {'SST','Chl'};

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
LAtab = LAtab(lid,:);
LFtab = LFtab(lid,:);
LPtab = LPtab(lid,:);
LDtab = LDtab(lid,:);

sigA = (LAtab(:,2) <= 0.05);
sigF = (LFtab(:,2) <= 0.05);
sigP = (LPtab(:,2) <= 0.05);
sigD = (LDtab(:,2) <= 0.05);

AtabS = LAtab(sigA,:);
FtabS = LFtab(sigF,:);
PtabS = LPtab(sigP,:);
DtabS = LDtab(sigD,:);

clear sigA sigF sigP sigD LAtab LFtab LPtab LDtab

%%  ---------------- drivers --------------------------
load([spath,'LMEs_corr_catch_chlyrs15_inputs_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

dtex = {'TP','TB','Det','ZmLoss','SST','Chl'};

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
LAtab = LAtab(lid,:);
LFtab = LFtab(lid,:);
LPtab = LPtab(lid,:);
LDtab = LDtab(lid,:);

sigA = (LAtab(:,2) <= 0.05);
sigF = (LFtab(:,2) <= 0.05);
sigP = (LPtab(:,2) <= 0.05);
sigD = (LDtab(:,2) <= 0.05);

AtabD = LAtab(sigA,:);
FtabD = LFtab(sigF,:);
PtabD = LPtab(sigP,:);
DtabD = LDtab(sigD,:);

clear sigA sigF sigP sigD LAtab LFtab LPtab LDtab

%%  ---------------- constfish --------------------------
load([spath,'LMEs_corr_catch_chlyrs15_inputs_feisty_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

ftex = {'TP','TB','Det','ZmLoss','SST','Chl','Biom','Prod'};

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
sigA = (LAtab(:,2) <= 0.05);
sigF = (LFtab(:,2) <= 0.05);
sigP = (LPtab(:,2) <= 0.05);
sigD = (LDtab(:,2) <= 0.05);

AtabF = LAtab(sigA,:);
FtabF = LFtab(sigF,:);
PtabF = LPtab(sigP,:);
DtabF = LDtab(sigD,:);

clear sigA sigF sigP sigD LAtab LFtab LPtab LDtab

%%  ---------------- obsfish --------------------------
load([spath,'LMEs_corr_catch_chlyrs15_inputs_obsfish2015_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

otex = {'TP','TB','Det','ZmLoss','SST','Chl','Biom','Prod','Yield'};

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
sigA = (LAtab(:,2) <= 0.05);
sigF = (LFtab(:,2) <= 0.05);
sigP = (LPtab(:,2) <= 0.05);
sigD = (LDtab(:,2) <= 0.05);

AtabO = LAtab(sigA,:);
FtabO = LFtab(sigF,:);
PtabO = LPtab(sigP,:);
DtabO = LDtab(sigD,:);

clear sigA sigF sigP sigD LAtab LFtab LPtab LDtab

%%

%All fishes
matA(1,1) = sum(AtabS(:,4)==1);
matA(2,1) = sum(AtabS(:,4)==2);

matA(1,2) = sum(AtabD(:,4)==5);
matA(2,2) = sum(AtabD(:,4)==6);
matA(3,2) = sum(AtabD(:,4)==1);
matA(4,2) = sum(AtabD(:,4)==2);
matA(5,2) = sum(AtabD(:,4)==3);
matA(6,2) = sum(AtabD(:,4)==4);

matA(1,3) = sum(AtabF(:,4)==5);
matA(2,3) = sum(AtabF(:,4)==6);
matA(3,3) = sum(AtabF(:,4)==1);
matA(4,3) = sum(AtabF(:,4)==2);
matA(5,3) = sum(AtabF(:,4)==3);
matA(6,3) = sum(AtabF(:,4)==4);
matA(7,3) = sum(AtabF(:,4)==7);
matA(8,3) = sum(AtabF(:,4)==8);

matA(1,4) = sum(AtabO(:,4)==5);
matA(2,4) = sum(AtabO(:,4)==6);
matA(3,4) = sum(AtabO(:,4)==1);
matA(4,4) = sum(AtabO(:,4)==2);
matA(5,4) = sum(AtabO(:,4)==3);
matA(6,4) = sum(AtabO(:,4)==4);
matA(7,4) = sum(AtabO(:,4)==7);
matA(8,4) = sum(AtabO(:,4)==8);
matA(9,4) = sum(AtabO(:,4)==9);

%Forage
matF(1,1) = sum(FtabS(:,4)==1);
matF(2,1) = sum(FtabS(:,4)==2);

matF(1,2) = sum(FtabD(:,4)==5);
matF(2,2) = sum(FtabD(:,4)==6);
matF(3,2) = sum(FtabD(:,4)==1);
matF(4,2) = sum(FtabD(:,4)==2);
matF(5,2) = sum(FtabD(:,4)==3);
matF(6,2) = sum(FtabD(:,4)==4);

matF(1,3) = sum(FtabF(:,4)==5);
matF(2,3) = sum(FtabF(:,4)==6);
matF(3,3) = sum(FtabF(:,4)==1);
matF(4,3) = sum(FtabF(:,4)==2);
matF(5,3) = sum(FtabF(:,4)==3);
matF(6,3) = sum(FtabF(:,4)==4);
matF(7,3) = sum(FtabF(:,4)==7);
matF(8,3) = sum(FtabF(:,4)==8);

matF(1,4) = sum(FtabO(:,4)==5);
matF(2,4) = sum(FtabO(:,4)==6);
matF(3,4) = sum(FtabO(:,4)==1);
matF(4,4) = sum(FtabO(:,4)==2);
matF(5,4) = sum(FtabO(:,4)==3);
matF(6,4) = sum(FtabO(:,4)==4);
matF(7,4) = sum(FtabO(:,4)==7);
matF(8,4) = sum(FtabO(:,4)==8);
matF(9,4) = sum(FtabO(:,4)==9);

%Lg Pel
matP(1,1) = sum(PtabS(:,4)==1);
matP(2,1) = sum(PtabS(:,4)==2);

matP(1,2) = sum(PtabD(:,4)==5);
matP(2,2) = sum(PtabD(:,4)==6);
matP(3,2) = sum(PtabD(:,4)==1);
matP(4,2) = sum(PtabD(:,4)==2);
matP(5,2) = sum(PtabD(:,4)==3);
matP(6,2) = sum(PtabD(:,4)==4);

matP(1,3) = sum(PtabF(:,4)==5);
matP(2,3) = sum(PtabF(:,4)==6);
matP(3,3) = sum(PtabF(:,4)==1);
matP(4,3) = sum(PtabF(:,4)==2);
matP(5,3) = sum(PtabF(:,4)==3);
matP(6,3) = sum(PtabF(:,4)==4);
matP(7,3) = sum(PtabF(:,4)==7);
matP(8,3) = sum(PtabF(:,4)==8);

matP(1,4) = sum(PtabO(:,4)==5);
matP(2,4) = sum(PtabO(:,4)==6);
matP(3,4) = sum(PtabO(:,4)==1);
matP(4,4) = sum(PtabO(:,4)==2);
matP(5,4) = sum(PtabO(:,4)==3);
matP(6,4) = sum(PtabO(:,4)==4);
matP(7,4) = sum(PtabO(:,4)==7);
matP(8,4) = sum(PtabO(:,4)==8);
matP(9,4) = sum(PtabO(:,4)==9);

%Dem
matD(1,1) = sum(DtabS(:,4)==1);
matD(2,1) = sum(DtabS(:,4)==2);

matD(1,2) = sum(DtabD(:,4)==5);
matD(2,2) = sum(DtabD(:,4)==6);
matD(3,2) = sum(DtabD(:,4)==1);
matD(4,2) = sum(DtabD(:,4)==2);
matD(5,2) = sum(DtabD(:,4)==3);
matD(6,2) = sum(DtabD(:,4)==4);

matD(1,3) = sum(DtabF(:,4)==5);
matD(2,3) = sum(DtabF(:,4)==6);
matD(3,3) = sum(DtabF(:,4)==1);
matD(4,3) = sum(DtabF(:,4)==2);
matD(5,3) = sum(DtabF(:,4)==3);
matD(6,3) = sum(DtabF(:,4)==4);
matD(7,3) = sum(DtabF(:,4)==7);
matD(8,3) = sum(DtabF(:,4)==8);

matD(1,4) = sum(DtabO(:,4)==5);
matD(2,4) = sum(DtabO(:,4)==6);
matD(3,4) = sum(DtabO(:,4)==1);
matD(4,4) = sum(DtabO(:,4)==2);
matD(5,4) = sum(DtabO(:,4)==3);
matD(6,4) = sum(DtabO(:,4)==4);
matD(7,4) = sum(DtabO(:,4)==7);
matD(8,4) = sum(DtabO(:,4)==8);
matD(9,4) = sum(DtabO(:,4)==9);

%%
lname = cellstr(num2str(lid));
cname = {'Sat','Sat-OBGC','Sat-OBGC-FishC','Sat-OBGC-FishO'};
rname = {'SST','Chl','TP','TB','Det','ZmLoss','Biom','Prod','Yield'};

% cnam, lname, tanom2
TabA = array2table(matA,"RowNames",rname,"VariableNames",cname);

TabF = array2table(matF,"RowNames",rname,"VariableNames",cname);

TabP = array2table(matP,"RowNames",rname,"VariableNames",cname);

TabD = array2table(matD,"RowNames",rname,"VariableNames",cname);

%%
writetable(TabA,[spath,'Num_LMEs_corr_catch_chlyrs15_inputs_feisty_maxcorr_posfood_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(TabF,[spath,'Num_LMEs_corr_catch_chlyrs15_inputs_feisty_maxcorr_posfood_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(TabP,[spath,'Num_LMEs_corr_catch_chlyrs15_inputs_feisty_maxcorr_posfood_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(TabD,[spath,'Num_LMEs_corr_catch_chlyrs15_inputs_feisty_maxcorr_posfood_D.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'Num_LMEs_corr_catch_chlyrs15_inputs_feisty_maxcorr_posfood.mat'],...
    'matF','matP','matD','matA',...
    'TabF','TabP','TabD','TabA','lid');


