% Save LME anoms of driver and fn types
% To use in mult linear regression

clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";

%% FOSI input forcing

%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adety','atb','atp','azlosy','azoo');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

tanom = {'TP','TB','Det','Zmeso','ZmLoss'};

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

% Anoms with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'_lme_ann_mean_anoms.mat'],...
    'aa','ab','ad','af','ap','as','am','al')

%%
yr = [1:68]';
yname = cellstr(num2str(yr));

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = aa(:,1);
lid = find(~isnan(AA));

lname = cellstr(num2str(lid));

% Drivers
TPtab1 = array2table(atp(lid,:),"RowNames",lname);
TPtab1.Properties.VariableNames = yname;

TBtab1 = array2table(atb(lid,:),"RowNames",lname);
TBtab1.Properties.VariableNames = yname;

DEtab1 = array2table(adety(lid,:),"RowNames",lname);
DEtab1.Properties.VariableNames = yname;

ZBtab1 = array2table(azoo(lid,:),"RowNames",lname);
ZBtab1.Properties.VariableNames = yname;

ZLtab1 = array2table(azlosy(lid,:),"RowNames",lname);
ZLtab1.Properties.VariableNames = yname;

% Fish
Stab1 = array2table(as(lid,:),"RowNames",lname);
Stab1.Properties.VariableNames = yname;

Mtab1 = array2table(am(lid,:),"RowNames",lname);
Mtab1.Properties.VariableNames = yname;

Ltab1 = array2table(al(lid,:),"RowNames",lname);
Ltab1.Properties.VariableNames = yname;

Ftab1 = array2table(af(lid,:),"RowNames",lname);
Ftab1.Properties.VariableNames = yname;

Ptab1 = array2table(ap(lid,:),"RowNames",lname);
Ptab1.Properties.VariableNames = yname;

Dtab1 = array2table(ad(lid,:),"RowNames",lname);
Dtab1.Properties.VariableNames = yname;

Atab1 = array2table(aa(lid,:),"RowNames",lname);
Atab1.Properties.VariableNames = yname;

Btab1 = array2table(ab(lid,:),"RowNames",lname);
Btab1.Properties.VariableNames = yname;

%%
writetable(TPtab1,[fpath,'LMEs_anoms_TP.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(TBtab1,[fpath,'LMEs_anoms_TB.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(DEtab1,[fpath,'LMEs_anoms_Det.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(ZBtab1,[fpath,'LMEs_anoms_Zmeso.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(ZLtab1,[fpath,'LMEs_anoms_ZmLoss.csv'],...
    'Delimiter',',','WriteRowNames',true);

writetable(Stab1,[fpath,'LMEs_anoms_S.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Mtab1,[fpath,'LMEs_anoms_M.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ltab1,[fpath,'LMEs_anoms_L.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ftab1,[fpath,'LMEs_anoms_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[fpath,'LMEs_anoms_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[fpath,'LMEs_anoms_D.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Atab1,[fpath,'LMEs_anoms_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Btab1,[fpath,'LMEs_anoms_B.csv'],...
    'Delimiter',',','WriteRowNames',true);


