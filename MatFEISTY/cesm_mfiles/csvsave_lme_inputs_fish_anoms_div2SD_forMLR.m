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

%% Divide anoms by 2SD
atp = atp ./ (2*std(atp,0,2));
atb = atb ./ (2*std(atb,0,2));
adety = adety ./ (2*std(adety,0,2));
azoo = azoo ./ (2*std(azoo,0,2));
azlosy = azlosy ./ (2*std(azlosy,0,2));

%Fish
as = as ./ (2*std(as,0,2));
am = am ./ (2*std(am,0,2));
al = al ./ (2*std(al,0,2));
af = af ./ (2*std(af,0,2));
ap = ap ./ (2*std(ap,0,2));
ad = ad ./ (2*std(ad,0,2));
aa = aa ./ (2*std(aa,0,2));
ab = ab ./ (2*std(ab,0,2));

%% Transpose
atp = atp';
atb = atb';
adety = adety';
azoo = azoo';
azlosy = azlosy';

%Fish
as = as';
am = am';
al = al';
af = af';
ap = ap';
ad = ad';
aa = aa';
ab = ab';

%%
yr = [1:68]';
yname = cellstr(num2str(yr));

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = aa(1,:);
lid = find(~isnan(AA));

for l=1:length(lid)
    lname{l} = ['LME' num2str(lid(l))];
end

%% Drivers
TPtab1 = array2table([yr,atp(:,lid)],"VariableNames",['Year' lname]);

TBtab1 = array2table([yr,atb(:,lid)],"VariableNames",['Year' lname]);

DEtab1 = array2table([yr,adety(:,lid)],"VariableNames",['Year' lname]);

ZBtab1 = array2table([yr,azoo(:,lid)],"VariableNames",['Year' lname]);

ZLtab1 = array2table([yr,azlosy(:,lid)],"VariableNames",['Year' lname]);

% Fish
Stab1 = array2table([yr,as(:,lid)],"VariableNames",['Year' lname]);

Mtab1 = array2table([yr,am(:,lid)],"VariableNames",['Year' lname]);

Ltab1 = array2table([yr,al(:,lid)],"VariableNames",['Year' lname]);

Ftab1 = array2table([yr,af(:,lid)],"VariableNames",['Year' lname]);

Ptab1 = array2table([yr,ap(:,lid)],"VariableNames",['Year' lname]);

Dtab1 = array2table([yr,ad(:,lid)],"VariableNames",['Year' lname]);

Atab1 = array2table([yr,aa(:,lid)],"VariableNames",['Year' lname]);

Btab1 = array2table([yr,ab(:,lid)],"VariableNames",['Year' lname]);

%%
writetable(TPtab1,[fpath,'LMEs_anoms_div2SD_TP.csv'],...
    'Delimiter',',');
writetable(TBtab1,[fpath,'LMEs_anoms_div2SD_TB.csv'],...
    'Delimiter',',');
writetable(DEtab1,[fpath,'LMEs_anoms_div2SD_Det.csv'],...
    'Delimiter',',');
writetable(ZBtab1,[fpath,'LMEs_anoms_div2SD_Zmeso.csv'],...
    'Delimiter',',');
writetable(ZLtab1,[fpath,'LMEs_anoms_div2SD_ZmLoss.csv'],...
    'Delimiter',',');

writetable(Stab1,[fpath,'LMEs_anoms_div2SD_S.csv'],...
    'Delimiter',',');
writetable(Mtab1,[fpath,'LMEs_anoms_div2SD_M.csv'],...
    'Delimiter',',');
writetable(Ltab1,[fpath,'LMEs_anoms_div2SD_L.csv'],...
    'Delimiter',',');
writetable(Ftab1,[fpath,'LMEs_anoms_div2SD_F.csv'],...
    'Delimiter',',');
writetable(Ptab1,[fpath,'LMEs_anoms_div2SD_P.csv'],...
    'Delimiter',',');
writetable(Dtab1,[fpath,'LMEs_anoms_div2SD_D.csv'],...
    'Delimiter',',');
writetable(Atab1,[fpath,'LMEs_anoms_div2SD_A.csv'],...
    'Delimiter',',');
writetable(Btab1,[fpath,'LMEs_anoms_div2SD_B.csv'],...
    'Delimiter',',');


