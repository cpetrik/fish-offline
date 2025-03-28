% Save LME anoms of nu of fn types
% To use in mult linear regression

clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
%fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

% Anoms with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms.mat'],...
    'aa','ad','af','ap')

%% Divide anoms by 2SD
%Fish
af = af ./ (2*std(af,0,2));
ap = ap ./ (2*std(ap,0,2));
ad = ad ./ (2*std(ad,0,2));
aa = aa ./ (2*std(aa,0,2));

%% Transpose
af = af';
ap = ap';
ad = ad';
aa = aa';

%%
yr = [1:68]';
yname = cellstr(num2str(yr));

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = aa(1,:);
lid = find(~isnan(AA));

for l=1:length(lid)
    lname{l} = ['LME' num2str(lid(l))];
end

%% Fish
Ftab1 = array2table([yr,af(:,lid)],"VariableNames",['Year' lname]);

Ptab1 = array2table([yr,ap(:,lid)],"VariableNames",['Year' lname]);

Dtab1 = array2table([yr,ad(:,lid)],"VariableNames",['Year' lname]);

Atab1 = array2table([yr,aa(:,lid)],"VariableNames",['Year' lname]);

%%
writetable(Ftab1,[fpath,'LMEs_anoms_div2SD_Fnu.csv'],...
    'Delimiter',',');
writetable(Ptab1,[fpath,'LMEs_anoms_div2SD_Pnu.csv'],...
    'Delimiter',',');
writetable(Dtab1,[fpath,'LMEs_anoms_div2SD_Dnu.csv'],...
    'Delimiter',',');
writetable(Atab1,[fpath,'LMEs_anoms_div2SD_Anu.csv'],...
    'Delimiter',',');

