% Save LME anoms of nu of fn types
% To use in mult linear regression

clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

% Anoms with linear trend removed
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms.mat'],...
    'aa','ad','af','ap')

%%
yr = [1:68]';
yname = cellstr(num2str(yr));

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = aa(:,1);
lid = find(~isnan(AA));

lname = cellstr(num2str(lid));

% Fish
Ftab1 = array2table(af(lid,:),"RowNames",lname);
Ftab1.Properties.VariableNames = yname;

Ptab1 = array2table(ap(lid,:),"RowNames",lname);
Ptab1.Properties.VariableNames = yname;

Dtab1 = array2table(ad(lid,:),"RowNames",lname);
Dtab1.Properties.VariableNames = yname;

Atab1 = array2table(aa(lid,:),"RowNames",lname);
Atab1.Properties.VariableNames = yname;

%%
writetable(Ftab1,[fpath,'LMEs_anoms_Fnu.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[fpath,'LMEs_anoms_Pnu.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[fpath,'LMEs_anoms_Dnu.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Atab1,[fpath,'LMEs_anoms_Anu.csv'],...
    'Delimiter',',','WriteRowNames',true);

