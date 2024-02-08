% Save LME anoms of FishMIP catch of fn types
% To use in mult linear regression

clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
sims = {'v15_All_fish03_';'v15_climatol_';'v15_varTemp_';'v15_varFood_'};
mod = sims{1};

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

load([fpath 'FEISTY_FOSI_',mod,'lme_gam_rec_ann_mean_anoms.mat'],...
    'agf','agp','agd','aga',...
    'arf','arp','ard','ara','units');

%% Divide anoms by 2SD
%gam
agf = agf ./ (2*std(agf,0,2));
agp = agp ./ (2*std(agp,0,2));
agd = agd ./ (2*std(agd,0,2));
aga = aga ./ (2*std(aga,0,2));
%rec
arf = arf ./ (2*std(arf,0,2));
arp = arp ./ (2*std(arp,0,2));
ard = ard ./ (2*std(ard,0,2));
ara = ara ./ (2*std(ara,0,2));

%% Transpose
agf = agf';
agp = agp';
agd = agd';
aga = aga';

arf = arf';
arp = arp';
ard = ard';
ara = ara';

%%
yr = [1948:2015]';
yname = cellstr(num2str(yr));

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
lid = [1:22,24:32,34:61,63:66];

for l=1:length(lid)
    lname{l} = ['LME' num2str(lid(l))];
end

%% Fish
gFtab1 = array2table([yr,agf(:,lid)],"VariableNames",['Year' lname]);
gPtab1 = array2table([yr,agp(:,lid)],"VariableNames",['Year' lname]);
gDtab1 = array2table([yr,agd(:,lid)],"VariableNames",['Year' lname]);
gAtab1 = array2table([yr,aga(:,lid)],"VariableNames",['Year' lname]);

rFtab1 = array2table([yr,arf(:,lid)],"VariableNames",['Year' lname]);
rPtab1 = array2table([yr,arp(:,lid)],"VariableNames",['Year' lname]);
rDtab1 = array2table([yr,ard(:,lid)],"VariableNames",['Year' lname]);
rAtab1 = array2table([yr,ara(:,lid)],"VariableNames",['Year' lname]);

%%
writetable(gFtab1,[fpath,'LMEs_FishMIP_gam_anoms_div2SD_F.csv'],...
    'Delimiter',',');
writetable(gPtab1,[fpath,'LMEs_FishMIP_gam_anoms_div2SD_P.csv'],...
    'Delimiter',',');
writetable(gDtab1,[fpath,'LMEs_FishMIP_gam_anoms_div2SD_D.csv'],...
    'Delimiter',',');
writetable(gAtab1,[fpath,'LMEs_FishMIP_gam_anoms_div2SD_A.csv'],...
    'Delimiter',',');

writetable(rFtab1,[fpath,'LMEs_FishMIP_rec_anoms_div2SD_F.csv'],...
    'Delimiter',',');
writetable(rPtab1,[fpath,'LMEs_FishMIP_rec_anoms_div2SD_P.csv'],...
    'Delimiter',',');
writetable(rDtab1,[fpath,'LMEs_FishMIP_rec_anoms_div2SD_D.csv'],...
    'Delimiter',',');
writetable(rAtab1,[fpath,'LMEs_FishMIP_rec_anoms_div2SD_A.csv'],...
    'Delimiter',',');

