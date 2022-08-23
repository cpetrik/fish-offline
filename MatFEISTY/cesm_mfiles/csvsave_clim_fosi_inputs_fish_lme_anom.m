% save lme ts anomalies as csv files for R

clear all 
close all

%% FOSI input forcing
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat']);
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat']);

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
ppath = [pp cfile '/corrs/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

sims = {'v15_All_fish03_';'v15_climatol_';'v15_varFood_';'v15_varTemp_'};
mod = sims{1};

load([fpath 'FEISTY_FOSI_',mod,'lme_ann_mean_anoms.mat']) % Anoms with linear trend removed

%% put into one matrix per LME

lid = [54,1:2,10,3,5:7];
lname = {'CHK','EBS','GAK','HI','CCE','GMX','SE','NE'};

CHanom(:,1) = [1948:2015]';
BSanom(:,1) = [1948:2015]';
AKanom(:,1) = [1948:2015]';
HIanom(:,1) = [1948:2015]';
CCanom(:,1) = [1948:2015]';
MXanom(:,1) = [1948:2015]';
SEanom(:,1) = [1948:2015]';
NEanom(:,1) = [1948:2015]';

%% forcing inputs
iname = {'Tp','Tb','Det','LZbiom','LZloss'};

CHanom(:,2) = atp(lid(1),:);
CHanom(:,3) = atb(lid(1),:);
CHanom(:,4) = adet(lid(1),:);
CHanom(:,5) = azoo(lid(1),:);
CHanom(:,6) = azlos(lid(1),:);

BSanom(:,2) = atp(lid(2),:);
BSanom(:,3) = atb(lid(2),:);
BSanom(:,4) = adet(lid(2),:);
BSanom(:,5) = azoo(lid(2),:);
BSanom(:,6) = azlos(lid(2),:);

AKanom(:,2) = atp(lid(3),:);
AKanom(:,3) = atb(lid(3),:);
AKanom(:,4) = adet(lid(3),:);
AKanom(:,5) = azoo(lid(3),:);
AKanom(:,6) = azlos(lid(3),:);

HIanom(:,2) = atp(lid(4),:);
HIanom(:,3) = atb(lid(4),:);
HIanom(:,4) = adet(lid(4),:);
HIanom(:,5) = azoo(lid(4),:);
HIanom(:,6) = azlos(lid(4),:);

CCanom(:,2) = atp(lid(5),:);
CCanom(:,3) = atb(lid(5),:);
CCanom(:,4) = adet(lid(5),:);
CCanom(:,5) = azoo(lid(5),:);
CCanom(:,6) = azlos(lid(5),:);

MXanom(:,2) = atp(lid(6),:);
MXanom(:,3) = atb(lid(6),:);
MXanom(:,4) = adet(lid(6),:);
MXanom(:,5) = azoo(lid(6),:);
MXanom(:,6) = azlos(lid(6),:);

SEanom(:,2) = atp(lid(7),:);
SEanom(:,3) = atb(lid(7),:);
SEanom(:,4) = adet(lid(7),:);
SEanom(:,5) = azoo(lid(7),:);
SEanom(:,6) = azlos(lid(7),:);

NEanom(:,2) = atp(lid(8),:);
NEanom(:,3) = atb(lid(8),:);
NEanom(:,4) = adet(lid(8),:);
NEanom(:,5) = azoo(lid(8),:);
NEanom(:,6) = azlos(lid(8),:);

%% FEISTY outputs grouped
gname = {'S','M','L','F','P','D','A','B'};

CHanom(:,7)  = as(lid(1),:);
CHanom(:,8)  = am(lid(1),:);
CHanom(:,9)  = al(lid(1),:);
CHanom(:,10) = af(lid(1),:);
CHanom(:,11) = ap(lid(1),:);
CHanom(:,12) = ad(lid(1),:);
CHanom(:,13) = aa(lid(1),:);
CHanom(:,14) = ab(lid(1),:);

BSanom(:,7)  = as(lid(2),:);
BSanom(:,8)  = am(lid(2),:);
BSanom(:,9)  = al(lid(2),:);
BSanom(:,10) = af(lid(2),:);
BSanom(:,11) = ap(lid(2),:);
BSanom(:,12) = ad(lid(2),:);
BSanom(:,13) = aa(lid(2),:);
BSanom(:,14) = ab(lid(2),:);

AKanom(:,7)  = as(lid(3),:);
AKanom(:,8)  = am(lid(3),:);
AKanom(:,9)  = al(lid(3),:);
AKanom(:,10) = af(lid(3),:);
AKanom(:,11) = ap(lid(3),:);
AKanom(:,12) = ad(lid(3),:);
AKanom(:,13) = aa(lid(3),:);
AKanom(:,14) = ab(lid(3),:);

HIanom(:,7)  = as(lid(4),:);
HIanom(:,8)  = am(lid(4),:);
HIanom(:,9)  = al(lid(4),:);
HIanom(:,10) = af(lid(4),:);
HIanom(:,11) = ap(lid(4),:);
HIanom(:,12) = ad(lid(4),:);
HIanom(:,13) = aa(lid(4),:);
HIanom(:,14) = ab(lid(4),:);

CCanom(:,7)  = as(lid(5),:);
CCanom(:,8)  = am(lid(5),:);
CCanom(:,9)  = al(lid(5),:);
CCanom(:,10) = af(lid(5),:);
CCanom(:,11) = ap(lid(5),:);
CCanom(:,12) = ad(lid(5),:);
CCanom(:,13) = aa(lid(5),:);
CCanom(:,14) = ab(lid(5),:);

MXanom(:,7)  = as(lid(6),:);
MXanom(:,8)  = am(lid(6),:);
MXanom(:,9)  = al(lid(6),:);
MXanom(:,10) = af(lid(6),:);
MXanom(:,11) = ap(lid(6),:);
MXanom(:,12) = ad(lid(6),:);
MXanom(:,13) = aa(lid(6),:);
MXanom(:,14) = ab(lid(6),:);

SEanom(:,7)  = as(lid(7),:);
SEanom(:,8)  = am(lid(7),:);
SEanom(:,9)  = al(lid(7),:);
SEanom(:,10) = af(lid(7),:);
SEanom(:,11) = ap(lid(7),:);
SEanom(:,12) = ad(lid(7),:);
SEanom(:,13) = aa(lid(7),:);
SEanom(:,14) = ab(lid(7),:);

NEanom(:,7)  = as(lid(8),:);
NEanom(:,8)  = am(lid(8),:);
NEanom(:,9)  = al(lid(8),:);
NEanom(:,10) = af(lid(8),:);
NEanom(:,11) = ap(lid(8),:);
NEanom(:,12) = ad(lid(8),:);
NEanom(:,13) = aa(lid(8),:);
NEanom(:,14) = ab(lid(8),:);

%% table to save
Tch = array2table(CHanom,'VariableNames',['Year',iname,gname]);
Tbs = array2table(BSanom,'VariableNames',['Year',iname,gname]);
Tak = array2table(AKanom,'VariableNames',['Year',iname,gname]);
Thi = array2table(HIanom,'VariableNames',['Year',iname,gname]);
Tcc = array2table(CCanom,'VariableNames',['Year',iname,gname]);
Tmx = array2table(MXanom,'VariableNames',['Year',iname,gname]);
Tse = array2table(SEanom,'VariableNames',['Year',iname,gname]);
Tne = array2table(NEanom,'VariableNames',['Year',iname,gname]);

writetable(Tch,[fpath 'CHK_fosi_inputs_fish_annual_means_anomalies.csv'],'WriteRowNames',false);
writetable(Tbs,[fpath 'EBS_fosi_inputs_fish_annual_means_anomalies.csv'],'WriteRowNames',false);
writetable(Tak,[fpath 'GAK_fosi_inputs_fish_annual_means_anomalies.csv'],'WriteRowNames',false);
writetable(Thi,[fpath 'HI_fosi_inputs_fish_annual_means_anomalies.csv'],'WriteRowNames',false);
writetable(Tcc,[fpath 'CCE_fosi_inputs_fish_annual_means_anomalies.csv'],'WriteRowNames',false);
writetable(Tmx,[fpath 'GMX_fosi_inputs_fish_annual_means_anomalies.csv'],'WriteRowNames',false);
writetable(Tse,[fpath 'SE_fosi_inputs_fish_annual_means_anomalies.csv'],'WriteRowNames',false);
writetable(Tne,[fpath 'NE_fosi_inputs_fish_annual_means_anomalies.csv'],'WriteRowNames',false);



