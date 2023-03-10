% Plot max coeffs of climate-driver or climate-fish regressions
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

load([spath,'LMEs_regressDJF_climate_div2SD_maxcoefs.mat'])

%%  ---------------------------------------------------------
cnam = {'coef','p','lag','iclimate','climate'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

% 1= AMO, 2=AO, 3=NAO, 4=Nino34, 5=PDO
mtype = {'o','s','^','v','d','p'}; %for lag
mcol = [238/255 102/255 119/255;... %red
    34/255 136/255 51/255;...   %green
    170/255 51/255 119/255;...  %purple
    51/255 187/255 238/255;...  %cyan
    0/255 68/255 136/255];    %blue
    

%colorblind friendly
% cb=[34/255 136/255 51/255;...   %green
%     153/255 153/255 51/255;...  %olive
%     51/255 187/255 238/255;...  %cyan
%     0/255 68/255 136/255;...    %blue
%     238/255 102/255 119/255;... %red
%     170/255 51/255 119/255;...  %purple
%     0 0 0;...                   %black
%     0.25 0.25 0.25;...             %dk grey
%     0.50 0.50 0.50;...             % grey
%     0.75 0.75 0.75];               %lt grey

%%
figure(1)
plot(0:67,zeros(68,1),'--k'); hold on;
for i=1:length(lid)
    L=lid(i);
    if (LAtab(i,2) <= 0.05)
        plot(L,LAtab(i,1),mtype{(LAtab(i,3))+1},'Color',mcol(LAtab(i,4),:),'MarkerFaceColor',mcol(LAtab(i,4),:));
        hold on
    else
        plot(L,LAtab(i,1),mtype{(LAtab(i,3))+1},'Color',mcol(LAtab(i,4),:));
        hold on
    end
end
xlim([0 67])
ylim([-0.8 0.8])
xlabel('LME')
ylabel('Coeff')
%legend(tanom2(:,1))



