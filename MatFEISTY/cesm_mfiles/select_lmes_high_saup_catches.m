% SAUP catch by LME to find LMEs with high catch by fn type
% To select case studies for DPLE ms
% ms says "These LMEs were chosen for each functional type as having both 
% (1) a large proportion of total biomass contributed by that functional 
% type, and (2) the LME mean catch of that functional type was greater than 
% the median catch of that functional type across all LMEs.” 

clear 
close all

%% SAUP
spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/SAUP/';

% use weighted catches
% use top 10 catch yrs
load([spath 'SAUP_LME_Catch_top10_Stock.mat']); %units: MT/km2
 
%%
Dmed = median(Dlme_mcatch10,'omitnan');
Fmed = median(Flme_mcatch10,'omitnan');
Pmed = median(Plme_mcatch10,'omitnan');

[did,~] = find(Dlme_mcatch10 > Dmed);
[fid,~] = find(Flme_mcatch10 > Fmed);
[pid,~] = find(Plme_mcatch10 > Pmed);

%% Table
Dlme(:,1) = did;
Dlme(:,2) = Dlme_mcatch10(did);

Flme(:,1) = fid;
Flme(:,2) = Flme_mcatch10(fid);

Plme(:,1) = pid;
Plme(:,2) = Plme_mcatch10(pid);

Fstat = array2table(Flme,'VariableNames',{'LME','Catch'});
writetable(Fstat,[spath 'LME_SAUP_high_catch_F.csv'],'Delimiter',',',...
    'WriteRowNames',false)

Dstat = array2table(Dlme,'VariableNames',{'LME','Catch'});
writetable(Dstat,[spath 'LME_SAUP_high_catch_D.csv'],'Delimiter',',',...
    'WriteRowNames',false)

Pstat = array2table(Plme,'VariableNames',{'LME','Catch'});
writetable(Pstat,[spath 'LME_SAUP_high_catch_P.csv'],'Delimiter',',',...
    'WriteRowNames',false)

save([sPmed ...
    path 'LME_SAUP_stats_high_catch_types.mat'])

