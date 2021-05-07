% Use seasonal climatologies of phys & BGC to create
% experimental time-series
% sqrt-transform abund
% No transform environ (temp)

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CORE-forced/';

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Climatologies
load([fpath 'cobalt_biome_core_climatol_1950_2007.mat']);

%% Scale factors
%1: raw, 2:anom, 3:log10
load([fpath 'cobalt_biome_core_anom_1950_2007.mat'],...
    'tp_sc');

%% 5. Interpolate onto an annual cycle with steps equivalent to the number of
%model time steps in a year.
stepsperyear = 365;
interpx=linspace(0,12,stepsperyear);

%Interpolate variables and scale factors onto timesteps
tp_clim_spl = spline(0:1:12,[tp_clim(:,12),tp_clim],interpx);

%Make sure variance is seasonal
tp_sc_spl = spline(0:1:12,[tp_sc(:,12),tp_sc],interpx);

%% Repeat this one year over all model run years
simulationyears = 150;
tp_clim_spl_rep=repmat(tp_clim_spl',[simulationyears,1]);

tp_sc_spl_rep=repmat(tp_sc_spl',[simulationyears,1]);

%% Generate red noise of known color
iter = stepsperyear * simulationyears;

%Set the color of noise to simulate
al = 0:0.15:1.5;

noise = NaN*ones(length(al)*20,iter);
k=0;
for i=1:length(al)
    alpha = al(i);
    %20 random draws of each alpha (100 too big for storage)
    for j=1:20
        k=k+1;
        rednoise = ftk(iter,alpha);
        %Scale to have a mean of zero and std=1;
        rednoise_normdist=(rednoise-mean(rednoise))./std(rednoise);
        
        noise(k,:) = rednoise_normdist;
    end
end

clear rednoise_normdist rednoise

%% 7. Scale by the scale factor
%one for each variable and biome
nal = size(noise,1);
tot = nal * 8;
st = 1:nal:tot;
en = nal:nal:tot;
noise_sc_tp = NaN*ones(tot,iter);

for b=1:8
    scaled = noise .* tp_sc_spl_rep(:,b)';
    noise_sc_tp(st(b):en(b),:) = scaled;
end

%% Combine repeating climatology and noise, and square
tp_clim_spl_rep_backtrans=(noise_sc_tp + repelem(tp_clim_spl_rep',220,1));

save([fpath 'cobalt_biome_core_noise_tp.mat'],'tp_clim_spl_rep_backtrans');%, '-v7.3');

%% Check abund never negative
hist(tp_clim_spl_rep_backtrans(:))
%K=K+0.1; %Offset so it wasn’t ever negative

