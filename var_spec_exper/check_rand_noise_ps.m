% Calc powerspectra of rand generated noise of given color
% Used TK95 in R package RobPer

clear all
close all

%% 
load('color_noise_rand11.mat')

[ns,nt] = size(noise);

noise = noise';
ps = NaN*ones(ns,1);

%% calc powerspec
for i = 1:ns 
    i
 
    xi = noise(:, i);
    inan = ~isnan(xi);
    %R = log(xi(inan)); %Already log-transformed
    R = (xi(inan));
    T = length(R);
    t = (1:T)';
    b = Theil_Sen_Regress(t,R);
    int = nanmedian(R-b.*t);
    tH = b*t + int;
    dR = R - tH;
    [freq1,p1,b1,int1] = powerspectra(dR,12,0);
    ps(i) = b1;
    clear R T t b int tH dR freq1 p1 b1 int1
    
end
