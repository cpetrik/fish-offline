%Filename: powerspectra_TS.m
%
%Creator: Andrew Barton (abarton@princeton.edu)
%
%Date: September 29, 2015
%
%Purpose: Open and read annual anomaly time series and calculate power
%spectral slopes. I prefer it to the periodogram function on Matlab
%because of the Blakman and Tukey cosine bell taper, which tends to give
%smoother slopes for the power spectra.
%
%Input: 
%x=anomaly time series. The long-term trend should be removed.
%f = sampling interval in months. Use this to convert the frequency axis
%to readable scale. 12 for monthly data, 1 for annual data. 
%If fignum~=0, do a plot at number specified.
%
%Output: 
%f=frequency, in cycles per entire time series length.
%p=power
%b=slope
%int=intercept
%
%Notes: Updated by C. Petrik 3/22/21 to use TheilSen.m instead of 
% Theil_Sen_Regress.m because much faster
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [freq,p,b,int]=powerspectra_TS(x,f,fignum)

%Total number of samples
n=length(x);

%Tapering
ptaper=0.1;
m=floor(n*ptaper);
w=0.5*( 1 - cos( pi * ( 1 : 2 : 2 * m - 1) / (2 * m) ) ); %Blackman and Tukey (1959) cosine bell taper
x=[w ones( 1, n - 2 * m ) fliplr( w )]' .* x;

%Calculate FFT and power
xdft=fft(x);
xdft=xdft(1:n/2);  %Take only the positive half
p=(1/n).*abs(xdft).^2;  %Normalize by the number of samples, square to get power
p(1)=0.5*(p(2)+p(end));  

freq=((0:n/2-1)/n)';
p=p(2:end);
freq=freq(2:end);

%Get slope. Theil-sen slope is probably be more robust, but regress is
%quicker.
%b=regress(log10(p), [ones(length(w),1), log10(w)] );
%b=b(2);
%Don't include the longest period in the fit
fitnum=1;
% b=Theil_Sen_Regress_wduplicates(log10(freq(fitnum:end)),log10(p(fitnum:end)));
% b=Theil_Sen_Regress(log10(freq(fitnum:end)),log10(p(fitnum:end)));
% int=nanmedian(log10(p(fitnum:end))-b.*log10(freq(fitnum:end)));
[b,int] = TheilSen([log10(freq(fitnum:end)),log10(p(fitnum:end))]);

%Plot out spectra.
if (fignum~=0)
    figure(fignum)
    plot(log10(freq),log10(p))
    set( gca, 'XTick',log10(1./([ 500 200 100 50 20 10 5 2 1 6/12 2/12].*f)), ...
        'XTickLabel', { '500a'; '200a';'100a'; '50a'; '20a'; '10a'; '5a';
         '2a'; '1a'; '6m'; '2m' } )
    hold on;
    plot(log10(freq),int+b.*log10(freq),'r')
    title('Periodogram Using FFT')
    xlabel('Period')
    ylabel('Power/Frequency (B/Hz)')
    grid on;
    hold off;
    
    
end