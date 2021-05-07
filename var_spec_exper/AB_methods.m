%Set the color of noise to simulate
    alpha=1.0; %Noise color
    %Generate red noise of known color
    %I copy function ftk below
    rednoise=ftk(iter,alpha);
    

    %Scale to have a mean of zero and std=1;
    rednoise_normdist=(rednoise-mean(rednoise))./std(rednoise);
    

    %Interpolate onto an annual cycle with steps equivalent to the number of model time steps in a year.
    interpx=linspace(0,12,stepsperyear);
    

   %Interpolate PAR and nitrate data onto timesteps
   par_clim_spline=spline(0:1:12,[par_clim_month(12,2);par_clim_month(:,2)],interpx);
   scalefact_sqrt_vary_spline=spline(0:1:12,[scalefact_sqrt_vary(12);scalefact_sqrt_vary],interpx);  %Make sure variance is seasonal
   maxres_sqrt_vary_spline=spline(0:1:12,[maxres_sqrt_vary(12);maxres_sqrt_vary],interpx);
   

    %Scale PAR to be from 0-1
   par_clim_spline_scaled=(par_clim_spline-min(par_clim_spline))./(max(par_clim_spline)-min(par_clim_spline));
   

   %Scale PAR back to nitrate
   par_clim_spline_Nscaled=par_clim_spline_scaled*(max(nitrate_sqrt_clim)-min(nitrate_sqrt_clim))+min(nitrate_sqrt_clim); 
  

   %Repeat this one year over all model run years
   par_clim_spline_Nscaled_repmat=repmat(par_clim_spline_Nscaled',[simulationyears,1]);
   maxres_sqrt_vary_spline_repmat=repmat(maxres_sqrt_vary_spline',[simulationyears,1]);
   maxres_sqrt_vary_spline_repmat=repmat(maxres_sqrt_vary_spline',[simulationyears,1]);
    

   %Scale by the scale factor. This produces unrealistically high K values
    %for log and ln transforms at high and low alpha. Sqrt looks OK.
    rednoise_scaled_sqrt=rednoise_normdist.*scalefact_sqrt;
    

    %Combine repeating climatology and noise, and square
    nitrate_backtrans_sqrt=(rednoise_scaled_sqrt+par_clim_spline_Nscaled_repmat').^2;
   

    K=nitrate_backtrans_sqrt';
    K=K+0.1; %Offset so it wasn’t ever negative

%% Here is the function ftk:


%Create power law noise of a known exponent (a) and for a time series of
%length n.
%
% From: Timmer, J.; Koenig, M. 1995. On generating power law noise. Astronomy and Astrophysics 300, 707
%
% pink noise        a = 1
% red noise         a = 1.5
 

 

function x=ftk(n,a)
 

w=linspace(0,pi,n/2+1);
w= w(2:end-1);
p=1./w.^a;
  

rl=sqrt(0.5*p).*randn(1,n/2-1);
im=sqrt(0.5*p).*randn(1,n/2-1);
  

x=real(ifft(complex([randn(1),rl,randn(1),rl(end-1:1)], ...
    [0,im,0,-im(end-1:1)]),n));
