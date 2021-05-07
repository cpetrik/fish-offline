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
