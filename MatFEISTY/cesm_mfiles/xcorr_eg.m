n=1:15;
x=0.84.^n;
y=circshift(x,5);
[c,lags] = xcorr(y,x);
stem(lags,c)