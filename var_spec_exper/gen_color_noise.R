
rm(list=ls())

library(RobPer)

setwd("/Users/cpetrik/Dropbox/Princeton/MAPP-METF/NCAR3/var_spec_exper/")

x0 <- TK95(N=1800, alpha = 0)
x1 <- TK95(N=1800, alpha = 1)
x2 <- TK95(N=1800, alpha = 2)

tt1 <- seq(along=x0)

# Plot Lomb-Scargle periodogram with log-axes:
y <- x0
tt <- tt1
# Show time series:
plot(tt,y, type="l", main="Power Law Noise", xlab="t", ylab="y")
# Plot Fourier periodogram with log-axes:
temp <- spectrum(y, plot=FALSE)
plot(log(temp$freq), log(temp$spec), main="log-log-Fourier periodogram",
     xlab="log(frequency)", ylab="log(periodogram)")
# A line with slope -alpha for comparison
abline(a=7, b=0, col="red")
text(-2, 12, expression(alpha==1.5), col="red")


y <- x1
# Show time series:
plot(tt,y, type="l", main="Power Law Noise", xlab="t", ylab="y")
# Plot Fourier periodogram with log-axes:
temp <- spectrum(y, plot=FALSE)
plot(log(temp$freq), log(temp$spec), main="log-log-Fourier periodogram",
     xlab="log(frequency)", ylab="log(periodogram)")
# A line with slope -alpha for comparison
abline(a=7, b=-1, col="red")
text(-2, 12, expression(alpha==1), col="red")


y <- x2
# Show time series:
plot(tt,y, type="l", main="Power Law Noise", xlab="t", ylab="y")
# Plot Fourier periodogram with log-axes:
temp <- spectrum(y, plot=FALSE)
plot(log(temp$freq), log(temp$spec), main="log-log-Fourier periodogram",
     xlab="log(frequency)", ylab="log(periodogram)")
# A line with slope -alpha for comparison
abline(a=7, b=-2, col="red")
text(-2, 12, expression(alpha==2), col="red")


### Create 11 time series of random noise slope 0 to -1.5
al <- seq(0,1.5,by=0.15)
xall <- matrix(NA, nrow = 11, ncol = 1800)
for(i in 1:11)
{
  x <- TK95(N=1800, alpha = al[i])
  xall[i,] <- x
}

write.table(xall,"color_noise_rand11.csv",sep=",",row.names=F)




### EXAMPLE FROM DOCUMENTATION
set.seed(31)
# Generate power law noise with exponent alpha=1.5:
y <- TK95(N=2000, alpha=1.5)
tt <- seq(along=y)

# Show time series:
plot(tt,y, type="l", main="Power Law Noise", xlab="t", ylab="y")

# Plot Fourier periodogram with log-axes:
temp <- spectrum(y, plot=FALSE)
plot(log(temp$freq), log(temp$spec), main="log-log-Fourier periodogram",
     xlab="log(frequency)", ylab="log(periodogram)")

# A line with slope -alpha for comparison
abline(a=8, b=-1.5, col="red")
text(-2, 12, expression(alpha==1.5), col="red")

