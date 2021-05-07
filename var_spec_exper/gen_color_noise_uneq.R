
rm(list=ls())

library(RobPer)

setwd("/Users/cpetrik/Dropbox/Princeton/MAPP-METF/NCAR3/var_spec_exper/")

#tt <- 1:365
tt <- sampler(ttype="unif", ps=1, ncycles=365, npoints=365)

x0 <- TK95_uneq(tt, alpha = 0)
x1 <- TK95_uneq(tt, alpha = 1)
x2 <- TK95_uneq(tt, alpha = 2)

plot(tt,x0)
plot(tt,x1)
plot(tt,x2)

# Plot Lomb-Scargle periodogram with log-axes:
y <- x0
temp <- RobPer(cbind(tt,y,1), weighting=FALSE, model="sine", regression="L2",
               periods=2000/seq(2, 1000, 2))
plot(log(seq(2, 1000, 2)/2000), log(temp), main="log-log-Fourier periodogram",
     xlab="log(frequency)", ylab="log(periodogram)")
abline(a=-6, b=0, col="red")

y <- x1
temp <- RobPer(cbind(tt,y,1), weighting=FALSE, model="sine", regression="L2",
               periods=2000/seq(2, 1000, 2))
plot(log(seq(2, 1000, 2)/2000), log(temp), main="log-log-Fourier periodogram",
     xlab="log(frequency)", ylab="log(periodogram)")
abline(a=-8, b=-1, col="red")

y <- x2
temp <- RobPer(cbind(tt,y,1), weighting=FALSE, model="sine", regression="L2",
               periods=2000/seq(2, 1000, 2))
plot(log(seq(2, 1000, 2)/2000), log(temp), main="log-log-Fourier periodogram",
     xlab="log(frequency)", ylab="log(periodogram)")
abline(a=-10, b=-2, col="red")


### Create 11 time series of random noise slope 0 to -1.5
#doesn't work
#xall <- TK95_uneq(tt, alpha = seq(0,1.5,by=0.15))
al <- seq(0,1.5,by=0.15)
xall <- matrix(NA, nrow = 11, ncol = 1800)
for(i in 1:11)
{
  tt <- sampler(ttype="unif", ps=1, ncycles=150, npoints=1800)
  x <- TK95_uneq(tt, alpha = al[i])
  xall[i,] <- x
}

write.table(xall,"color_noise_rand11.csv",sep=",",row.names=F)




### EXAMPLE FROM DOCUMENTATION
# Compare with example in TK95 to see that the power law is much more clear in
# equally sampled data! - draws random sample
set.seed(31)
tt <- sampler(ttype="unif", ps=1, ncycles=2000, npoints=2000)
# Generate power law noise with exponent alpha=1.5:
y <- TK95_uneq(tt, alpha=1.5)

# Show time series:
plot(tt,y, type="l", main="Irregular Power Law Noise", xlab="t", ylab="y")

# Plot Lomb-Scargle periodogram with log-axes:
temp <- RobPer(cbind(tt,y,1), weighting=FALSE, model="sine", regression="L2",
               periods=2000/seq(2, 1000, 2))
plot(log(seq(2, 1000, 2)/2000), log(temp), main="log-log-Fourier periodogram",
     xlab="log(frequency)", ylab="log(periodogram)")
title(main= "Power Law not so obvious", cex.main=0.8, line=0.5)

# A line with slope -alpha for comparison
abline(a=-10, b=-1.5, col="red")
text(-5, -1.5, expression(alpha==1.5), col="red")
