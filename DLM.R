SP1990.2015 <- read.csv("C:/Users/asdf/Desktop/SP1990-2015.csv", header=FALSE)
sp <- ts(SP1990.2015$V2, frequency = 52, start = 1990)
Lsp <- log(sp)
plot(sp,ylab = "sp500")
plot(Lsp, xlab = "", ylab = "log sp500")
acf(sp)
acf(Lsp)
SP1990 <- read.csv("C:/Users/asdf/Desktop/SP1990.csv", header=FALSE)
SP2016 <- read.csv("C:/Users/asdf/Desktop/SP2016.csv", header=FALSE)

sp<-num_earthquakes
sp<-Magnitudes
build<-function(u){
  rw<-dlmModPoly(order = 1, dV = u[1],dW = u[2])
  return(rw)
}
init <- c(1e-7,1)

level0<-sp[1]
slope0<-mean(diff(sp))
build<-function(u){
  trend<-dlmModPoly(dV=1e-7, dW=exp(u[1:2]),m0=c(level0, slope0),C0=2*diag(2))
  gap<-dlmModARMA(ar=u[4:5],sigma2=exp(u[3]))
return(trend+gap)}
init<-c(-3,-1,-3,.4,.4)



outmle<-dlmMLE(sp,init,build)
outmle$par
SpPoly <-dlmModPoly(order = 1, dV= 60.51726,dW =557.41563)
SpFilt <- dlmFilter(sp, SpPoly)
#a little more comprehensive summary of the filtering results
str(SpFilt, 1)
n <- length(sp)
plot(sp, type='o', col = c("darkgrey"),  xlab = "", ylab = "", lwd = 4)
lines(dropFirst(SpFilt$m), lty = "longdash", lwd = 4)
SpSmooth <- dlmSmooth(SpFilt)
attach(SpSmooth)
n<-length(sp)

###these two should be the same, as they are both \theta_n| y_{1:n}
drop(dlmSvd2var(U.S[[n + 1]], D.S[n + 1,]))
drop(dlmSvd2var(U.C[[n + 1]], D.C[n + 1,]))

##these two are different: the first one is Var(theta_t|y_{1:n}}
##and the second is Var(theta_t|y_{1:t}} 
drop(dlmSvd2var(U.S[[n / 2 + 1]], D.S[n / 2 + 1,]))
drop(dlmSvd2var(U.C[[n / 2 + 1]], D.C[n / 2 + 1,]))

smoothed_var_record=rep(0, n+1)
for(i in 1:(n+1)){
  smoothed_var_record[i]=dlmSvd2var(U.S[[i]], D.S[i, ])
}

##C_0
smoothed_var_record[1]
##plot the rest
plot(smoothed_var_record[2:(n+1)],type='b')

#####Smoothing with W=755 (small signal) and W=7550 (large signal)
#NileSmooth1 <- dlmSmooth(NileFilt1)
#NileSmooth2 <- dlmSmooth(NileFilt2)

#########plot the filtering results
#par(mfrow=c(1,2))
#attach(NileSmooth1)
hwid <- qnorm(0.025, lower = FALSE) *  sqrt(unlist(dlmSvd2var(U.S, D.S)))
smooth <- cbind(s, as.vector(s) + hwid %o% c(-1, 1))

plot(dropFirst(smooth), plot.type = "s", type = "l",
     lty = c(1, 5, 5), ylab = "Level", xlab = "",
     ylim = range(sp))
lines(sp, type = "o", col = "darkgrey")
legend("bottomleft", col = c("darkgrey", rep("black", 2)),
       lty = c(1, 1, 5), pch = c(1, NA, NA), bty = "n",
       legend = c("data", "smoothed level",
                  "95% probability limits"))

a <- window(cbind(sp, SpFilt$f),
            start = 1966, end = 2016)
plot(a[, 1], type = 'o', col = "darkgrey",
     xlab = "", ylab = "")
lines(a[, 2], lty = "longdash")
leg <- c("data", "one-step-ahead forecast")
legend("topleft", legend = leg, col = c("darkgrey", "black"),
       lty = c("solid", "longdash"),
       pch = c(1, NA), bty = "n")
par(mfrow=c(1,1))
