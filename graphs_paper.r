
rm(list=ls())

#### ------- Chi-square pdf
dChisqScore <-function(x,a,b,lambda){
  # pdf of the log score based on Gaussian additive model
  # i.e. s(f,Y) = a + b * Chi-square
  # with Chi-square is 1 dof non central chi-square with  ncp = lambda
  input<-(x-a)/b
  out<-dchisq(x=input, df=1, ncp = lambda, log = FALSE)
  out<-as.vector(out)/b
  return(out)	}

#### ------- Figure 1 - PIT histogram
# PIT histograms
mu0 <- 0
sig0 <- 2
c <- 1
mu1 <- mu0
sig1 <- sig0
mu2 <- mu0 + 1
sig2 <- sig0 - .3
x0 <- rnorm(n=10000,mean=0,sd=1)
x <- mu0 + sig0*x0
y = mu0 + sqrt(sig0^2+c^2)*x0

uY1 = pnorm(q=y,mean=mu1,sd=sig1)
uY2 = pnorm(q=y,mean=mu2,sd=sig2)
uX1 = pnorm(q=x,mean=mu1,sd=sig1)
uX2 = pnorm(q=x,mean=mu2,sd=sig2)
#
layout(matrix(1:2,1,2,byrow=TRUE))
xb <- seq(0,1,length.out=25)
par(mar=c(4.5, 5.5, 4.5, 1.5))
hist(uX1,probability=TRUE,breaks=xb,ylim=c(0,1.7),xlab='Probability Integral Transform',ylab='Frequency',main='Perfect verification data',col=rgb(0.1,0.1,0.1,0.5),cex=2,cex.lab=2,cex.axis=1.4,xaxs="i",cex.main=2)
hist(uX2,probability=TRUE,breaks=xb,add=TRUE,col=rgb(0.8,0.8,0.8,0.5))
legend(x='topright',c('Ideal forecast','Imperfect forecast'),col=c(rgb(0.1,0.1,0.1,0.5),rgb(0.8,0.8,0.8,0.5)),pch=15,cex=1.7,bty='n')
abline(h=1,lwd=2)
#
hist(uY1,probability=TRUE,breaks=xb,ylim=c(0,1.6),xlab='Probability Integral Transform',ylab='',main='Imperfect verification data',col=rgb(0.1,0.1,0.1,0.5),cex=2,cex.lab=2,cex.axis=1.4,xaxs="i",cex.main=2)
hist(uY2,probability=TRUE,breaks=xb,add=TRUE,col=rgb(0.8,0.8,0.8,0.5))
abline(h=1,lwd=2)


#### ------- Figure 2 - Time series of wind speed data
# station 8 from (Bessac et al., 2018), x=(nwp25k, nwp5k, obs) data starts jan, 2 at 6am
# nwp25k: NWP outputs at 25km, nwp5k: NWP outputs at 5km, obs: Observation
load('station8_C1_Jan12_ws.RData')
sig1 <- 0.9425318  # NWP estimated uncertainty from (Bessac et al., 2018)
lab <- unlist(lapply(21:30,function(i) paste('Jan',i)))
tw <- 456:504
h <- 1:length(tw)
par(mar=c(4.5, 5, 3, 2))
plot(h,x[tw,3],typ='l',lty=3,ylim=range(x[tw,],8.5),lwd=1.5,xaxs='i',xlab='Time (hour)',ylab='Wind speed (m/s)',xaxt='n',cex.lab=1.7,cex=1.3,cex.axis=1.5)
polygon(c(h,rev(h)),c(c(x[tw,1]+1*sig1),rev(x[tw,1]-1*sig1)),col="grey63",border=NA)
polygon(c(h,rev(h)),c(c(x[tw,3]+1*sqrt(0.25)),rev(x[tw,3]-1*sqrt(0.25))),col="grey83",border=NA)
lines(h,x[tw,3],lwd=1.5,lty=1)
lines(h,x[tw,1],lwd=1.5,lty=2)
axis(1,at=seq(12,length(tw),by=24),labels=lab[1:length(seq(12,length(tw),by=24))],cex.axis=1.3)
legend(x='topleft',c('Observation','Prediction model','Measurement uncertainty','Prediction uncertainty'),lty=c(1,2,NA,NA),col=c(1,1,'grey83','grey63'),pch=c(NA,NA,15,15),lwd=2,bty='n',cex=1.4)



#### ------- Figure 3 - Log-score distributions for additive model
# true hidden state
mu0 <- 0
sig0 <- 2
x0 <- rnorm(n=10000,mean=0,sd=1)
x <- mu0 + sig0*x0

## --- forecast 1 - ideal forecast
mup1 <- mu0
sp1 <- sig0

## --- forecast 2 - biased forecast
mup2 <- mu0 + 1
sp2 <- sig0 + 1

## ------------------ figure 3 - graphs of densities 
sigma = sp1
mu = mup1
#
xl <- 3.1
x <- seq(from=1.5, to=xl, length=90) # forecast 0
xF <- seq(from=1.5, to=xl, length=90) # forecast 0
cst <- log(sigma)+(log(2*pi)/2)
TrueMean0 <- cst+((mu0-mu)^2+sig0^2)/(2*sigma^2)
#
seq_c = c(sqrt(.5),sqrt(1),sqrt(3))
rNew = NULL
r0 = NULL
for (i in 1:length(seq_c)){
  omega = seq_c[i]
  p0<- sig0^2/(sig0^2+omega^2)
  p0<-p0^2
  # OBSERVED LOG SCORE PDF
  a0 <- cst 
  b0 <- (sig0^2+omega^2)/(2*sigma^2)
  lambda0 <- (mu0-mu)^2/(sig0^2+omega^2)
  r0[[i]] <- dChisqScore(x,a0,b0,lambda0)
  # OUR NEW LOG SCORE PDF
  aNew <- cst + (1/(2*sigma^2)) *omega^2 * sqrt(p0)
  bNew <- b0*p0
  lambdaNew <- lambda0/p0
  rNew[[i]] <- dChisqScore(x,aNew,bNew,lambdaNew) }
#####
sigma = sp2
mu = mup2
cst <- log(sigma)+(log(2*pi)/2)
TrueMean1 <- cst+((mu0-mu)^2+sig0^2)/(2*sigma^2)
for (i in 1:length(seq_c)){
  omega = seq_c[i]
  p0<- sig0^2/(sig0^2+omega^2)
  p0<-p0^2
  # OBSERVED LOG SCORE PDF
  a0 <- cst 
  b0 <- (sig0^2+omega^2)/(2*sigma^2)
  lambda0 <- (mu0-mu)^2/(sig0^2+omega^2)
  r0[[(i+3)]] <- dChisqScore(x,a0,b0,lambda0) 
  # OUR NEW LOG SCORE PDF
  aNew <- cst + (1/(2*sigma^2)) *omega^2 * sqrt(p0)
  bNew <- b0*p0
  lambdaNew <- lambda0/p0
  rNew[[(i+3)]] <- dChisqScore(x,aNew,bNew,lambdaNew)}


# left panel: ideal forecast ; right panel: imperfect forecast
op <- par(mar = c(5,5,4,1) + 0.1)
layout(matrix(1:2,1,2))
plot(x,r0[[2]],type="l",xlab="Score value",ylab="PDF of the log-score",col="green",xlim=range(x,xF,xl),ylim=range(rNew,r0),cex=1.4,cex.lab=1.4,lwd=1.5,xaxs='i')
grid(lwd=2)
#
lines(x,rNew[[1]],col="red",lwd=1.5,lty=1)
polygon(c(1,x,xl),c(0,rNew[[1]],0),col=rgb(red=0.4,green=.0,blue=0.0,alpha=.2),border=F)
lines(x,rNew[[2]],col="red",lwd=1.5,lty=2)
polygon(c(1,x,xl),c(0,rNew[[2]],0),col=rgb(red=0.4,green=.0,blue=0.0,alpha=.2),border=F)
lines(x,rNew[[3]],col="red",lwd=1.5,lty=3)
polygon(c(1,x,xl),c(0,rNew[[3]],0),col=rgb(red=0.4,green=.0,blue=0.0,alpha=.2),border=F)
#
lines(x,r0[[1]],col="green",lwd=1.5,lty=1)
polygon(c(1,x,xl),c(0,r0[[1]],0),col=rgb(red=0.,green=.2,blue=0.0,alpha=.1),border=F)
lines(x,r0[[2]],col="green",lwd=1.5,lty=2)
polygon(c(1,x,xl),c(0,r0[[2]],0),col=rgb(red=0.,green=.2,blue=0.0,alpha=.1),border=F)
lines(x,r0[[3]],col="green",lwd=1.5,lty=3)
polygon(c(1,x,xl),c(0,r0[[3]],0),col=rgb(red=0.,green=.2,blue=0.0,alpha=.1),border=F)
#
abline(v=TrueMean0,col="orange", lwd=2)
legend(x='topright',c(expression(s[0](f[0],Y)),expression(s[v](f[0],Y)),expression(omega^2~'= 0.5'),expression(omega^2~'= 1'),expression(omega^2~'= 3')),lty=c(1,1,1,2,3),col=c(3,2,1,1,1),cex=1.6,lwd=1.5,bty='n')
axis(side=3,at=TrueMean0,label=expression('E('~s[0]~'('~f[0]~',X))'),tick=TRUE,cex.axis=1.5,line=0)
####
plot(x,r0[[5]],type="l",xlab="Score value",ylab="",col="green",xlim=range(x,xF,xl),ylim=range(rNew,r0),cex=1.4,cex.lab=1.4,lwd=1.5,xaxs='i')
grid(lwd=2)
#
lines(x,rNew[[4]],col="red",lwd=1.5,lty=1)
polygon(c(1,x,xl),c(0,rNew[[4]],0),col=rgb(red=0.4,green=.0,blue=0.0,alpha=.2),border=F)
lines(x,rNew[[5]],col="red",lwd=1.5,lty=2)
polygon(c(1,x,xl),c(0,rNew[[5]],0),col=rgb(red=0.4,green=.0,blue=0.0,alpha=.2),border=F)
lines(x,rNew[[6]],col="red",lwd=1.5,lty=3)
polygon(c(1,x,xl),c(0,rNew[[6]],0),col=rgb(red=0.4,green=.0,blue=0.0,alpha=.2),border=F)
#
lines(x,r0[[4]],col="green",lwd=1.5,lty=1)
polygon(c(1,x,xl),c(0,r0[[4]],0),col=rgb(red=0.,green=.2,blue=0.0,alpha=.1),border=F)
lines(x,r0[[5]],col="green",lwd=1.5,lty=2)
polygon(c(1,x,xl),c(0,r0[[5]],0),col=rgb(red=0.,green=.2,blue=0.0,alpha=.1),border=F)
lines(x,r0[[6]],col="green",lwd=1.5,lty=3)
polygon(c(1,x,xl),c(0,r0[[6]],0),col=rgb(red=0.,green=.2,blue=0.0,alpha=.1),border=F)
#
abline(v=TrueMean1,col="orange", lwd=2)
legend(x='topright',c(expression(s[0](f,Y)),expression(s[v](f,Y)),expression(omega^2~'= 0.5'),expression(omega^2~'= 1'),expression(omega^2~'= 3')),lty=c(1,1,1,2,3),col=c(3,2,1,1,1),cex=1.6,lwd=1.5,bty='n')
axis(side=3,at=TrueMean1,label=expression('E('~s[0]~'(f,X))'),tick=TRUE,cex.axis=1.5,line=0)




#### ------- Figure 4 - CRPS distribution for multiplicative model
### ------ CRPS with multiple a 
library('invgamma')

ns = 1000
## truth signal 
alpha0 = 7 #6.5
beta0 = 2 #.7
x0 = rgamma(n=ns,shape=alpha0,rate=beta0)

## ideal forecast
alpha = alpha0
beta = beta0
fx = dgamma(x0,shape=alpha,rate=beta)
Fbarx = pgamma(x0,shape=alpha,rate=beta,lower.tail=FALSE)
c00 = x0 + 2*(x0*fx/beta + (alpha/beta - x0)*Fbarx) - (alpha/beta + gamma(.5+alpha)/(beta*gamma(.5)*gamma(alpha)))
r0 = density(c00)

## imperfect forecast
alpha1 = alpha0-3 
beta1 = beta0-1
fx = dgamma(x0,shape=alpha1,rate=beta1)
Fbarx = pgamma(x0,shape=alpha1,rate=beta1,lower.tail=FALSE)
c01 = x0 + 2*(x0*fx/beta1 + (alpha1/beta1 - x0)*Fbarx) - (alpha1/beta1 + gamma(.5+alpha1)/(beta1*gamma(.5)*gamma(alpha1)))
r01 = density(c01)

seq_c = c(1,3)
rY = NULL
rNew = NULL
for (i in 1:length(seq_c)){
  c = seq_c[i]
  a = 6 + c
  b = 8
  err = rinvgamma(n=ns,shape=a,rate=b)
  y = x0*err
  fy = dgamma(y,shape=alpha,rate=beta)
  Fbary = pgamma(y,shape=alpha,rate=beta,lower.tail=FALSE)
  c0 = y + 2*(y*fy/beta + (alpha/beta - y)*Fbary) - (alpha/beta + gamma(.5+alpha)/(beta*gamma(.5)*gamma(alpha)))
  #
  fxbar = c()
  Fbarxbar = c()
  for (j in 1:ns){
    u = rgamma(n=ns,shape=(alpha0+a),rate=(beta0+(b/y[j])))
    fU = dgamma(u,shape=alpha,rate=beta)
    fxbar[j] = mean((u/beta)*fU)
    FbarU = pgamma(u,shape=alpha,rate=beta,lower.tail=FALSE)
    Fbarxbar[j] = mean((alpha/beta - u)*FbarU)	}
  cv = (alpha0+a)/(beta0+b/y) +  2*(fxbar+Fbarxbar) - (alpha/beta + gamma(.5+alpha)/(beta*gamma(.5)*gamma(alpha)))
  #
  rY[[i]] = density(c0)
  rNew[[i]] = density(cv)    }
###
### imperfect forecast
for (i in 1:length(seq_c)){
  c = seq_c[i]
  a = 6 + c
  b = 8
  err = rinvgamma(n=ns,shape=a,rate=b)
  y = x0*err
  fy = dgamma(y,shape=alpha1,rate=beta1)
  Fbary = pgamma(y,shape=alpha1,rate=beta1,lower.tail=FALSE)
  c0 = y + 2*(y*fy/beta1 + (alpha1/beta1 - y)*Fbary) - (alpha1/beta1 + gamma(.5+alpha1)/(beta1*gamma(.5)*gamma(alpha1)))
  #
  fxbar = c()
  Fbarxbar = c()
  for (j in 1:ns){
    u = rgamma(n=ns,shape=(alpha0+a),rate=(beta0+(b/y[j])))
    fU = dgamma(u,shape=alpha1,rate=beta1)
    fxbar[j] = mean((u/beta1)*fU)
    FbarU = pgamma(u,shape=alpha1,rate=beta1,lower.tail=FALSE)
    Fbarxbar[j] = mean((alpha1/beta1 - u)*FbarU)	}
  cv = (alpha0+a)/(beta0+b/y) +  2*(fxbar+Fbarxbar) - (alpha1/beta1 + gamma(.5+alpha1)/(beta1*gamma(.5)*gamma(alpha1)))
  #
  rY[[(i+2)]] = density(c0)
  rNew[[(i+2)]] = density(cv)    }


# left panel: ideal forecast ; right panel: imperfect forecast
op <- par(mar = c(5,5,4,1) + 0.1)
layout(matrix(1:2,1,2))
plot(r0$x,r0$y,type="l",xlab="Score value",ylab="PDFs of CRPS",xlim=c(0,3.5),ylim=range(r01$y,r0$y,rNew[[2]]$y,rNew[[1]]$y,rNew[[3]]$y,rNew[[4]]$y),cex=1.6,cex.lab=1.6,lwd=1.5,col="blue")
grid(lwd=2)
polygon(c(1,r0$x,r0$x[512]),c(0,r0$y,0),col=rgb(red=0.,green=.0,blue=0.4,alpha=.2),border=F)
#
lines(rNew[[1]]$x,rNew[[1]]$y,col="red",lwd=1.5,lty=1)
polygon(c(1,rNew[[1]]$x,rNew[[1]]$x[512]),c(0,rNew[[1]]$y,0),col=rgb(red=0.4,green=.0,blue=0.0,alpha=.2),border=F)
lines(rNew[[2]]$x,rNew[[2]]$y,col="red",lwd=1.5,lty=2)
polygon(c(1,rNew[[2]]$x,rNew[[2]]$x[512]),c(0,rNew[[2]]$y,0),col=rgb(red=0.4,green=.0,blue=0.0,alpha=.2),border=F)
#
lines(rY[[1]]$x,rY[[1]]$y,col="green",lwd=1.5,lty=1)
polygon(c(1,rY[[1]]$x,rY[[1]]$x[512]),c(0,rY[[1]]$y,0),col=rgb(red=0.,green=.2,blue=0.0,alpha=.1),border=F)
lines(rY[[2]]$x,rY[[2]]$y,col="green",lwd=1.5,lty=2)
polygon(c(1,rY[[2]]$x,rY[[2]]$x[512]),c(0,rY[[2]]$y,0),col=rgb(red=0.,green=.2,blue=0.0,alpha=.1),border=F)
#
abline(v=mean(c00),col="orange", lwd=2)
legend(x='topright',c(expression(c[0]~'('~f[0]~',X)'),expression(c[0]~'('~f[0]~',Y)'),expression(c[v]~'('~f[0]~',Y)'),expression(a~'= 7'),expression(a~'= 9')),lty=c(1,1,1,2,3),col=c(4,3,2,1,1,1),cex=1.3,lwd=1.5,bty='n')
axis(side=3,at=mean(c00),label=expression('E('~c[0]~'('~f[0]~',X))'),tick=TRUE,cex.axis=1.5)
##
##
plot(r01$x,r01$y,type="l",xlab="Score value",ylab="PDFs of CRPS",xlim=c(0,3.5),ylim=range(r01$y,r0$y,rNew[[2]]$y,rNew[[1]]$y,rNew[[3]]$y,rNew[[4]]$y),cex=1.6,cex.lab=1.6,lwd=1.5,col="blue")
grid(lwd=2)
polygon(c(1,r01$x,r01$x[512]),c(0,r01$y,0),col=rgb(red=0.,green=.0,blue=0.4,alpha=.2),border=F)
#
lines(rNew[[3]]$x,rNew[[3]]$y,col="red",lwd=1.5,lty=1)
polygon(c(1,rNew[[3]]$x,rNew[[3]]$x[512]),c(0,rNew[[3]]$y,0),col=rgb(red=0.4,green=.0,blue=0.0,alpha=.2),border=F)
lines(rNew[[4]]$x,rNew[[6]]$y,col="red",lwd=1.5,lty=2)
polygon(c(1,rNew[[4]]$x,rNew[[4]]$x[512]),c(0,rNew[[4]]$y,0),col=rgb(red=0.4,green=.0,blue=0.0,alpha=.2),border=F)
#
lines(rY[[3]]$x,rY[[3]]$y,col="green",lwd=1.5,lty=1)
polygon(c(1,rY[[3]]$x,rY[[3]]$x[512]),c(0,rY[[3]]$y,0),col=rgb(red=0.,green=.2,blue=0.0,alpha=.1),border=F)
lines(rY[[4]]$x,rY[[4]]$y,col="green",lwd=1.5,lty=2)
polygon(c(1,rY[[4]]$x,rY[[4]]$x[512]),c(0,rY[[4]]$y,0),col=rgb(red=0.,green=.2,blue=0.0,alpha=.1),border=F)
#
abline(v=mean(c00),col="orange", lwd=2)
legend(x='topright',c(expression(c[0]~'(f,X)'),expression(c[0]~'(f,Y)'),expression(c[v]~'(f,Y)'),expression(a~'= 7'),expression(a~'= 9')),lty=c(1,1,1,2,3),col=c(4,3,2,1,1,1),cex=1.3,lwd=1.5,bty='n')
axis(side=3,at=mean(c00),label=expression('E('~c[0]~'(f,X))'),tick=TRUE,cex.axis=1.5)



#### ------- Figure 5 - Power analysis

# compute the power for varying predictive mean
source('poweranalysisVSmean_logscore.r')

# compute the power for varying predictive standard deviation
source('poweranalysisVSsd_logscore.r')

setEPS()
postscript("figure5.eps",height = 18, width = 28)
layout(matrix(1:6,2,3,byrow = TRUE))
par(mar = c(10,10,7,2),mgp=c(4,1.8,0))
plot(mu,reject_sx_mu[[1]],typ='o',pch=20,xlab='',ylab='',main=expression('Power of hypothesis test - '~omega^2~'= 0.25'),xaxs="i",xaxt='n',cex=1.8,lwd=2,cex.main=4.5,cex.axis=3.5,ylim=c(0,1))
lines(mu,reject_sy_mu[[1]],col=2,pch=20,lwd=2,cex=1.5,typ='o')
lines(mu,reject_sc_mu[[1]],col=3,pch=20,lwd=2,cex=1.5,typ='o')
abline(v=mu0,lty=2)
axis(side=1,at=-3:3,label=c(-3,-2,-1,expression(mu[0]),1, 2,3),tick=FALSE,cex.axis=3.7)
mtext(side=1,text=expression('Predictive mean '~mu),line=7,cex=3.3)
mtext(side=2,text='Power',line=6,cex=3.3)
legend(x='bottomleft',c(expression(s[0]~'(.,X)'),expression(s[0]~'(.,Y)'),expression(s[v]~'(.,Y)')),col=c(1:4),pch=20,lwd=1.7,cex=4.3,bty='n')
#
plot(mu,reject_sx_mu[[2]],typ='o',pch=20,xlab='',ylab='',main=expression('Power of hypothesis test - '~omega^2~'=0.5'),xaxs="i",xaxt='n',cex=1.5,lwd=2,cex.main=4.5,cex.axis=3.5,ylim=c(0,1))
lines(mu,reject_sy_mu[[2]],col=2,pch=20,lwd=2,cex=1.5,typ='o')
lines(mu,reject_sc_mu[[2]],col=3,pch=20,lwd=2,cex=1.5,typ='o')
abline(v=mu0,lty=2)
axis(side=1,at=-3:3,label=c(-3,-2,-1,expression(mu[0]),1, 2,3),tick=FALSE,cex.axis=3.7)
mtext(side=1,text=expression('Predictive mean '~mu),line=7,cex=3.3)
#
plot(mu,reject_sx_mu[[3]],typ='o',pch=20,xlab='',ylab='',main=expression('Power of hypothesis test - '~omega^2~'= 1'),xaxs="i",xaxt='n',cex=1.5,lwd=2,cex.main=4.5,cex.axis=3.5,ylim=c(0,1))
lines(mu,reject_sy_mu[[3]],col=2,pch=20,lwd=2,cex=1.5,typ='o')
lines(mu,reject_sc_mu[[3]],col=3,pch=20,lwd=2,cex=1.5,typ='o')
abline(v=mu0,lty=2)
axis(side=1,at=-3:3,label=c(-3,-2,-1,expression(mu[0]),1, 2,3),tick=FALSE,cex.axis=3.7)
mtext(side=1,text=expression('Predictive mean '~mu),line=7,cex=3.3)
##
###
plot(sig,reject_sx_sd[[1]],xlim=c(0.5,3),typ='o',pch=20,xlab='',ylab='',main=expression('Power of hypothesis test - '~omega^2~'= 0.25'),xaxs="i",xaxt='n',cex=1.5,lwd=2,cex.main=4.5,cex.axis=3.5,ylim=c(0,1))
lines(sig,reject_sy_sd[[1]],col=2,pch=20,lwd=2,cex=1.5,typ='o')
lines(sig,reject_sc_sd[[1]],col=3,pch=20,lwd=2,cex=1.5,typ='o')
abline(v=sig0,lty=2)
mtext(side=1,text=expression('Predictive standard deviation '~sigma),line=7,cex=3.3)
mtext(side=2,text='Power',line=6,cex=3.3)
axis(side=1,at=c(.5,1.5,2,2.5,3.5),label=c(.5,1.5,expression(sigma[0]),2.5,3.5),tick=FALSE,cex.axis=3.7)
legend(x='bottomleft',c(expression(s[0]~'(.,X)'),expression(s[0]~'(.,Y)'),expression(s[v]~'(.,Y)')),col=c(1:4),pch=20,lwd=1.7,cex=4.3,bty='n')
#
plot(sig,reject_sx_sd[[2]],typ='o',pch=20,xlab='',ylab='',main=expression('Power of hypothesis test - '~omega^2~'= 0.5'),xaxs="i",xaxt='n',cex=1.5,lwd=2,cex.main=4.5,cex.axis=3.5,ylim=c(0,1))
lines(sig,reject_sy_sd[[2]],col=2,pch=20,lwd=2,cex=1.5,typ='o')
lines(sig,reject_sc_sd[[2]],col=3,pch=20,lwd=2,cex=1.5,typ='o')
abline(v=sig0,lty=2)
axis(side=1,at=c(.5,1.5,2,2.5,3.5),label=c(.5,1.5,expression(sigma[0]),2.5,3.5),tick=FALSE,cex.axis=3.7)
mtext(side=1,text=expression('Predictive standard deviation '~sigma),line=7,cex=3.3)
#
plot(sig,reject_sx_sd[[3]],typ='o',pch=20,xlab='',ylab='',main=expression('Power of hypothesis test - '~omega^2~'= 1'),xaxs="i",xaxt='n',cex=1.5,lwd=2,cex.main=4.5,cex.axis=3.5,ylim=c(0,1))
lines(sig,reject_sy_sd[[3]],col=2,pch=20,lwd=2,cex=1.5,typ='o')
lines(sig,reject_sc_sd[[3]],col=3,pch=20,lwd=2,cex=1.5,typ='o')
abline(v=sig0,lty=2)
axis(side=1,at=c(.5,1.5,2,2.5,3.5),label=c(.5,1.5,expression(sigma[0]),2.5,3.5),tick=FALSE,cex.axis=3.7)
mtext(side=1,text=expression('Predictive standard deviation '~sigma),line=7,cex=3.3)
dev.off()



#### ------- Figure 6 and 7

source('meanscore_wasserstein_surfaces.r')


#### ------- Figure 8 - Wind data stats and graphs
# X has the distribution of Yobs fitted from joint model and Y is the raw observation
# Wind data: January 2012 of region C1 in (Bessac et al. 2018)
load('BoxCoxTransfObs_Jan12_C1.RData')
s0 <- 5          # station with median MSE 
y <- xx2[,s0]    # box-cox transformed obs data (Yobs)

## From joint GP model:
## NWP forecast
mu1 <- 2.243813
sig1 <- 0.9425318
## Truth X: Yobs 
mu0 <- mean(y)
sig0 <- 1.230149  #(mu_obs from joint GP)
## Verification data Y with level of noise 
c <- 0.5
ybar <- (mu0*c^2)/(sig0^2+c^2) + (y*sig0^2)/(sig0^2+c^2)

## Log-scores 
# ideal mean score
Es1 <- log(sig1) + (sig0^2 + (mu0-mu1)^2)/(2*sig1^2) + log(2*pi)/2     
print('Ideal score')
print(c(Es1))
#
sy <- NULL  # log-score used-in-practise score
sy[[1]] <- log(sig1) + ((y - mu1)^2)/(2*sig1^2) + log(2*pi)/2
print('Scores used in practice')
print(c(mean(sy[[1]]),sd(sy[[1]])))
#
sc <- NULL  # corrected log-score
sc[[1]] <- log(sig1) + ((sig0^2*c^2)/(sig0^2+c^2) + (ybar - mu1)^2)/(2*sig1^2) + log(2*pi)/2
print('Corrected scores')
print(c(mean(sc[[1]]),sd(sc[[1]])))


## CRPS
c0 <- NULL
cy <- NULL
cv <- NULL
#
# ideal CRPS
c0[[1]] <- mu0 - (mu1 + sig1/sqrt(pi)) + 2 *(sqrt(sig1^2+sig0^2)*exp(-(mu0-mu1)^2/(2*(sig1^2+sig0^2)))/sqrt(2*pi) - (mu0 - mu1)*pnorm(q=(mu0-mu1)/sqrt(sig1^2+sig0^2),mean=0,sd=1,lower.tail=FALSE))
print('Ideal CRPS')
print(round(c(c0[[1]]),2))
#
# CRPS computed in practice
cy[[1]] <- y - (mu1+sig1/sqrt(pi)) + 2*sig1*dnorm(x=((y-mu1)/sig1),mean=0,sd=1) -  2*sig1*((y-mu1)/sig1)*pnorm(q=((y-mu1)/sig1),mean=0,sd=1,lower.tail=FALSE)
print('CRPS used in practice')
print(round(c(mean(cy[[1]]),sd(cy[[1]])),2))
#
# corrected CRPS
ybar <- (mu0*c^2)/(sig0^2+c^2) + (y*sig0^2)/(sig0^2+c^2)
ra <- (sig0^2*c^2)/(sig0^2+c^2)
a <- ybar
b <- sqrt(ra)
cv[[1]] <- ybar - (mu1 +sig1/sqrt(pi)) + 2 *( sqrt(sig1^2+b^2)*exp(-(a-mu1)^2/(2*(sig1^2+b^2)))/sqrt(2*pi) - (a - mu1)*pnorm(q=(ybar-mu1)/sqrt(sig1^2+b^2),mean=0,sd=1,lower.tail=FALSE) )  
print('corrected CRPS')
print(round(c(mean(cv[[1]]),sd(cv[[1]])),2))


## plot Figure 8
d1 <- density(sy[[1]])
d2 <- density(sc[[1]])
d3 <- density(cy[[1]])
d4 <- density(cv[[1]])

op <- par(mar = c(5,5,4,1) + 0.1)
layout(matrix(1:2,1,2))
plot(d1$x,d1$y,type="l",xlab="Score value",ylab="PDFs of log-score",xlim=c(0,9),ylim=range(d1$y,d2$y),cex=1.6,cex.lab=1.6,lwd=1.5,col="green")
polygon(c(1,d1$x,d1$x[512]),c(0,d1$y,0),col=rgb(red=0.,green=.2,blue=0.0,alpha=.1),border=F)
grid(lwd=2)
lines(d2$x,d2$y,col="red",lwd=1.5,lty=1)
polygon(c(1,d2$x,d2$x[512]),c(0,d2$y,0),col=rgb(red=0.4,green=.0,blue=0.0,alpha=.2),border=F)
abline(v=Es1,col="orange", lwd=2)
legend(x='topright',c(expression(s[0]~'(f,Y)'),expression(s[v]~'(f,Y)')),lty=c(1,1),col=c(3,2),cex=1.6,lwd=1.5,bty='n')
axis(side=3,at=Es1,label=expression('E('~s[0]~'(f,X))'),tick=FALSE,cex.axis=1.3,line=-1)
mtext(side=3,'Log-score distribution',line=1.8,cex=1.7)
##
##
plot(d3$x,d3$y,type="l",xlab="Score value",ylab="PDFs of CRPS",xlim=c(-0.2,3.5),ylim=range(d3$y,d4$y),cex=1.6,cex.lab=1.6,lwd=1.5,col="green")
polygon(c(1,d3$x,d3$x[512]),c(0,d3$y,0),col=rgb(red=0.,green=.2,blue=0.0,alpha=.1),border=F)
grid(lwd=2)
lines(d4$x,d4$y,col="red",lwd=1.5,lty=1)
polygon(c(1,d4$x,d4$x[512]),c(0,d4$y,0),col=rgb(red=0.4,green=.0,blue=0.0,alpha=.2),border=F)
abline(v=mean(c0[[1]]),col="orange", lwd=2)
legend(x='topright',c(expression(c[0]~'(f,Y)'),expression(c[v]~'(f,Y)')),lty=c(1,1),col=c(3,2),cex=1.6,lwd=1.5,bty='n')
axis(side=3,at=mean(c0[[1]]),label=expression('E('~c[0]~'(f,X))'),tick=FALSE,cex.axis=1.3,line=-1)
mtext(side=3,'CRPS distribution',line=1.8,cex=1.7)



