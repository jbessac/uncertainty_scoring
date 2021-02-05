rm(list=ls())

library('fields')
library('pals')
library('transport')

#### ------- graph 3 - Log-score distributions for additive model 
## ------------------ computing scores
n <- 10000
# true hidden state 
mu0 <- 0
sig0 <- 2
x0 <- rnorm(n=10000,mean=0,sd=1)
x <- mu0 + sig0*x0

# corrupted observations
# observational noise level 
omega0 <- sqrt(0.5)
omega1 <- sqrt(1)
y0 <- x + omega0*rnorm(n=n,mean=0,sd=1)
y1 <- x + omega1*rnorm(n=n,mean=0,sd=1)
den0 <- sig0^2 + omega0^2
den1 <- sig0^2 + omega1^2
ybar0 <- (omega0^2*mu0+y0*sig0^2)/den0
ybar1 <- (omega1^2*mu0+y1*sig0^2)/den1

## --- forecast with varying mean and varying variance 
mup <- seq(mu0-3,mu0+3,by=.05)
sp <- seq(.5,3,by=.05)

ra1 <- (sig0^2*omega1^2)/(sig0^2+omega1^2)
b1 <- sqrt(ra1)

sx <- array(0,c(length(mup),length(sp),n))
sy1 <- array(0,c(length(mup),length(sp),n))
sc1 <- array(0,c(length(mup),length(sp),n))
for (i in 1:length(mup)){
  for (j in 1:length(sp)){
    # ideal log-score
    sx[i,j,] <- log(2*pi)/2 + log(sp[j]) + (x-mup[i])^2/(2*sp[j]^2) 
    # log-score used in practice
    sy1[i,j,] <- log(2*pi)/2 + log(sp[j]) + (y1-mup[i])^2/(2*sp[j]^2)
    # corrected log-score
    sc1[i,j,] <- log(2*pi)/2 + log(sp[j]) + ((omega1^2*sig0^2)/den1 + (ybar1-mup[i])^2)/(2*sp[j]^2)
  }
}


## -------- log-score stats
mu_sx <- apply(X=sx,FUN=mean,MARGIN=1:2)
mu_sx0 <- mu_sx - min(mu_sx)
mu_sy <- apply(X=sy1,FUN=mean,MARGIN=1:2)
mu_sc <- apply(X=sc1,FUN=mean,MARGIN=1:2)

was_sx <- apply(X=sx,MARGIN=1:2,FUN=function(x){wasserstein1d(x,sx[61,31,])})
was_sy <- apply(X=sy1,MARGIN=1:2,FUN=function(x){wasserstein1d(x,sy1[61,31,])})
was_sc <- apply(X=sc1,MARGIN=1:2,FUN=function(x){wasserstein1d(x,sc1[61,31,])})

## relative differences 
mu_zy <- (mu_sx - mu_sy)/mu_sx

zlim0 <- range(mu_zy,na.rm=TRUE)
rbPal <- colorRampPalette(c('blue', 'gray96', 'red')) 
nc <- 100 
max_absolute_value <- max(abs(zlim0)) 
brks <- seq(-max_absolute_value,max_absolute_value,length.out=nc+1)
Col <- rbPal(nc)
n_in_class <- hist(mu_zy, breaks=brks, plot=F)$counts>0
col_to_include <- min(which(n_in_class==T)):max(which(n_in_class==T))
brks_to_include <- min(which(n_in_class==T)):(max(which(n_in_class==T))+1)

zlim1 <- range(mu_sx0,was_sx,was_sy,was_sc,na.rm=TRUE)
z1 <- c(mu_sx0,was_sx,was_sy,was_sc)
nc1 <- 24
Col1 <- parula(nc1)
brks1 <- c(seq(0,2,by=.2),3:10,seq(13,max(z1),by=3),max(z1))

zlevels1 <- c(10,7,4,3,2.5,2,1.5,1,.5,.25,.1)
zlevels0 <- c(-.25,-.2,-.15,-.1,-0.05,0,.05,.1,.2,.25)
xlab <- expression('Predictive mean'~mu)
ylab <- expression('Predictive s.d.'~sigma)
scl <- 2.3
#
setEPS()
postscript("figure7.eps",height = 5, width = 13)
par(mar=c(25,20,10,25),mai=c(.9,1,.8,1.9))
layout(matrix(1:2,1,2))
image(mup,sp,(mu_sx0),main=expression(E(s[0](f,X))-E(s[0](f[0],X))),xlab='',ylab='',col=Col1,breaks=brks1,cex=scl,cex.axis=scl,cex.lab=scl,cex.main=1.8,xaxt='n',yaxt='n',zlim=zlim1)
image.plot(mup,sp,(mu_sx0),col=Col1,breaks=brks1,legend.only=TRUE,legend.shrink=.6,axis.args=list(cex.axis=1.8))
contour(mup,sp,(mu_sx0),col='white',add=TRUE,levels=zlevels1,labcex=1.4,method='edge')
axis(side=1,at=-3:3,label=c(-3,-2,-1,expression(mu[0]),1, 2,3),tick=FALSE,cex.axis=scl)
axis(side=2,at=seq(0.5,3,by=.5),label=c(0.5,1.0,1.5,expression(sigma[0]),2.5,3.0),tick=FALSE,cex.axis=scl)
mtext(ylab,side=2,line=3,cex=1.8,cex.lab=1.8)
mtext(xlab,side=1,line=3,cex=1.8,cex.lab=1.8)
abline(h=sig0)
abline(v=mu0)
#
image(mup,sp,(mu_zy),main=expression(frac(E(s[0](f,X))-E(s[0](f,Y)),E(s[0](f,X)))),xlab='',ylab='',col=Col[col_to_include],breaks=brks[brks_to_include],cex=scl,cex.axis=scl,cex.lab=scl,cex.main=1.8,zlim=zlim0,xaxt='n',yaxt='n')
contour(mup,sp,(mu_zy),col='black',add=TRUE,levels=zlevels0,labcex=1.4,method='edge')
image.plot(mup,sp,(mu_zy),col=Col[col_to_include],breaks=brks[brks_to_include],legend.shrink=.8,legend.line=-1,legend.only=TRUE,axis.args=list(cex.axis=1.8,line=0),zlim=zlim0)
axis(side=1,at=-3:3,label=c(-3,-2,-1,expression(mu[0]),1, 2,3),tick=FALSE,cex.axis=scl)
axis(side=2,at=seq(0.5,3,by=.5),label=c(0.5,1.0,1.5,expression(sigma[0]),2.5,3.0),tick=FALSE,cex.axis=scl)
mtext(xlab,side=1,line=3,cex=1.8,cex.lab=1.8)
abline(h=sig0)
abline(v=mu0)
dev.off()
####

####
setEPS()
postscript("figure8.eps",height = 5, width = 15)
par(mar=c(25,20,4,35),mai=c(.7,.7,.75,1.1))
layout(matrix(1:3,1,3))
image(mup,sp,(was_sx),main=expression(W[1](s[0](f[0],X),s[0](f,X))),xlab='',ylab='',col=Col1,breaks=brks1,cex=scl,cex.axis=scl,cex.lab=scl,cex.main=scl,xaxt='n',yaxt='n',zlim=zlim1)
contour(mup,sp,(was_sx),col='white',add=TRUE,levels=zlevels1,labcex=1.4,method='edge')
axis(side=1,at=-3:3,label=c(-3,-2,-1,expression(mu[0]),1, 2,3),tick=FALSE,cex.axis=scl)
axis(side=2,at=seq(0.5,3,by=.5),label=c(0.5,1.0,1.5,expression(sigma[0]),2.5,3.0),tick=FALSE,cex.axis=scl)
mtext(xlab,side=1,line=4,cex=1.8,cex.lab=1.8)
mtext(ylab,side=2,line=3,cex=1.8,cex.lab=1.8)
abline(h=sig0)
abline(v=mu0)
#
image(mup,sp,(was_sy),main=expression(W[1](s[0](f[0],Y),s[0](f,Y))),xlab='',ylab='',col=Col1,breaks=brks1,cex=scl,cex.axis=scl,cex.lab=scl,cex.main=scl,zlim=zlim1,xaxt='n',yaxt='n')
contour(mup,sp,(was_sy),col='white',add=TRUE,levels=zlevels1,labcex=1.4,method='edge')
axis(side=1,at=-3:3,label=c(-3,-2,-1,expression(mu[0]),1, 2,3),tick=FALSE,cex.axis=scl)
axis(side=2,at=seq(0.5,3,by=.5),label=c(0.5,1.0,1.5,expression(sigma[0]),2.5,3.0),tick=FALSE,cex.axis=scl)
mtext(xlab,side=1,line=4,cex=1.8,cex.lab=1.8)
abline(h=sig0)
abline(v=mu0)
#
image(mup,sp,(was_sc),main=expression(W[1](s[v](f[0],Y),s[v](f,Y))),xlab='',ylab='',col=Col1,breaks=brks1,cex=scl,cex.axis=scl,cex.lab=scl,cex.main=scl,zlim=zlim1,xaxt='n',yaxt='n')
image.plot(mup,sp,(was_sc),col=Col1,breaks=brks1,legend.only=TRUE,axis.args=list(cex.axis=scl),zlim=zlim1)
contour(mup,sp,(was_sc),col='white',add=TRUE,levels=zlevels1,labcex=1.4,method='edge')
axis(side=1,at=-3:3,label=c(-3,-2,-1,expression(mu[0]),1, 2,3),tick=FALSE,cex.axis=scl)
axis(side=2,at=seq(0.5,3,by=.5),label=c(0.5,1.0,1.5,expression(sigma[0]),2.5,3.0),tick=FALSE,cex.axis=scl)
mtext(xlab,side=1,line=4,cex=1.8,cex.lab=1.8)
abline(h=sig0)
abline(v=mu0)
dev.off()

