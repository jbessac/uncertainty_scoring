
rm(list=ls())

### --- true parameters
mu0 <- 0
sig0 <- 2

### ------------------------------------------------- statistical test and associated power

mu <- seq(mu0-3,mu0+3,by=.05)

Nsim <- 10000
reject_sx_mu <- NULL ; accept_sx <- NULL
reject_sy_mu <- NULL ; accept_sy <- NULL
reject_sc_mu <- NULL ; accept_sc <- NULL

c <- sqrt(.25)
rej_sx <- rep(0,length(mu)) ; acc_sx <- rep(0,length(mu))
rej_sy <- rep(0,length(mu)) ; acc_sy <- rep(0,length(mu))
rej_sc <- rep(0,length(mu)) ; acc_sc <- rep(0,length(mu))
for (s in 1:Nsim){
    x0 <- rnorm(n=1000,mean=0,sd=1)
    for (i in 1:length(mu)){
        x <- mu0 + sig0*x0
        y <- mu0 + sqrt(sig0^2+c^2)*x0
        ybar <- (mu0*c^2)/(sig0^2+c^2) + (y*sig0^2)/(sig0^2+c^2)
        Es0 <- log(sig0) + 1/2 + log(2*pi)/2
        sx <- (log(sig0) + ((x - mu[i])^2)/(2*sig0^2) + log(2*pi)/2)
        sy <- (log(sig0) + ((y - mu[i])^2)/(2*sig0^2) + log(2*pi)/2)
        sc <- (log(sig0) + ((sig0^2*c^2)/(sig0^2+c^2) + (ybar - mu[i])^2)/(2*sig0^2) + log(2*pi)/2)
        ICix <- Es0 - 1.96*sd(sx)/sqrt(length(x)) ; ICsx <- Es0 + 1.96*sd(sx)/sqrt(length(x))
        ICiy <- Es0 - 1.96*sd(sy)/sqrt(length(x)) ; ICsy <- Es0 + 1.96*sd(sy)/sqrt(length(x))
        ICic <- Es0 - 1.96*sd(sc)/sqrt(length(x)) ; ICsc <- Es0 + 1.96*sd(sc)/sqrt(length(x))
        if (( mean(sx) > ICix )&( mean(sx) < ICsx )){ acc_sx[i] <- acc_sx[i] +1 }
        else { rej_sx[i] <- rej_sx[i] +1 }
        if (( mean(sy) > ICiy )&( mean(sy) < ICsy )){ acc_sy[i] <- acc_sy[i] +1 }
        else { rej_sy[i] <- rej_sy[i] +1 }
        if (( mean(sc) > ICic )&( mean(sc) < ICsc )){ acc_sc[i] <- acc_sc[i] +1 }
        else { rej_sc[i] <- rej_sc[i] +1 }  }}
accept_sx[[1]] <- acc_sx/Nsim ; reject_sx_mu[[1]] <- rej_sx/Nsim
accept_sy[[1]] <- acc_sy/Nsim ; reject_sy_mu[[1]] <- rej_sy/Nsim
accept_sc[[1]] <- acc_sc/Nsim ; reject_sc_mu[[1]] <- rej_sc/Nsim

c <- sqrt(0.5)
rej_sx <- rep(0,length(mu)) ; acc_sx <- rep(0,length(mu))
rej_sy <- rep(0,length(mu)) ; acc_sy <- rep(0,length(mu))
rej_sc <- rep(0,length(mu)) ; acc_sc <- rep(0,length(mu))
rej_med_sx <- rep(0,length(mu)) ; acc_med_sx <- rep(0,length(mu))
rej_med_sy <- rep(0,length(mu)) ; acc_med_sy <- rep(0,length(mu))
rej_med_sc <- rep(0,length(mu)) ; acc_med_sc <- rep(0,length(mu))
for (s in 1:Nsim){
    for (i in 1:length(mu)){
        x0 <- rnorm(n=1000,mean=0,sd=1)
        x <- mu0 + sig0*x0
        y <- mu0 + sqrt(sig0^2+c^2)*x0
        ybar <- (mu0*c^2)/(sig0^2+c^2) + (y*sig0^2)/(sig0^2+c^2)
        Es0 <- log(sig0) + 1/2 + log(2*pi)/2
        sx <- log(sig0) + ((x - mu[i])^2)/(2*sig0^2) + log(2*pi)/2
        sy <- log(sig0) + ((y - mu[i])^2)/(2*sig0^2) + log(2*pi)/2
        sc <- log(sig0) + ((sig0^2*c^2)/(sig0^2+c^2) + (ybar - mu[i])^2)/(2*sig0^2) + log(2*pi)/2
        ICix <- Es0 - 1.96*sd(sx)/sqrt(length(x)) ; ICsx <- Es0 + 1.96*sd(sx)/sqrt(length(x))
        ICiy <- Es0 - 1.96*sd(sy)/sqrt(length(x)) ; ICsy <- Es0 + 1.96*sd(sy)/sqrt(length(x))
        ICic <- Es0 - 1.96*sd(sc)/sqrt(length(x)) ; ICsc <- Es0 + 1.96*sd(sc)/sqrt(length(x))
        # power of the mean score 
        if (( mean(sx) > ICix )&( mean(sx) < ICsx )){ acc_sx[i] <- acc_sx[i] +1 }
        else { rej_sx[i] <- rej_sx[i] +1 }
        if (( mean(sy) > ICiy )&( mean(sy) < ICsy )){ acc_sy[i] <- acc_sy[i] +1 }
        else { rej_sy[i] <- rej_sy[i] +1 }
        if (( mean(sc) > ICic )&( mean(sc) < ICsc )){ acc_sc[i] <- acc_sc[i] +1 }
        else { rej_sc[i] <- rej_sc[i] +1 }
        }}
accept_sx[[2]] <- acc_sx/Nsim ; reject_sx_mu[[2]] <- rej_sx/Nsim
accept_sy[[2]] <- acc_sy/Nsim ; reject_sy_mu[[2]] <- rej_sy/Nsim
accept_sc[[2]] <- acc_sc/Nsim ; reject_sc_mu[[2]] <- rej_sc/Nsim
#

c <- sqrt(1)
rej_sx <- rep(0,length(mu)) ; acc_sx <- rep(0,length(mu))
rej_sy <- rep(0,length(mu)) ; acc_sy <- rep(0,length(mu))
rej_sc <- rep(0,length(mu)) ; acc_sc <- rep(0,length(mu))
for (s in 1:Nsim){
    for (i in 1:length(mu)){
        x0 <- rnorm(n=1000,mean=0,sd=1)
        x <- mu0 + sig0*x0
        y <- mu0 + sqrt(sig0^2+c^2)*x0
        ybar <- (mu0*c^2)/(sig0^2+c^2) + (y*sig0^2)/(sig0^2+c^2)
        Es0 <- log(sig0) + 1/2 + log(2*pi)/2
        sx <- log(sig0) + ((x - mu[i])^2)/(2*sig0^2) + log(2*pi)/2
        sy <- log(sig0) + ((y - mu[i])^2)/(2*sig0^2) + log(2*pi)/2
        sf <- log(sig0) + ((y - mu[i])^2 - c^2)/(2*sig0^2) + log(2*pi)/2
        sc <- log(sig0) + ((sig0^2*c^2)/(sig0^2+c^2) + (ybar - mu[i])^2)/(2*sig0^2) + log(2*pi)/2
        ICix <- Es0 - 1.96*sd(sx)/sqrt(length(x)) ; ICsx <- Es0 + 1.96*sd(sx)/sqrt(length(x))
        ICiy <- Es0 - 1.96*sd(sy)/sqrt(length(x)) ; ICsy <- Es0 + 1.96*sd(sy)/sqrt(length(x))
        ICic <- Es0 - 1.96*sd(sc)/sqrt(length(x)) ; ICsc <- Es0 + 1.96*sd(sc)/sqrt(length(x))
        if (( mean(sx) > ICix )&( mean(sx) < ICsx )){ acc_sx[i] <- acc_sx[i] +1 }
        else { rej_sx[i] <- rej_sx[i] +1 }
        if (( mean(sy) > ICiy )&( mean(sy) < ICsy )){ acc_sy[i] <- acc_sy[i] +1 }
        else { rej_sy[i] <- rej_sy[i] +1 }
        if (( mean(sc) > ICic )&( mean(sc) < ICsc )){ acc_sc[i] <- acc_sc[i] +1 }
        else { rej_sc[i] <- rej_sc[i] +1 }  }}
accept_sx[[3]] <- acc_sx/Nsim ; reject_sx_mu[[3]] <- rej_sx/Nsim
accept_sy[[3]] <- acc_sy/Nsim ; reject_sy_mu[[3]] <- rej_sy/Nsim
accept_sc[[3]] <- acc_sc/Nsim ; reject_sc_mu[[3]] <- rej_sc/Nsim
