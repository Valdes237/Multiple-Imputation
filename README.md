# Multiple-Imputation



#####################################################################
### Comparison: BD and PI
### Quantities of Interes: E(X3) and the betas
### Diagnostics: Bias and Coverage (a=0.05)
### Iterations: 1000
### n: 500 
### missing: 50% MCAR
###
###############################################################################

### The 'MASS' package has a function for drawing from the mvN
library(MASS)

set.seed(1000)

### true values
trueVal <- c(12.8, 5, 0.6, 0.5)

nIter <- 1000
n <- 500
ImpPI <- matrix(nrow=nIter, ncol=12)

### how long does it take?
zeit1 <- Sys.time()

for (i in 1:nIter) {	
  X1 <- rnorm(n,8,3)
  X2 <- 10 - 0.5*X1 + rnorm(n,0,3)
  X3 <- 5 + 0.6*X1 + 0.5*X2 + rnorm(n,0,sqrt(2))
  data1 <- as.data.frame(cbind(X1,X2,X3))
  misind <- sample(1:n,round(n/2))  # fehlende Werte werden generiert
  obsind <- which(!is.na(data1$X3)) # position der beobachten Werte von X3 in data1
  is.na(data1$X3[misind]) <- TRUE  # missung wird in data1$X3 ersetzt
  PI.dat <- data1
  regmod <- lm(X3~X1+X2, data = PI.dat, subset = obsind)
  beta.hat <- coef(regmod)
  sigma.hat <- summary(regmod)$sigma
  Design <- cbind(1,PI.dat$X1, PI.dat$X2)
  n.obs <- length(obsind)
  sig.sq.tilde <- (n.obs-3) * sigma.hat^2/rchisq(1, n.obs-3)
  ## Difference to 'cpd':
  ## Random draw for sigma^2 from a scaled-inverse Chi^2-distribution
  beta.tilde <- mvrnorm(1, beta.hat,
                        solve(t(Design[obsind, ])%*%
                                Design[obsind, ]) * sig.sq.tilde)  # ist gleich wie
  ## beta.tilde <- mvrnorm(1, beta.hat,
                      #  solve(t(Design)%*%
                      #           Design) * sig.sq.tilde) ##
  ## Random draw for beta from a mv. normal distribution
  X3.tilde <- rnorm(nrow(Design),Design%*%beta.tilde, sqrt(sig.sq.tilde)) # vollst?ndige Beobachtungen
                                                                          # von X3
  ## Using draws for beta and sigma^2 rather than the estimators
  ## from the regression
  PI.dat$X3[misind] <- X3.tilde[misind] # die fehlende Werte werden hier durch 
                                        # neue Werte von X3.tilde ersetz
  ##########################################################################
  PI.m <- meanEst(PI.dat)
  PI.r <- modEst(PI.dat)
  ImpPI[i, ] <- c(PI.m,PI.r)
}

save("ImpBD", "ImpCC", "Impcpd", "ImpPI",
     file = file.path(path,"Data/Simulation2.RData"))

zeit2 <- Sys.time()

cat("Overall run-time:",difftime(zeit2,zeit1,units="secs"),"seconds!\n")

### BD mean Coverage
BDcoverE3 <- factor(ImpBD[ ,2], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  Sys.sleep(0.001)
  plot(BDcoverE3[1:i], main = "BD-Coverage E(X3)",col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=4)
}
dev.off()

### Comparison: Coverage level of 'cpd'
cpdCoverX3 <- colMeans(Impcpd)[2]

### PI mean Coverage
PI.coverE3 <- factor(ImpPI[ ,2], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  Sys.sleep(0.001)
  plot(PI.coverE3[1:i], main = "PI-Coverage E(X3)",col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=4)
  abline(h=round(cpdCoverX3*nIter),lty=3,col="orange",lwd=3)
}
dev.off()

### Still too low -- what about beta1:

### BD beta1 Coverage
BDcoverB1 <- factor(ImpBD[ ,8], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  Sys.sleep(0.001)
  plot(BDcoverB1[1:i], main = expression(paste("BD-Coverage ",beta[1])),
       col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=4)
}
dev.off()

cpdCoverB1 <- colMeans(Impcpd)[8]

### PI beta1 Coverage
PI.coverB1 <- factor(ImpPI[ ,8], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  Sys.sleep(0.001)
  plot(PI.coverB1[1:i], main = expression(paste("PI-Coverage ",beta[1])),
       col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=4)
  abline(h=round(cpdCoverB1*nIter),lty=3,col="orange",lwd=3)
}
dev.off()

### E(X3), Coverage and Bias
Overview2 <- rbind(colMeans(ImpBD), colMeans(ImpCC), colMeans(Impcpd),
                   colMeans(ImpPI) )
colnames(Overview2) <- c("E(x3)","C.E(X3)","B.E(X3)", "alpha", "C.alpha",
                         "B.alpha","beta1", "C.beta1", "B.beta1",
                         "beta2", "C.beta2", "B.beta2")
rownames(Overview2) <- c("before deletion", "complete cases",
                         "cond. pred.-draws", "Post.-Imput.-draws")

round(Overview2, 4)

### coverage-wise we did not make any progress, but we did on a conceptual
### level! We just have to think it through completely: One draw does not 
### tell us much about the distribution -- we have to make several (M) draws.
### And then we have to find a smart way to combine the results:
### W: We average the variances from the M versions
### B: We derive the variance of the M estimators
### T: We combine W and B
###
### Yes but...if we increase the variance again...why can we not just use
### the complete cases???!!!

### Good point! We will compare the width of the CIs now too!


load(file = file.path(path,"Data/Simulation2.RData"))
library(MASS)


###################################################################

MI.inference <- function(thetahat, varhat.thetahat, alpha = 0.05)
{
  M <- length(thetahat)
  if (length(varhat.thetahat) != M) {
    stop ("Different length for 'thetahat' and 'varhat.thetahat'!\n")}
  lambda <- 1 - (alpha/2)
  MIestimate <- mean(thetahat)
  B <- var(thetahat)  
  W <- mean(varhat.thetahat)
  total <- W + (1 + 1/M) * B
  DF <- (M - 1) * (1 + W/((1 + 1/M) * B))^2
  CI.low <- MIestimate - qt(lambda, DF) * sqrt(total)
  CI.up <- MIestimate + qt(lambda, DF) * sqrt(total)
  list(MI.Est = MIestimate, MI.Var = total, CI.low = CI.low,
       CI.up = CI.up, BVar = B, WVar = W)
}


modEst2 <- function(data, alpha=0.05, Wert=trueVal[-1],...) {
  model <- lm(X3~X1+X2, data=data, na.action = "na.omit")
  beta.hat <- model$coefficients
  se.beta.hat <- summary(model)$coefficients[ ,2]
  ## CI lower bound
  CI.low <- beta.hat - qnorm(1-alpha/2)*se.beta.hat
  ## CI upper bound
  CI.hi <- beta.hat + qnorm(1-alpha/2)*se.beta.hat
  contained <- as.numeric(CI.low <= Wert & CI.hi >= Wert)
  bias <- beta.hat-Wert
  CI.range <- CI.hi-CI.low
  Imp <- vector(length = 12)
  Imp[seq(1,by = 4,length.out = 3)] <- beta.hat
  Imp[seq(2,by = 4,length.out = 3)] <- contained
  Imp[seq(3,by = 4,length.out = 3)] <- bias
  Imp[seq(4,by = 4,length.out = 3)] <- CI.range
  return(Imp)
}

meanEst2 <- function(data, alpha=0.05, Wert=trueVal[1],...) {
  mean.X3 <- mean(data$X3, na.rm=TRUE)
  ## CI lower bound
  CI.low <- mean.X3 - qnorm(1-alpha/2)*sqrt(var(data$X3,
                                                na.rm=TRUE)/sum(!is.na(data$X3)))
  ## CI upper bound
  CI.hi <- mean.X3 + qnorm(1-alpha/2)*sqrt(var(data$X3,
                                               na.rm=TRUE)/sum(!is.na(data$X3)))
  contained <- as.numeric(CI.low <= Wert & CI.hi >= Wert)
  bias <- mean.X3-Wert
  CI.range <- CI.hi-CI.low
  Imp <- c(mean.X3, contained, bias, CI.range)
  return(Imp)
}


MIEst <- function(MI.mean=MI.mean, MI.a=MI.a, MI.b1=MI.b1, MI.b2=MI.b2,
                  Wert=trueVal,...){
  enth <- vector(length = 4)
  enth[1] <- as.numeric(MI.mean$CI.low <= Wert[1] & MI.mean$CI.up >= Wert[1])
  enth[2] <- as.numeric(MI.a$CI.low <= Wert[2] & MI.a$CI.up >= Wert[2])
  enth[3] <- as.numeric(MI.b1$CI.low <= Wert[3] & MI.b1$CI.up >= Wert[3])
  enth[4] <- as.numeric(MI.b2$CI.low <= Wert[4] & MI.b2$CI.up >= Wert[4])
  MI.Est <- vector(length = 4)
  MI.Est[1] <- MI.mean$MI.Est
  MI.Est[2] <- MI.a$MI.Est
  MI.Est[3] <- MI.b1$MI.Est
  MI.Est[4] <- MI.b2$MI.Est
  bias <- vector(length = 4)
  bias[1] <- MI.mean$MI.Est-Wert[1]
  bias[2] <- MI.a$MI.Est-Wert[2]
  bias[3] <- MI.b1$MI.Est-Wert[3]
  bias[4] <- MI.b2$MI.Est-Wert[4]
  CI.range <- vector(length = 4)
  CI.range[1] <- MI.mean$CI.up-MI.mean$CI.low
  CI.range[2] <- MI.a$CI.up-MI.a$CI.low
  CI.range[3] <- MI.b1$CI.up-MI.b1$CI.low
  CI.range[4] <- MI.b2$CI.up-MI.b2$CI.low
  Imp <- vector(length = 16)
  Imp[seq(1,by = 4,length.out = 4)]  <- MI.Est
  Imp[seq(2,by = 4,length.out = 4)]  <- enth
  Imp[seq(3,by = 4,length.out = 4)]  <- bias
  Imp[seq(4,by = 4,length.out = 4)]  <- CI.range
  return(Imp)
}


set.seed(2000)

### true values
trueVal <- c(12.8, 5, 0.6, 0.5)

nIter <- 1000
n <- 500
M <- 10 # 10 Imputations
ImpBD2 <- ImpCC2 <- ImpMI <- matrix(nrow=nIter, ncol=16)

### how long does it take?
zeit1 <- Sys.time()

for (i in 1:nIter) {	
  X1 <- rnorm(n,8,3)
  X2 <- 10 - 0.5*X1 + rnorm(n,0,3)
  X3 <- 5 + 0.6*X1 + 0.5*X2 + rnorm(n,0,sqrt(2))
  data1 <- as.data.frame(cbind(X1,X2,X3))
  BD.dat <- data1 # BD (before deletion)
  n <- nrow(BD.dat)
  misind <- sample(1:n,round(n/2))
  obsind <- which(!is.na(data1$X3))
  is.na(data1$X3[misind]) <- TRUE
  CC.dat <- data1 # CC (complete cases)
  MI.hilfdat <- matrix(nrow=M, ncol=8)
  for (m in 1:M) {
    MI.dat <- data1
    regmod <- lm(X3~X1+X2, data = MI.dat, subset = obsind)
    beta.hat <- coef(regmod)
    sigma.hat <- summary(regmod)$sigma
    Design <- cbind(1,MI.dat$X1, MI.dat$X2)
    n.obs <- length(obsind)
    sig.sq.tilde <- (n.obs-3) * sigma.hat^2/rchisq(1, n.obs-3)
    ## Differences to 'cpd':
    ## Random draws for sigma^2 from the scaled inverse Chi^2-distribution
    beta.tilde <- mvrnorm(1, beta.hat,
                          solve(t(Design[obsind, ])%*%
                                  Design[obsind, ]) * sig.sq.tilde)
    ## Random draw for beta from a mv. normal distribution
    X3.tilde <- rnorm(nrow(Design),Design%*%beta.tilde, sqrt(sig.sq.tilde))
    MI.dat$X3[misind] <- X3.tilde[misind] 
    MI.hilfdat[m, 1] <- mean(MI.dat$X3)
    MI.hilfdat[m, 2] <- var(MI.dat$X3)/n
    model <- lm(X3~X1+X2, data=MI.dat, na.action = "na.omit")
    beta.hat <- model$coefficients
    se.beta.hat <- summary(model)$coefficients[ ,2]
    MI.hilfdat[m, seq(3,by=2, length.out=3)] <- beta.hat
    MI.hilfdat[m, seq(4,by=2, length.out=3)] <- se.beta.hat^2
  }
  ##########################################################################
  ## MI Analysis / combining rules
  MI.mean <- MI.inference(MI.hilfdat[ ,1], MI.hilfdat[ ,2])
  MI.a <- MI.inference(MI.hilfdat[ ,3], MI.hilfdat[ ,4])
  MI.b1 <- MI.inference(MI.hilfdat[ ,5], MI.hilfdat[ ,6])
  MI.b2 <- MI.inference(MI.hilfdat[ ,7], MI.hilfdat[ ,8])
  ## MI CIs finished -- now we have to fill the result matrix
  ImpMI[i, ] <- MIEst(MI.mean,MI.a,MI.b1,MI.b2)
  BD.m <- meanEst2(BD.dat)
  BD.r <- modEst2(BD.dat)
  ImpBD2[i, ] <- c(BD.m,BD.r)
  CC.m <- meanEst2(CC.dat)
  CC.r <- modEst2(CC.dat)
  ImpCC2[i, ] <- c(CC.m,CC.r)
}

save("ImpBD2", "ImpCC2", "ImpMI", "MI.inference", "MIEst", "modEst2",
     "meanEst2", file = file.path(path,"Data/Simulation3.RData"))

zeit2 <- Sys.time()

cat("Overall run-time:",difftime(zeit2,zeit1,units="secs"),"seconds!\n")


### again for comparison: the coverage level of 'cpd'
cpdCoverX3 <- colMeans(Impcpd)[2]

### PI mean Coverage
MI.coverE3 <- factor(ImpMI[ ,2], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  plot(MI.coverE3[1:i], main = "MI-Coverage E(X3)",col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=4)
  abline(h=round(cpdCoverX3*nIter),lty=3,col="orange",lwd=3)
}
dev.off()

### now we are checking all paramters
### MI alpha Coverage
cpdCover.a <- colMeans(Impcpd)[5]
MI.cover.a <- factor(ImpMI[ ,6], labels = c("not within CI", "within CI"))

windows()
for (i in 1:nIter){
  plot(MI.cover.a[1:i], main = expression(paste("MI-Coverage ",alpha)),
       col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=4)
  abline(h=round(cpdCover.a*nIter),lty=3,col="orange",lwd=3)
}
dev.off()

### MI beta1 Coverage
cpdCover.b1 <- colMeans(Impcpd)[8]
MI.cover.b1 <- factor(ImpMI[ ,10], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  plot(MI.cover.b1[1:i], main = expression(paste("MI-Coverage ",beta[1])),
       col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=4)
  abline(h=round(cpdCover.b1*nIter),lty=3,col="orange",lwd=3)
}
dev.off()

### MI beta2 Coverage
cpdCover.b2 <- colMeans(Impcpd)[11]
MI.cover.b2 <- factor(ImpMI[ ,14], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  plot(MI.cover.b2[1:i], main = expression(paste("MI-Coverage ",beta[2])),
       col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=4)
  abline(h=round(cpdCover.b2*nIter),lty=3,col="orange",lwd=3)
}
dev.off()

### E(X3), Coverage, Bias and range of CIs
Overview3 <- rbind(colMeans(ImpBD2), colMeans(ImpCC2), colMeans(ImpMI) )
colnames(Overview3) <- c("E(x3)","C.E(X3)","B.E(X3)", "S.E(X3)", "alpha",
                         "C.alpha", "B.alpha", "S.alpha", "beta1", "C.beta1",
                         "B.beta1", "S.beta1", "beta2", "C.beta2",
                         "B.beta2","S.beta2")
rownames(Overview3) <- c("before deletion", "complete cases", "MI")

round(Overview3, 4)

rangC <- Overview3[,grep("S.",colnames(Overview3))]
matplot(log(rangC),type="b")

############################################################################
### modified regression model: X1 = b0 + b1X2 + b2X3


modEst3 <- function(data, alpha=0.05, Wert=trueVal[-1],...) {
  model <- lm(X1~X2+X3, data=data, na.action = "na.omit")
  beta.hat <- model$coefficients
  se.beta.hat <- summary(model)$coefficients[ ,2]
  ## CI lower bound
  CI.low <- beta.hat - qnorm(1-alpha/2)*se.beta.hat
  ## CI upper bound
  CI.hi <- beta.hat + qnorm(1-alpha/2)*se.beta.hat
  contained <- as.numeric(CI.low <= Wert & CI.hi >= Wert)
  bias <- beta.hat-Wert
  CI.range <- CI.hi-CI.low
  Imp <- vector(length = 12)
  Imp[seq(1,by = 4,length.out = 3)] <- beta.hat
  Imp[seq(2,by = 4,length.out = 3)] <- contained
  Imp[seq(3,by = 4,length.out = 3)] <- bias
  Imp[seq(4,by = 4,length.out = 3)] <- CI.range
  return(Imp)
}


### true values for new regression model are a little bit tedious...

### Var(X2)=0.5^2*Var(X1)+3^2
### Var(X3)=0.6^2*Var(X1)+0.5^2*Var(X2)+2*0.6*0.5*Cov(X1,X2)
### -0.5=Cov(X1,X2)/Var(X1)
covX1X2 <- -4.5
VX3 <- 0.6^2*9+0.5^2*11.25+2*0.6*0.5*covX1X2+2
varX <- c(9, 11.25, VX3)

### b1 = Cov(X1,X3)/Var(X1)-b2*Cov(x1,X2)/Var(X1)
### --> Cov(x1,X3)=b1*Var(X1)+b2*Cov(x1,X2)
covX1X3 <- 0.6*varX[1]+0.5*covX1X2

### b2 = Cov(X2,X3)/Var(X2)-b1*Cov(x1,X2)/Var(X2)
### --> Cov(x2,X3)=b2*Var(X2)+b1*Cov(x1,X2)
covX2X3 <- 0.5*varX[2]+0.6*covX1X2

### Covariance Matrix
Sigma <- matrix(c(varX[1],covX1X2,covX1X3,
                  covX1X2,varX[2],covX2X3,
                  covX1X3,covX2X3,varX[3]),nrow=3)

### now: X1 dependent variable
### beta1 and beta2
beta12 <- as.numeric(c(Sigma[1,c(2,3)])%*%solve(Sigma[-1,-1])) 
beta0 <- 8-beta12[1]*6-beta12[2]*12.8

### true values
trueVal <- c(12.8, beta0, beta12)

#################################################################################

nIter <- 1000
n <- 500
M <- 10 # 10 Imputations
ImpBD3 <- ImpCC3 <- ImpMI2 <- matrix(nrow=nIter, ncol=16)

### how long does it take?
zeit1 <- Sys.time()

for (i in 1:nIter) {	
  X1 <- rnorm(n,8,3)
  X2 <- 10 - 0.5*X1 + rnorm(n,0,3)
  X3 <- 5 + 0.6*X1 + 0.5*X2 + rnorm(n,0,sqrt(2))
  data1 <- as.data.frame(cbind(X1,X2,X3))
  BD.dat <- data1 # BD (before deletion)
  n <- nrow(BD.dat)
  misind <- sample(1:n,round(n/2))
  obsind <- which(!is.na(data1$X3))
  is.na(data1$X3[misind]) <- TRUE
  CC.dat <- data1 # CC (complete cases)
  MI.hilfdat <- matrix(nrow=M, ncol=8)
  for (m in 1:M) {
    MI.dat <- data1
    regmod <- lm(X3~X1+X2, data = MI.dat, subset = obsind)
    beta.hat <- coef(regmod)
    sigma.hat <- summary(regmod)$sigma
    Design <- cbind(1,MI.dat$X1, MI.dat$X2)
    n.obs <- length(obsind)
    sig.sq.tilde <- (n.obs-3) * sigma.hat^2/rchisq(1, n.obs-3)
    beta.tilde <- mvrnorm(1, beta.hat,
                          solve(t(Design[obsind, ])%*%
                                  Design[obsind, ]) * sig.sq.tilde)
    X3.tilde <- rnorm(nrow(Design),Design%*%beta.tilde, sqrt(sig.sq.tilde))
    MI.dat$X3[misind] <- X3.tilde[misind] 
    MI.hilfdat[m, 1] <- mean(MI.dat$X3)
    MI.hilfdat[m, 2] <- var(MI.dat$X3)/n
    model <- lm(X1~X2+X3, data=MI.dat, na.action = "na.omit")
    beta.hat <- model$coefficients
    se.beta.hat <- summary(model)$coefficients[ ,2]
    MI.hilfdat[m, seq(3,by=2, length.out=3)] <- beta.hat
    MI.hilfdat[m, seq(4,by=2, length.out=3)] <- se.beta.hat^2
  }
  ##########################################################################
  ## MI-Analysis / combining rules
  MI.mean <- MI.inference(MI.hilfdat[ ,1], MI.hilfdat[ ,2])
  MI.a <- MI.inference(MI.hilfdat[ ,3], MI.hilfdat[ ,4])
  MI.b1 <- MI.inference(MI.hilfdat[ ,5], MI.hilfdat[ ,6])
  MI.b2 <- MI.inference(MI.hilfdat[ ,7], MI.hilfdat[ ,8])
  
  ImpMI2[i, ] <- MIEst(MI.mean,MI.a,MI.b1,MI.b2)
  BD.m <- meanEst2(BD.dat)
  BD.r <- modEst3(BD.dat)
  ImpBD3[i, ] <- c(BD.m,BD.r)
  CC.m <- meanEst2(CC.dat)
  CC.r <- modEst3(CC.dat)
  ImpCC3[i, ] <- c(CC.m,CC.r)
}

save("ImpBD3", "ImpCC3", "ImpMI2", "MI.inference", "MIEst", "modEst2",
     "meanEst2", file = file.path(path,"Data/Simulation4.RData"))

zeit2 <- Sys.time()

cat("Overall run-time:",difftime(zeit2,zeit1,units="secs"),
    "seconds!\n")


### MI mean Coverage
MI.coverE3 <- factor(ImpMI2[ ,2], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  plot(MI.coverE3[1:i], main = "MI-Coverage E(X3)",col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=4)
}
dev.off()

### MI alpha Coverage
MI.cover.a <- factor(ImpMI2[ ,6], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  plot(MI.cover.a[1:i], main = expression(paste("MI-Coverage ",alpha)),
       col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=4)
}
dev.off()

### MI beta1 Coverage
MI.cover.b1 <- factor(ImpMI2[ ,10], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  plot(MI.cover.b1[1:i], main = expression(paste("MI-Coverage ",beta[1])),
       col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=4)
}
dev.off()

### MI beta2 Coverage
MI.cover.b2 <- factor(ImpMI2[ ,14], labels = c("not within CI", "within CI"))
windows()
for (i in 1:nIter){
  plot(MI.cover.b2[1:i], main = expression(paste("MI-Coverage ",beta[2])),
       col=c("red","green"),
       ylim = c(0,1.01*nIter))
  abline(h=round(0.95*nIter),lty=2,col="blue",lwd=4)
}
dev.off()

### E(X3), Coverage, Bias and range of CIs
Overview4 <- rbind(colMeans(ImpBD3), colMeans(ImpCC3), colMeans(ImpMI2) )
colnames(Overview4) <- c("E(x3)","C.E(X3)","B.E(X3)", "S.E(X3)", "alpha",
                         "C.alpha", "B.alpha", "S.alpha", "beta1", "C.beta1",
                         "B.beta1", "S.beta1", "beta2", "C.beta2",
                         "B.beta2","S.beta2")
rownames(Overview4) <- c("before deletion", "complete cases", "MI")

round(Overview4, 4)

rangC <- Overview4[,grep("S.",colnames(Overview4))]
matplot(log(rangC),type="b")
