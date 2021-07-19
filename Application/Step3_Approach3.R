###Application with MGI data cont age cat bmi
############Application
library(statmod)
library(robustbase)
#library(ggplot2)
#library(ggpubr)
library(reshape2)
#library(gridExtra)
#library(grid)
library(R.utils)
library(data.table)
library(sampling)
library(dplyr)
#library(arsenal)
library(betareg)

myremlscore = function (y, X, Z, trace = FALSE, tol = 1e-05, maxit = 40)
{
  conv=1
  n <- length(y)
  p <- dim(X)[2]
  q <- dim(Z)[2]
  const <- n * log(2 * pi)
  fitm <- lm.fit(X, y)
  if (fitm$qr$rank < p)
    stop("X is of not of full column rank")
  Q <- qr.Q(fitm$qr)
  h <- as.vector(Q^2 %*% array(1, c(p, 1)))
  d <- fitm$residuals^2
  wd <- 1 - h
  zd <- log(d/(1 - h)) + 1.27
  fitd <- lm.wfit(Z, zd, wd)
  gam <- ifelse(is.na(fitd$coef), 0, fitd$coef)
  g <- fitd$fitted.values
  phi <- exp(g)
  wm <- 1/phi
  fitm <- lm.wfit(X, y, wm)
  d <- fitm$residuals^2
  dev <- sum(d/phi) + sum(log(phi)) + const + 2 * log(prod(abs(diag(fitm$qr$qr))))
  iter <- 0
  if (trace)
    cat("Iter =", iter, ", Dev =", dev, " Gamma", gam, "\n")
  Q2 <- array(0, c(n, p * (p + 1)/2))
  repeat {
    iter <- iter + 1
    Q <- qr.qy(fitm$qr, diag(1, nrow = n, ncol = p))
    j0 <- 0
    for (k in 0:(p - 1)) {
      Q2[, (j0 + 1):(j0 + p - k)] <- Q[, 1:(p - k)] * Q[,
                                                        (k + 1):p]
      j0 <- j0 + p - k
    }
    if (p > 1)
      Q2[, (p + 1):(p * (p + 1)/2)] <- sqrt(2) * Q2[, (p +
                                                         1):(p * (p + 1)/2)]
    h <- drop(Q2[, 1:p] %*% array(1, c(p, 1)))
    Q2Z <- t(Q2) %*% Z
    ZVZ <- (t(Z) %*% vecmat(1 - 2 * h, Z) + t(Q2Z) %*% Q2Z)/2
    maxinfo <- max(diag(ZVZ))
    if (iter == 1) {
      lambda <- abs(mean(diag(ZVZ)))/q
      I <- diag(q)
    }
    zd <- (d - (1 - h) * phi)/phi
    dl <- crossprod(Z, zd)/2
    gamold <- gam
    devold <- dev
    lev <- 0
    repeat {
      lev <- lev + 1
      R <- chol(ZVZ + lambda * I)
      dgam <- backsolve(R, backsolve(R, dl, transpose = TRUE))
      gam <- gamold + dgam
      phi <- as.vector(exp(Z %*% gam))
      wm <- 1/phi
      fitm <- lm.wfit(X, y, wm)
      d <- fitm$residuals^2
      dev <- sum(d/phi) + sum(log(phi)) + const + 2 * log(prod(abs(diag(fitm$qr$qr))))
      if (dev < devold - 1e-15)
        break
      if (lambda/maxinfo > 1e+15) {
        gam <- gamold
        warning("Too much damping - convergence tolerance not achievable")
        break
      }
      lambda <- 2 * lambda
      if (trace)
        cat("Damping increased to", lambda, "\n")
    }
    if (trace)
      cat("Iter =", iter, ", Dev =", dev, " Gamma", gam,
          "\n")
    if (lambda/maxinfo > 1e+15)
      break
    if (lev == 1)
      lambda <- lambda/10
    if (crossprod(dl, dgam) < tol)
      break
    if (iter > maxit) {
      conv=0     ## not convergent
      warning("reml: Max iterations exceeded")
      break
    }
  }
  cov.gam <- chol2inv(chol(ZVZ))
  se.gam <- sqrt(diag(cov.gam))
  cov.beta <- chol2inv(qr.R(fitm$qr))
  se.beta <- sqrt(diag(cov.beta))
  list(beta = fitm$coef, se.beta = se.beta, gamma = gam, se.gam = se.gam,
       mu = fitm$fitted, phi = phi, deviance = dev, h = h, cov.beta = cov.beta,
       cov.gam = cov.gam, conv=conv)
}

###1. MGI data and Pilot data
mgi <- read.table("mgi_cont.txt", header = TRUE, sep = "\t", fill = TRUE)

mgi$sexc[mgi$sex=='Female'] <- 1
mgi$sexc[mgi$sex=='Male'] <- 2

mgi$racec[mgi$race=='Non-Hispanic White'] <- 1
mgi$racec[mgi$race=='Non-Hispanic Black'] <- 2
mgi$racec[mgi$race=='Other'] <- 3
mgi$racec[mgi$race=='Hispanic'] <- 4

mgi$smkc[mgi$smk=="Never"] <- 1
mgi$smkc[mgi$smk=="Former"] <- 2
mgi$smkc[mgi$smk=="Current"] <- 3


nh <- read.table("nhanes.txt", header = TRUE, sep = "\t", fill = TRUE)
nh$weight <- round(1/nh$pr,0)
nh$weight <- round(0.0016*(nh$weight),0)
nh$agec = data.table::between(nh$age, 18,39)+
  2*data.table::between(nh$age, 40,59)+
  3*data.table::between(nh$age, 60,150)

Target.sub <- nh[rep(seq.int(1,nrow(nh)), nh$weight), 1:10]

Target.sub$select <- 1
Target.sub$iftest <- NA
Target.sub$y <- NA
Target.sub$lambda1 <- NA
Target.sub$smkc[Target.sub$smk=="Never"] <- 1
Target.sub$smkc[Target.sub$smk=="Former"] <- 2
Target.sub$smkc[Target.sub$smk=="Current"] <- 3

Target.sub$e1 <- NA
Target.sub$e2 <- NA
setDT(Target.sub)
setkey(Target.sub, agec, race, sex)
Target.sub.dist <- Target.sub[, .N, keyby = list(agec, race,sex)]
Target.sub.dist$pr <- Target.sub.dist$N/dim(Target.sub)[1]
Target.sub.dist$N_mgi <- round(Target.sub.dist$pr*5000,0)

Target.sub$racec[Target.sub$race=='Non-Hispanic White'] <- 1
Target.sub$racec[Target.sub$race=='Non-Hispanic Black'] <- 2
Target.sub$racec[Target.sub$race=='Other'] <- 3
Target.sub$racec[Target.sub$race=='Hispanic'] <- 4

Target.sub$sexc[Target.sub$sex=='Female'] <- 1
Target.sub$sexc[Target.sub$sex=='Male'] <- 2

Target.sub = Target.sub[!is.na(Target.sub$bmi),]

# dat <- mgi
# budget=1000000
# n=330000000
# c0=10000
# c1=0.01
# npilot = 100

###2. Application function
my_application=function(dat, budget=500000, n=330000000, c0=10000, c1=0.01,npilot = 100){
  
  #pilot sample for approach 3
  idlength <- dim(dat)[1]
  pilot.ind <- sample(1:idlength, npilot)
  pilot <- dat[pilot.ind,]
  
  mgi.sample <- dat
  ne <- length(mgi.sample$lambda1)
  
  #estimate v.v.2
  Y <- pilot$y
  X <- model.matrix(~pilot$bmic+pilot$age+pilot$racec+pilot$sexc+pilot$smkc)
  
  Z = model.matrix(~pilot$bmic+pilot$age+pilot$racec+pilot$sexc+pilot$smkc)
  Zp = model.matrix(~mgi.sample$bmic+mgi.sample$age+mgi.sample$racec+mgi.sample$sexc+mgi.sample$smkc)
  Zp.2 = model.matrix(~Target.sub$bmic+Target.sub$age+Target.sub$racec+Target.sub$sexc+Target.sub$smkc)
  
  fitv = myremlscore(Y, X=X, Z=Z, maxit=10000)
  gamma = fitv$gamma

  v.v = exp(Zp%*%gamma) 
  v.v.2 = exp(Zp.2%*%gamma) 
  
  #Estimate Epp
  a <- 100
  b <- 5
  c2.v = a+b*mgi.sample$bmic+b*mgi.sample$age+b*mgi.sample$racec+b*mgi.sample$sexc+b*mgi.sample$smkc
  c2.v.2 = a+b*Target.sub$bmic+b*Target.sub$age+b*Target.sub$racec+b*Target.sub$sexc+b*Target.sub$smkc
  
  Epp = mean(sqrt(c2.v.2*v.v.2))
  
  #Calculate lambda2 for Theorem 1
  lambda1.v <- mgi.sample$lambda1
  Bprime <- budget - c0 - ne*c1
  lambda2.opt = pmin(sqrt(v.v/c2.v)*Bprime/(lambda1.v*n*Epp), 0.9999)#sqrt(v.v.2/c2.v)*Bprime/(lambda1.v*n*Epp)#
  lambda2.rd = pmin(Bprime/(lambda1.v*n*mean(c2.v)), 0.9999)
  
  #study sample
  mgi.sample$r2 = rbinom(ne,1,lambda2.opt)
  mgi.sample$r2.rd = rbinom(ne,1,lambda2.rd)
  sum(mgi.sample$r2)
  sum(mgi.sample$r2.rd)
  
  study <- mgi.sample[mgi.sample$r2==1,]
  study.rd <- mgi.sample[mgi.sample$r2.rd==1,]
  
  fite2 = glmrob(y~factor(bmic)+age+factor(racec)+factor(sex)+factor(smkc), data=study, family=binomial, control = glmrobMqle.control(maxit = 500))
  fite2.rd = glmrob(y~factor(bmic)+age+factor(racec)+factor(sex)+factor(smkc), data=study.rd, family=binomial, control = glmrobMqle.control(maxit = 500))
  
  mgi.sample$e2=predict(fite2, newdata=mgi.sample, type = 'response')
  mgi.sample$e2.rd=predict(fite2.rd, newdata=mgi.sample, type = 'response')
  
  fite3 = glm(e2 ~ age+factor(racec)+factor(sex), data=mgi.sample, family=binomial)
  fite3.rd = glm(e2.rd ~ age+factor(racec)+factor(sex), data=mgi.sample, family=binomial)
  
  mgi.sample$e3 = predict(fite3, newdata=mgi.sample, type = 'response')
  mgi.sample$e3.rd = predict(fite3.rd, newdata=mgi.sample, type = 'response')
  
  Target.sub$e3 = predict(fite3, newdata=Target.sub, type = 'response') 
  Target.sub$e3.rd = predict(fite3.rd, newdata=Target.sub, type = 'response') 
  
  #Estimate beta
  parta = (mgi.sample$r2)*(mgi.sample$y)/(mgi.sample$lambda1*lambda2.opt)
  parta.rd = (mgi.sample$r2.rd)*(mgi.sample$y)/(mgi.sample$lambda1*lambda2.rd)

  partb = (mgi.sample$r2-lambda2.opt)*(mgi.sample$e2)/(mgi.sample$lambda1*lambda2.opt)
  partb.rd = (mgi.sample$r2.rd-lambda2.rd)*(mgi.sample$e2.rd)/(mgi.sample$lambda1*lambda2.rd)

  partc = (mgi.sample$e3)/(mgi.sample$lambda1)
  partc.rd = (mgi.sample$e3.rd)/(mgi.sample$lambda1)
  
  beta.opt = (sum(parta)-sum(partb)-sum(partc)+sum(Target.sub$e3)*(n/nrow(Target.sub)))/n
  beta.rd = (sum(parta.rd)-sum(partb.rd)-sum(partc.rd)+sum(Target.sub$e3.rd)*(n/nrow(Target.sub)))/n
  beta.naive = mean(mgi$y[mgi.sample$r2.rd==1])
  beta <- c(beta.naive, beta.rd, beta.opt)
  return(beta)
}

###3. bootstrap for ci and re
index_length = 50000

##########b1
bootbeta <- matrix(ncol = 3, nrow = index_length)
for (i in 1:index_length) {
  print(i)
  set.seed(i)
  skip_to_next <- FALSE
  tryCatch({
    setTimeLimit(20)
    bootdat <- mgi[sample(1:length(mgi$lambda1), length(mgi$lambda1), replace = TRUE),]
    bootbeta_tmp <- my_application(dat = bootdat, npilot = 100, budget = 100000)
    bootbeta[i,] <- bootbeta_tmp
  },
  error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next } 
}

mean_beta=apply(bootbeta, 2, mean, na.rm=TRUE)
var_beta=apply(bootbeta, 2, var, na.rm=TRUE)
mean_beta
var_beta
var_beta/var_beta[2]


CI = c(mean_beta[3]-qnorm(0.975)*sqrt(var_beta[3]), mean_beta[3]+qnorm(0.975)*sqrt(var_beta[3]),
  mean_beta[4]-qnorm(0.975)*sqrt(var_beta[4]), mean_beta[4]+qnorm(0.975)*sqrt(var_beta[4]))

save(bootbeta, mean_beta, var_beta, CI, file = "app_finalscript5.RData")

load("app_finalscript.RData")
load("app_finalscript2.RData")
load("app_finalscript3.RData")
load("app_finalscript4.RData")
load("app_finalscript5.RData")

