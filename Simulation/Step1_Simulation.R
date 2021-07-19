####Final Simulation
#Final Main Update: July 18, 2021
#Guanghao Zhang and Xu Shi
#ghzhang@umich.edu

library(statmod)
library(robustbase)  
library(ggplot2)
library(reshape)
library(ggpubr)
library(MatchIt)

i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) #use array to 

## The following function was contributed by Xuesong Yu and Peter Gilbert, in Gilbert P, Yu X, Rotnitzky A (2013) 
## "Optimal Auxiliary-Covariate Based Sampling Design for 
## Semiparametric Efficient Estimation of a Mean or Mean Difference, with 
## Application to Clinical Trials"
## The function uses the remlscore function from the statmod library. 

myremlscore = function (y, X, Z, trace = FALSE, tol = 1e-05, maxit = 10000)
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

simulation = function(n, npilot, mn, sigma, alpha0, alpha1, alpha2, nu00, nu0, nu1, nu2, nu3, Bprime) {
  
  #data generation
  W0target = rnorm(n, mean=mn, sd=sigma)
  
  lambda1.true = plogis(W0target)
  #this is used to generate ne'<<ne
  #lambda1.true = ifelse(W0target>0.08, 0.9,0.1)
  #this is used to generate ne'<ne
  
  dt <- rbinom(n, 1, prob = lambda1.true)
  ne <- sum(dt)
  
  Bprime <- Bprime - 0.01*(sum(dt)-5000)
  W1target=rnorm(n, mean=mn, sd=sigma)
  Ytarget=rnorm(n, mean=alpha0 + alpha1*W0target + alpha2*W1target, sd=sqrt(exp(nu00+nu0*W0target+nu1*W0target^2+nu2*W1target+nu3*W1target^2)))
  Target <- data.frame(cbind(Y=Ytarget,W0=W0target,W1=W1target))
  
  R1=dt
  lambda.1.v = lambda1.true[R1==1]
  
  W0ehr=W0target[R1==1]
  W1ehr=W1target[R1==1]
  Yehr=Ytarget[R1==1]
  EHR <- data.frame(cbind(Y=Yehr,W0=W0ehr,W1=W1ehr,lambda1=lambda.1.v))
  
  #matching process for approach 4
  EHR.sub <- EHR
  EHR.sub$select <- 0
  Target.sub <- data.frame(cbind(Y=rnorm(nrow(EHR), mean=alpha0 + alpha1*W0target + alpha2*W1target, sd=sqrt(exp(nu00+nu0*W0target+nu1*W0target^2+nu2*W1target+nu3*W1target^2))),
                                 W0=rnorm(nrow(EHR), mean=mn, sd=sigma),
                                 W1=rnorm(nrow(EHR), mean=mn, sd=sigma)))
  Target.sub$select <- 1
  Target.sub$lambda1 <- 0
  dat.full <- rbind(EHR.sub, Target.sub)
  match.it <- matchit(select ~ W0,
                      data = dat.full, method = "nearest", ratio = 1, caliper = 0.05, replace = F) 
  EHR.new <- match.data(match.it, group = "control")[,1:ncol(EHR)] 
  W0ehr.new <- EHR.new$W0
  W1ehr.new <- EHR.new$W1
  lambda.1.v.new <- nrow(EHR.new)/n
  
  pilot.ind <- sample(1:ne, npilot)
  Pilot <- EHR[pilot.ind,]
  W0pilot <- Pilot$W0
  W1pilot <- Pilot$W1
  Ypilot <- Pilot$Y
  
  Pilot.new <- EHR.new[sample(1:nrow(EHR.new), npilot),]
  W0pilot.new <- Pilot.new$W0
  
  W1pilot.new <- Pilot.new$W1
  Ypilot.new <- Pilot.new$Y
  
  X=cbind(1, W0pilot, W1pilot)
  Z=cbind(1, W0pilot, W0pilot^2, W1pilot, W1pilot^2)
  X.new=cbind(1, W0pilot.new, W1pilot.new)
  Z.new=cbind(1, W0pilot.new, W0pilot.new^2, W1pilot.new, W1pilot.new^2)
  Zwrong=cbind(1, W0pilot, W1pilot)
  
  fit.var=myremlscore(Ypilot, X, Z, maxit=5000)
  gamma= fit.var$gamma
  nu00.hat=gamma[1,1]
  nu0.hat=gamma[2,1]
  nu1.hat=gamma[3,1]
  nu2.hat=gamma[4,1]
  nu3.hat=gamma[5,1]
  v.v.2=exp(nu00.hat+nu0.hat*W0ehr+nu1.hat*W0ehr^2+nu2.hat*W1ehr+nu3.hat*W1ehr^2) #variance model
  
  v.v.2.new=exp(nu00.hat+nu0.hat*W0ehr.new+nu1.hat*W0ehr.new^2+nu2.hat*W1ehr.new+nu3.hat*W1ehr.new^2)
  
  fit.var.wrong=myremlscore(Ypilot, X, Zwrong, maxit=5000)
  gamma.wrong= fit.var.wrong$gamma
  nu00.hat.wrong=gamma.wrong[1,1]
  nu0.hat.wrong=gamma.wrong[2,1]
  nu1.hat.wrong=gamma.wrong[3,1]
  v.v.2.wrong=exp(nu00.hat.wrong+nu0.hat.wrong*W0ehr+nu1.hat.wrong*W1ehr)
  
  v.v.2.true=exp(nu00+nu0*W0ehr+nu1*W0ehr^2+nu2*W1ehr+nu3*W1ehr^2)
  
  fit1 = glmrob(Y~W0, data=Pilot, family=gaussian, control = glmrobMqle.control(maxit = 500)) #mean model E(|W0)
  fit2 = glmrob(Y~W0+W1, data=Pilot, family=gaussian, control = glmrobMqle.control(maxit = 500)) #mean model E(|W1bar)
  
  fit1.wrong = glmrob(Y~W0+W0^2, data=Pilot, family=gaussian, control = glmrobMqle.control(maxit = 500))
  fit2.wrong = glmrob(Y~W0+W1+W0^2+W1^2, data=Pilot, family=gaussian, control = glmrobMqle.control(maxit = 500))
  
  Target$e1 = predict(fit1, newdata=Target)
  EHR$e1=predict(fit1, newdata=EHR)
  EHR$e2=predict(fit2, newdata=EHR)
  
  EHR.new$e1.new=predict(fit1, newdata=EHR.new)
  EHR.new$e2.new=predict(fit2, newdata=EHR.new)
  
  Target$e1.w = predict(fit1.wrong, newdata=Target)
  EHR$e1.w=predict(fit1.wrong, newdata=EHR)
  EHR$e2.w=predict(fit2.wrong, newdata=EHR)
  
  Target$e1.t = alpha0 + alpha1*W0target + alpha2*mn
  EHR$e1.t=alpha0 + alpha1*W0ehr + alpha2*mn
  EHR$e2.t=alpha0 + alpha1*W0ehr + alpha2*W1ehr
  
  c2.v = rep(100,ne) 
  c2.v.new = rep(100,nrow(EHR.new))
  
  Epp = mean(sqrt(c2.v*v.v.2))
  Epp.new = mean(sqrt(c2.v.new*v.v.2.new))
  EppW = mean(sqrt(c2.v*v.v.2.wrong))
  EppT = mean(sqrt(c2.v*v.v.2.true))
  
  lambda2.bar = Bprime/(lambda.1.v*n*mean(c2.v))
  lambda2.bar = ifelse(lambda2.bar>=1, 1, lambda2.bar)
  
  lambda2.opt = sqrt(v.v.2/c2.v)*Bprime/(lambda.1.v*n*Epp)
  lambda2.opt = ifelse(lambda2.opt>=1, 1, lambda2.opt)
  
  lambda2.opt.new = sqrt(v.v.2.new/c2.v.new)*(Bprime+0.01*(ne-nrow(EHR.new)))/(lambda.1.v.new*n*Epp.new)
  lambda2.opt.new = ifelse(lambda2.opt.new>=1, 1, lambda2.opt.new)
  
  lambda2.opt.wrong = sqrt(v.v.2.wrong/c2.v)*Bprime/(lambda.1.v*n*EppW)
  lambda2.opt.wrong = ifelse(lambda2.opt.wrong>=1, 1, lambda2.opt.wrong)
  
  lambda2.opt.true = sqrt(v.v.2.true/c2.v)*Bprime/(lambda.1.v*n*EppT)
  lambda2.opt.true = ifelse(lambda2.opt.true>=1, 1, lambda2.opt.true)
  
  Rsimple <- rbinom(ne,1,lambda2.bar)
  Ropt <- rbinom(ne,1,lambda2.opt)
  Ropt.new <- rbinom(nrow(EHR.new),1,lambda2.opt.new)
  RoptW <- rbinom(ne,1,lambda2.opt.wrong)
  RoptT <- rbinom(ne,1,lambda2.opt.true)
  
  EHR$lambda2bar <- lambda2.bar
  EHR$lambda2opt <- lambda2.opt
  EHR.new$lambda2opt.new <- lambda2.opt.new
  EHR$lambda2optW <- lambda2.opt.wrong
  EHR$lambda2optT <- lambda2.opt.true
  
  EHR$R2simple <- Rsimple
  EHR$R2opt <- Ropt
  EHR.new$R2opt.new <- Ropt.new
  EHR$R2optW <- RoptW
  EHR$R2optT <- RoptT
  
  obs.Y = Yehr[Rsimple==1]
  beta.hat.naive = mean(obs.Y)#(1)
  
  Uge = (EHR$e2-EHR$e1)/(EHR$lambda1)
  Uge.new = (EHR.new$e2.new-EHR.new$e1.new)/(lambda.1.v.new)
  UgeW = (EHR$e2.w-EHR$e1.w)/(EHR$lambda1)
  UgeT = (EHR$e2.t-EHR$e1.t)/(EHR$lambda1)
  
  Ugs.simple = (EHR$R2simple)*(EHR$Y-EHR$e2)/(EHR$lambda1*EHR$lambda2bar)
  Ugs = (EHR$R2opt)*(EHR$Y-EHR$e2)/(EHR$lambda1*EHR$lambda2opt)
  Ugs.new = (EHR.new$R2opt.new)*(EHR.new$Y-EHR.new$e2.new)/(lambda.1.v.new*EHR.new$lambda2opt.new)
  UgsERVW = (EHR$R2optW)*(EHR$Y-EHR$e2)/(EHR$lambda1*EHR$lambda2optW)
  UgsEWVR = (EHR$R2opt)*(EHR$Y-EHR$e2.w)/(EHR$lambda1*EHR$lambda2opt)
  UgsTT = (EHR$R2optT)*(EHR$Y-EHR$e2.t)/(EHR$lambda1*EHR$lambda2optT)
  
  beta.hat.simple= (sum(Target$e1)+sum(Uge)+sum(Ugs.simple))/n #(2)
  
  beta.hat.opt= (sum(Target$e1)+sum(Uge)+sum(Ugs))/n #(3-a) 
  beta.hat.optERVW= (sum(Target$e1)+sum(Uge)+sum(UgsERVW))/n#(3-b)
  beta.hat.optEWVR= (sum(Target$e1.w)+sum(UgeW)+sum(UgsEWVR))/n#(3-c)
  beta.hat.optTT= (sum(Target$e1.t)+sum(UgeT)+sum(UgsTT))/n#(3-d)
  beta.hat.opt.new= (sum(Target$e1)+sum(Uge.new)+sum(Ugs.new))/n #(3-e)
  
  all.beta=c(beta.hat.naive, beta.hat.simple, beta.hat.opt, beta.hat.optERVW, beta.hat.optEWVR, beta.hat.optTT, beta.hat.opt.new)
  return(all.beta)
}


nu0.finalseq <- c(0.97,0.82,-0.64)

for (j in 1:3) {
  set.seed(i)
  beta=simulation(n=10000,
                  npilot=200,
                  Bprime=49500,
                  mn=0.05,
                  sigma=sqrt(2),
                  alpha0=0.1,
                  alpha1=3,
                  alpha2=0.01,
                  nu00=-1.5,
                  nu0=nu0.finalseq[j],
                  nu1=0.2,
                  nu2=0.01,
                  nu3=0.01)
  save(beta, file=paste0("intermediate_simdata/sim_beta_bignp_", i, "_", j, ".Rdata"))
}

for (j in 1:3) {
  set.seed(i)
  beta=simulation(n=10000,
                  npilot=50,
                  Bprime=49500,
                  mn=0.05,
                  sigma=sqrt(2),
                  alpha0=0.1,
                  alpha1=3,
                  alpha2=0.01,
                  nu00=-1.5,
                  nu0=nu0.finalseq[j],
                  nu1=0.2,
                  nu2=0.01,
                  nu3=0.01)
  save(beta, file=paste0("intermediate_simdata/sim_beta_smallnp_", i, "_", j, ".Rdata"))
}
