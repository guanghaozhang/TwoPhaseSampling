### Plot for application study
rm(list=ls())

library(statmod)
library(robustbase)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(gridExtra)
library(grid)
library(R.utils)
library(data.table)
library(sampling)
library(dplyr)

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

###2. Plot function
my_application_plot=function(dat, budget=50000, n=330000000, c0=10000, c1=0.01, npilot=100){
  
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
  lambda2.opt = pmin(sqrt(v.v/c2.v)*Bprime/(lambda1.v*n*Epp), 0.9999)
  return(lambda2.opt)
}

###3. Obtain lambda2 and smooth curves
x1 <- my_application_plot(dat = mgi, budget=100000)
x2 <- my_application_plot(dat = mgi, budget=1000000)

fit.np.age.b1.c1.v1 <- smooth.spline(x1~mgi$age, spar = 0.5)
fit.np.age.b3.c1.v1 <- smooth.spline(x2~mgi$age, spar = 0.5)
fit.np.bmi.b1.c1.v1 <- smooth.spline(x1~mgi$bmic, spar = 1)
fit.np.bmi.b3.c1.v1 <- smooth.spline(x2~mgi$bmic, spar = 1)

###A glance at the curves
# jpeg('app_fig_ageb1c1v1_line.jpg')
# plot(fit.np.age.b1.c1.v1, type = "l")
# dev.off()
# 
# jpeg('app_fig_ageb3c1v1_line.jpg')
# plot(fit.np.age.b3.c1.v1, type = "l")
# dev.off()
# 
# 
# jpeg('app_fig_bmib1c1v1_line.jpg')
# plot(fit.np.bmi.b1.c1.v1, type = "l")
# dev.off()
# 
# 
# jpeg('app_fig_bmib3c1v1_line.jpg')
# plot(fit.np.bmi.b3.c1.v1, type = "l")
# dev.off()

###4. Data for plot
plot_age_data <- data.frame(x=c(fit.np.age.b1.c1.v1$x, fit.np.age.b3.c1.v1$x), 
                            y=c(fit.np.age.b1.c1.v1$y, fit.np.age.b3.c1.v1$y),
                            B=c(rep("b1",length(fit.np.age.b1.c1.v1$y)), rep("b2", length(fit.np.age.b3.c1.v1$y))))

# id_bmidata_1 <- sort(sample(1:length(fit.np.bmi.b1.c1.v1$y), 500))
# id_bmidata_2 <- sort(sample(1:length(fit.np.bmi.b3.c1.v1$y), 500))

# plot_bmi_data <- data.frame(x=c(fit.np.bmi.b1.c1.v1$x[id_bmidata_1], fit.np.bmi.b3.c1.v1$x[id_bmidata_1]), 
#                             y=c(fit.np.bmi.b1.c1.v1$y[id_bmidata_2], fit.np.bmi.b3.c1.v1$y[id_bmidata_2]),
#                             B=c(rep("b1",length(id_bmidata_1)), rep("b2", length(id_bmidata_2))))
plot_bmi_data <- data.frame(x=c(fit.np.bmi.b1.c1.v1$x, fit.np.bmi.b3.c1.v1$x), 
                            y=c(fit.np.bmi.b1.c1.v1$y, fit.np.bmi.b3.c1.v1$y),
                            B=c(rep("b1",length(fit.np.bmi.b1.c1.v1$y)), rep("b2", length(fit.np.bmi.b3.c1.v1$y))))
save(plot_age_data, plot_age_data, file = "app_plot_final.RData")


plot_age <- ggplot(data = plot_age_data, aes(x=x, y=y, linetype=B)) +
  geom_line() + theme_classic() + xlab("Age") + ylab(expression(paste(lambda[2]^"*",(bar(W)[1])))) + 
  scale_linetype_manual(values = c(1,3),
                        labels = c(expression(paste("B"==10^5,", ", hat(beta)["3"]==0.38,", ", "95% CI for ",hat(beta)["3"]=="(0.35,0.41)",", ", hat(beta)["4"]==99,", ", "95% CI for ",hat(beta)["4"]=="(99,99)",", ", "RE"["3vs1"]=="0.90",", ","RE"["3vs2"]=="0.96",", ","RE"["4vs2"]=="99")),
                                   expression(paste("B"==10^6,", ", hat(beta)["3"]==0.39,", ", "95% CI for ",hat(beta)["3"]=="(0.30,0.48)",", ", hat(beta)["4"]==99,", ", "95% CI for ",hat(beta)["4"]=="(99,99)",", ", "RE"["3vs1"]=="0.87",", ","RE"["3vs2"]=="0.97",", ","RE"["4vs2"]=="99")))) + 
  labs(subtitle = expression(paste(lambda[2]^"*",(bar(W)[1])," vs. Age"))) +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  theme(legend.title=element_blank()) + 
  theme(legend.text.align = 0) + scale_y_continuous(labels = scales::scientific)
#theme(legend.text=element_text(size=8)) +
#scale_y_continuous(labels = scales::number_format(accuracy = 0.0001, decimal.mark = '.'))

plot_bmi <- ggplot(data = plot_bmi_data, aes(x=x, y=y, linetype=B)) +
  geom_line() + theme_classic() + xlab("BMI") + ylab(expression(paste(lambda[2]^"*",(bar(W)[1])))) + 
  scale_linetype_manual(values = c(1,3),
                        labels = c(expression(paste("B"==10^5,", ", hat(beta)["3"]==0.38,", ", "95% CI for ",hat(beta)["3"]=="(0.35,0.41)",", ", hat(beta)["4"]==99,", ", "95% CI for ",hat(beta)["4"]=="(99,99)",", ", "RE"["3vs1"]=="0.90",", ","RE"["3vs2"]=="0.96",", ","RE"["4vs2"]=="99")),
                                   expression(paste("B"==10^6,", ", hat(beta)["3"]==0.39,", ", "95% CI for ",hat(beta)["3"]=="(0.30,0.48)",", ", hat(beta)["4"]==99,", ", "95% CI for ",hat(beta)["4"]=="(99,99)",", ", "RE"["3vs1"]=="0.87",", ","RE"["3vs2"]=="0.97",", ","RE"["4vs2"]=="99")))) + 
  labs(subtitle = expression(paste(lambda[2]^"*",(bar(W)[1])," vs. BMI"))) +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
  theme(legend.title=element_blank()) + 
  theme(legend.text.align = 0) + scale_y_continuous(labels = scales::scientific)
#theme(legend.text=element_text(size=8)) +
#scale_y_continuous(labels = scales::number_format(accuracy = 0.0001, decimal.mark = '.'))

#function to extract the legend of a ggplot; source:
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


p2_legend2 <- get_legend(plot_bmi)

#sum_plot <- ggarrange(plot_age, plot_bmi, ncol = 2, nrow = 1, heights = c(4, 4), common.legend = TRUE, legend = "bottom")
sum_plot2 <- grid.arrange(arrangeGrob(plot_age + theme(legend.position="none"), 
                                      plot_bmi + theme(legend.position="none"), nrow=1), 
                          p2_legend2, 
                          nrow=2)

ggsave("app_final_final.pdf", width = 8, height = 4, sum_plot2)
