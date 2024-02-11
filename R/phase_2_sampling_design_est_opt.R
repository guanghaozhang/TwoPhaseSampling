#' Estimate the mean of an expensive outcome from optimal two-phase sampling scheme
#' Description
#' @param y_pilot outcome from the pilot data
#' @param dat_pilot the entire pilot data consisting of both the outcome and the covariates
#' @param dat_target the target population simulated based on W_0, i.e., patient characteristics predictive of both the outcome and whether a patient seeks healthcare and thus shows up in the EHR system.
#' @param dat_ehr the entire EHR data
#' @param lambda1_ehr the estimated sampling probability from the target population to the EHR population
#' @param formula_mean_target_est an object of class "formula" (a symbolic description of the model to be fitted) for E(Y|W_0)
#' @param formula_mean_ehr_est an object of class "formula" (a symbolic description of the model to be fitted) for E(Y|W_0, W_1)
#' @param design_matrix_mean_pilot_var_est the pilot data design matrix for mean, i.e., E(Y|W_0, W_1, R_1 = 1),in a heteroscedastic regression model. This design matrix should follow \code{formula_mean_ehr_est}.
#' @param design_matrix_var_pilot_var_est the pilot data design matrix for the variance, i.e., Var(Y|W_0, W_1, R_1 = 1), in a heteroscedastic regression model. 
#' @param design_matrix_var_ehr_var_est the design matrix for predicting the conditional variance of the outcome, i.e., Var(Y|W_0, W_1, R_1 = 1), in the EHR data.
#' @param b_prime the total budget for selecting a study sample 
#' @param n sample size of the target population
#' @details
#' Suppose a study to estimate the mean or mean difference of an expensive outcome in a target population, with inexpensive auxiliary covariates predictive of the outcome available in patients' electronic health records (EHR).
#' We proposed to recruit a study sample from the EHR sample Thus patient recruitment using EHR data constitutes a two-phase sampling framework: the EHR sample is a subset of the target population, and the study sample is a subset
#' of the EHR sample. This two-phase sampling framework also ensembles a monotone missing data problem and therefore we applied existing missing data works to derive the efficient estimator for the mean outcome. We then derived 
#' the variance formula for the estimator, which is a function of the sampling probabilities from EHR sample to the study sample. We obtained the efficient sampling probabilities by minimizing the variance of the estimator.
#' @return A list containing:\tabular{ll}{
#'    \code{beta_hat} \tab estimated mean outcome \cr
#'    \tab \cr
#'    \code{lambda2} \tab A numeric vector of sampling probabilities for the phase-II study sample \cr
#' }
#' @examples
#' # example code
#' set.seed(1)
#' n = 10000
#'n_pilot = 200
#'mn = 0.05
#'sigma = sqrt(2)
#'alpha0 = 0.1
#'alpha1 = 3
#'alpha2 = 0.01 
#'nu00 = -1.5
#'nu0 = 0.97
#'nu1 = 0.2
#'nu2 = 0.01
#'nu3 = 0.01
#'b = 49500

#'#target population
#'w0_target = rnorm(n, mean = mn, sd = sigma)

#'lambda1_true = plogis(w0_target) #sampling probability from target population to EHR sample
#'#this is used to generate n_e'<<n_e
#'#to generate n_e'<n_e, use lambda1_true = ifelse(w0_target>0.08, 0.9,0.1)

#'r1 = rbinom(n, 1, prob = lambda1_true) #sampling process from target population to EHR sample

#'b_prime = b - 0.01*(sum(r1)-5000) #budge left to extract the study sample
#'w1_target = rnorm(n, mean = mn, sd = sigma)
#'y_target = rnorm(n, 
#'                 mean = alpha0 + alpha1*w0_target + alpha2*w1_target, 
#'                 sd = sqrt(exp(nu00 + nu0*w0_target + nu1*w0_target^2 + nu2*w1_target + nu3*w1_target^2)))
#'dat_target = data.frame(cbind(y = y_target,
#'                              w0 = w0_target,
#'                              w1 = w1_target))

#'#ehr sample
#'lambda1_ehr = lambda1_true[r1==1]

#'w0_ehr = w0_target[r1==1]
#'w1_ehr = w1_target[r1==1]
#'w0_ehr_sq = w0_ehr^2
#'w1_ehr_sq = w1_ehr^2
#'y_ehr = y_target[r1==1]
#'dat_ehr = data.frame(cbind(y = y_ehr,
#'                           w0 = w0_ehr,
#'                           w1 = w1_ehr,
#'                           lambda1 = lambda1_ehr))
#'design_matrix_var_ehr_var_est = model.matrix(~ w0_ehr + w0_ehr_sq + w1_ehr + w1_ehr_sq)
#'n_e = nrow(dat_ehr)

#'#pilot data
#'pilot_idx = sample(1:n_e, n_pilot)
#'dat_pilot = dat_ehr[pilot_idx,]
#'w0_pilot = dat_pilot$w0
#'w0_pilot_sq = dat_pilot$w0^2
#'w1_pilot = dat_pilot$w1
#'w1_pilot_sq = dat_pilot$w1^2
#'y_pilot = dat_pilot$y

#'formula_mean_target_est = 'y ~ w0'
#'formula_mean_ehr_est = 'y ~ w0 + w1'

#'design_matrix_mean_pilot_var_est = model.matrix(~ w0_pilot + w1_pilot)
#'design_matrix_var_pilot_var_est = model.matrix(~ w0_pilot + w0_pilot_sq + w1_pilot + w1_pilot_sq)

#'rslt = phase_2_sampling_design_est_opt(y_pilot = y_pilot, 
#'                                       dat_pilot = dat_pilot,
#'                                       dat_target = dat_target,
#'                                       dat_ehr = dat_ehr,
#'                                       lambda1_ehr = lambda1_ehr,
#'                                       formula_mean_target_est = formula_mean_target_est,
#'                                       formula_mean_ehr_est = formula_mean_ehr_est,
#'                                       design_matrix_mean_pilot_var_est = design_matrix_mean_pilot_var_est,
#'                                       design_matrix_var_pilot_var_est = design_matrix_var_pilot_var_est,
#'                                       design_matrix_var_ehr_var_est = design_matrix_var_ehr_var_est,
#'                                       b_prime = b_prime,
#'                                       n = n)
#' mean(y_target) #mean outcome in target population
#' mean(y_ehr) #mean outcome in EHR
#' rslt$beta_hat #final estimate is close to the mean outcome in target population
#'                                       
#' @export
#' @import robustbase
#' @import stats

phase_2_sampling_design_est_opt = function(y_pilot = y_pilot, 
                                           dat_pilot = dat_pilot,
                                           dat_target = dat_target,
                                           dat_ehr = dat_ehr,
                                           lambda1_ehr = lambda1_ehr,
                                           formula_mean_target_est = formula_mean_target_est,
                                           formula_mean_ehr_est = formula_mean_ehr_est,
                                           design_matrix_mean_pilot_var_est = design_matrix_mean_pilot_var_est,
                                           design_matrix_var_pilot_var_est = design_matrix_var_pilot_var_est,
                                           design_matrix_var_ehr_var_est = design_matrix_var_ehr_var_est,
                                           b_prime = b_prime,
                                           n = n){
  
  #estimate conditional variance
  fit_var = statmod::remlscore(y_pilot, design_matrix_mean_pilot_var_est, design_matrix_var_pilot_var_est)
  nu_hat = matrix(fit_var$gamma[,1], ncol = 1)
  var_hat = exp(design_matrix_var_ehr_var_est%*%nu_hat)
  
  #estimate outcome models
  fit_mean_w0 = glmrob(formula_mean_target_est, data = dat_pilot, family = "gaussian", control = glmrobMqle.control(maxit = 500)) #mean model E(y|w0)
  fit_mean_w1bar = glmrob(formula_mean_ehr_est, data = dat_pilot,family = "gaussian", control = glmrobMqle.control(maxit = 500)) #mean model E(y|w1bar)
  
  dat_target$e1 = predict(fit_mean_w0, newdata = dat_target)
  dat_ehr$e1 = predict(fit_mean_w0, newdata = dat_ehr)
  dat_ehr$e2 = predict(fit_mean_w1bar, newdata = dat_ehr)
  
  n_e = nrow(dat_ehr)
  c2 = rep(100,n_e) 
  
  Epp = mean(sqrt(c2*var_hat))
  
  lambda2 = sqrt(var_hat/c2)*b_prime/(lambda1_ehr*n*Epp)
  lambda2 = ifelse(lambda2 >= 1, 1, lambda2)
  
  dat_ehr$lambda2 <- lambda2
  dat_ehr$r2 <- rbinom(n_e, 1, lambda2)
  
  Uge = (dat_ehr$e2-dat_ehr$e1)/(dat_ehr$lambda1)
  
  Ugs = (dat_ehr$r2)*(dat_ehr$y-dat_ehr$e2)/(dat_ehr$lambda1*dat_ehr$lambda2)
  
  beta_hat = (sum(dat_target$e1)+sum(Uge)+sum(Ugs))/n
  
  return(list(beta_hat = beta_hat,
              lambda2 = lambda2))
}
