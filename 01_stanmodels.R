
#read in packages
sapply(c("survival","SemiCompRisks","rstan","loo","spBayesSurv",
         "tidyverse","readxl","lubridate","tidylog"),
       function(x){if(!suppressWarnings(require(package=x,character.only = T)))
       {install.packages(pkgs=x);require(package=x,character.only = T)}})

source("00_import_rosmap_data.R")
source("AFT_functions_spline.R")
source("BayesSurv_AFTtvstan.R")
ROSMAPtemp <- "[path]"

####Data Prep####
rosmap_baseline_sub <- rosmap_baseline[complete.cases(rosmap_baseline[,c("race_fctwhite","msex","marital_bl_fctmarried", "education_15plus","apoe4_any")]),]
rosmap_baseline_sub_ic <- rosmap_baseline_sub
rosmap_baseline_sub_ic$time_ad_dx <- ifelse(rosmap_baseline_sub$ad_diagnosis==1, 
                                            rosmap_baseline_sub$time_ad_dx,Inf)
rosmap_baseline_sub_ic$time_prev_ad_dx <- ifelse(rosmap_baseline_sub$ad_diagnosis==1, 
                                            rosmap_baseline_sub$time_prev_ad_dx,rosmap_baseline_sub$time_last)

Y_diag <- as.matrix(cbind(rosmap_baseline_sub$time_ad_dx_mid, #Left "edge" of censoring
                          ifelse(rosmap_baseline_sub$ad_diagnosis==1, 
                                 rosmap_baseline_sub$time_ad_dx_mid,Inf), #Right "edge" of censoring
                          rosmap_baseline_sub$time_bl)) #enrollment time

Y_death <- as.matrix(cbind(rosmap_baseline_sub$time_last, #Left "edge" of censoring
                           ifelse(rosmap_baseline_sub$death==1, 
                                  rosmap_baseline_sub$time_last,Inf), #Right "edge" of censoring
                           rosmap_baseline_sub$time_bl))  #enrollment time

Xmat <- as.matrix(rosmap_baseline_sub[,c("race_fctwhite","msex","marital_bl_fctmarried", "education_15plus","apoe4_any")])
Xmat_names <- colnames(Xmat)

Xmat_baselineonly <- as.matrix(rosmap_baseline_sub[,c("race_fctwhite","msex","marital_bl_fctmarried", "education_15plus")])
Xmat_baselineonly_names <- colnames(Xmat_baselineonly)

#for use when tv covariate is apoe, which is time-invariant
Xtv_vec <- rosmap_baseline_sub$apoe4_any

#for use when tv covariate is AD onset, which is time-varying
Xtv_time_diag <- rosmap_baseline_sub$time_ad_dx_mid

#### RUN STAN MODELS ####
stan_AFT_tvcov_piece_compiled <- stan_model(paste0(stanpath,"AFT_all_tvcov_piece.stan"))
stan_AFT_tvcov_spline_compiled <- stan_model(paste0(stanpath,"AFT_all_tvcov_spline.stan"))

stan_AFT_tbp_mvn_tvcov_piece_compiled <- stan_model(paste0(stanpath,"AFT_tbp_mvn_tvcov_piece.stan"))
stan_AFT_tbp_mvn_tvcov_spline_compiled <- stan_model(paste0(stanpath,"AFT_tbp_mvn_tvcov_spline.stan"))

#### AD ONSET ####
knots_piece <- c(0,7.5,15,22.5,30,Inf)

#parametric models
for(effect_label in c("invariant", "piece","spline",
  NULL)){
  for(baseline_temp in c("weibull", "lognormal",
                       NULL)){
    baseline_label <- switch(baseline_temp, "weibull"="wb","lognormal"="ln")
    effect_temp <- switch(effect_label, "invariant"="piecewise","piece"="piecewise","spline"="spline")
    nP_tv_temp <- switch(effect_label, "invariant"=1,"piece"=NULL,"spline"=4)
    knots_temp <- if(effect_label=="piece") knots_piece else NULL
    # control_temp <- NULL
    control_temp <- if(effect_label=="spline") list(adapt_delta=0.999) else NULL
    test_fit <- BayesSurv_AFTtvstan(Y = Y_diag,Xmat = Xmat_baselineonly,
                                        Xtv = Xtv_vec,
                                        baseline = baseline_temp,
                                        prior_centering = "none",
                                        prior_mean_centering = NULL, prior_prec_centering = NULL,
                                        prior_intercept = "none", m_intercept = 0, sd_intercept = 0,
                                        prior_scale = "gamma", a_scale=0.3, b_scale=0.05,
                                        prior_beta = "none", m_beta=0, sd_beta=10,
                                        prior_beta_tv = "none", m_beta_tv=0, sd_beta_tv=10,
                                        knots = knots_temp,
                                        tv_type = effect_temp, nP_tv = nP_tv_temp,
                                        tbp=FALSE,
                                        dirichlet_alpha_fixed = FALSE,
                                        dirichlet_alpha_data = 2,
                                        a_dirichlet_alpha = 1,
                                        b_dirichlet_alpha = 1,
                                        n_chains= 3, n_cores = 3, 
                                        seed = 1,
                                        n_sample = 10000, n_warmup=2000)    
    print(summary(test_fit$stan_fit,pars=test_fit$pars)$summary)
    saveRDS(test_fit, file = paste0(ROSMAPtemp,"ROSMAP_",baseline_label,"_apoe_",
                                              effect_label,".RDS"))
  }
}

##TBP model initialization fits

#use delta method to transform PH parameter estimates into corresponding ests for stan model
#intercept = log(kappa)/alpha, log(scale) = log(alpha), beta_aft = -beta_ph / alpha
get_wb_AFT <- function(fit){
  # browser()
  lkappa <- fit$estimate[1]
  lalpha <- fit$estimate[2]
  beta_ph <- fit$estimate[-(1:2)]
  ncov <- length(beta_ph)
  A <- rbind(
    c(-1/exp(lalpha),
      lkappa/exp(lalpha),
      numeric(ncov)),
    c(0,1,numeric(ncov)),
    cbind(0,0,diag(1/exp(lalpha),nrow=ncov,ncol=ncov))
  )
  vcovAFT_temp <- A %*% vcov(fit) %*% t(A)
  dimnames(vcovAFT_temp) <- NULL
  # rownames(vcovAFT_temp) <- colnames(vcovAFT_temp) <- c("intercept","lalpha",names(beta_ph))
  # vcovAFT <- vcovAFT_temp[c(1,(3:3+ncov),2),c(1,(3:3+ncov),2)]
  list(est=c(intercept=-lkappa/exp(lalpha),lalpha=lalpha,-beta_ph/exp(lalpha)),
       vcov=vcovAFT_temp)
}

freq_wb_diag <- 
  survreg(Surv(time = time_ad_dx_mid,event = ad_diagnosis) ~ race_fctwhite + msex + marital_bl_fctmarried + education_15plus + apoe4_any, 
          data = rosmap_baseline_sub, dist = "weibull")
#to fit a parametric model accounting for left truncation, use external package
#devtools::install_github("https://github.com/harrisonreeder/SemiCompRisksFreq")
freq_wb_diag_lt <- 
  SemiCompRisksFreq::FreqSurv_HReg2(time_bl | time_ad_dx_mid + ad_diagnosis ~ race_fctwhite + msex + marital_bl_fctmarried + education_15plus + apoe4_any, 
                                    data = rosmap_baseline_sub,hazard = "wb")

prior_mean_centering <- get_wb_AFT(freq_wb_diag_lt)$est[1:2]
freq_beta <- get_wb_AFT(freq_wb_diag_lt)$est[-(1:2)]
prior_prec_centering <- MASS::ginv(get_wb_AFT(freq_wb_diag_lt)$vcov[1:2,1:2])/10
init_func <- function(){
  list(centering=prior_mean_centering, 
       beta=freq_beta[1:4], #excludes covariate with potential tv effect, apoe (HARDCODED!!)
       w = rep(1/5,5), #hardcoding that J will be 5
       dirichlet_alpha_param=as.array(1))
}

#TBP prior with weibull baseline
J_temp <- 5
for(effect_label in c("invariant","piece","spline",
  NULL)){
  for(baseline_temp in c("weibull",
                         NULL)){
    baseline_label <- switch(baseline_temp, "weibull"="wb","lognormal"="ln")
    effect_temp <- switch(effect_label, "invariant"="piecewise","piece"="piecewise","spline"="spline")
    nP_tv_temp <- switch(effect_label, "invariant"=1,"piece"=NULL,"spline"=4)
    knots_temp <- if(effect_label=="piece") knots_piece else NULL
    control_temp <- list(adapt_delta=0.9)
    n_chains_temp <- 3
    seed_temp <- if(effect_label=="spline") 1 else 3
    test_fit_tbp <- BayesSurv_AFTtvstan(Y = Y_diag,Xmat = Xmat_baselineonly,
                                        Xtv = Xtv_vec,
                                        baseline = baseline_temp,
                                        prior_centering = "mvn",
                                        prior_mean_centering = prior_mean_centering,
                                        prior_prec_centering = prior_prec_centering,
                                        prior_intercept = "none", m_intercept = 0, sd_intercept = 0,
                                        prior_scale = "gamma", a_scale=0.3, b_scale=0.05,
                                        prior_beta = "none", m_beta=0, sd_beta=10,
                                        prior_beta_tv = "none", m_beta_tv=0, sd_beta_tv=10,
                                        knots = knots_temp,
                                        tv_type = effect_temp, nP_tv = nP_tv_temp,
                                        tbp=TRUE, J = J_temp,
                                        dirichlet_alpha_fixed = FALSE,
                                        dirichlet_alpha_data = 2,
                                        a_dirichlet_alpha = 1,
                                        b_dirichlet_alpha = 1,
                                        n_chains= n_chains_temp,n_cores = n_chains_temp, 
                                        seed = seed_temp, 
                                        init = init_func,
                                      n_sample = 10000, n_warmup=2000)
    print(summary(test_fit_tbp$stan_fit,pars=test_fit_tbp$pars)$summary)
    saveRDS(test_fit_tbp, file = paste0(ROSMAPtemp,"ROSMAP_",baseline_label,"_tbp_apoe_",
                                        effect_label,".RDS"))
  }
}

#### Death ####
knots_piece_death <- c(0,1,3,5,10,Inf)
knots_piece_death_scale <- c(0,1,3,5,10,Inf) / 10

#parametric models
for(baseline_temp in c('weibull','lognormal')){
  for(effect_label in c("invariant","piece","spline")){
    baseline_label <- switch(baseline_temp, "weibull"="wb","lognormal"="ln")
    effect_temp <- switch(effect_label, "invariant"="piecewise","piece"="piecewise","spline"="spline")
    nP_tv_temp <- switch(effect_label, "invariant"=1,"piece"=NULL,"spline"=4)
    knots_temp <- if(effect_label=="piece") knots_piece_death else NULL
    test_fit <- BayesSurv_AFTtvstan_tvcov(Y = Y_death,Xmat = Xmat,
                              Xtv_time = Xtv_time_diag, #ignored when tv_type="none"
                              baseline = baseline_temp,
                              prior_centering = "none",
                              prior_mean_centering = prior_mean_centering,
                              prior_prec_centering = prior_prec_centering,
                              prior_intercept = "none", m_intercept = 0, sd_intercept = 0,
                              prior_scale = "gamma", a_scale=0.3, b_scale=0.05,
                              prior_beta = "none", m_beta=0, sd_beta=10,
                              prior_beta_tv = "none", m_beta_tv=0, sd_beta_tv=10,
                              tv_type = effect_temp, nP_tv = nP_tv_temp, knots = knots_temp,
                              n_chains=3, n_sample = 10000, n_warmup=2000)
                              # n_chains=4, n_sample = 5000, n_warmup=3000)
    print(summary(test_fit$stan_fit,pars=test_fit$pars)$summary)
    saveRDS(test_fit, file = paste0(ROSMAPtemp,"ROSMAP_",baseline_label,"_deathtv_",
                                    effect_label,".RDS"))
  }
}

##TBP Initialization frequentist fit
#record-split to accommodate time-varying covariate
rosmap_baseline_sub_long <- 
  tmerge(data1=rosmap_baseline_sub,
         data2=rosmap_baseline_sub, #set data as both "data1", "data2"
         id=studyid, #"subject" is our identifier variable for each person
         death = event(time_last, death), #define "death" by time2 and event2 variables
         ad_tv_ind=tdc(time_ad_dx)) #create a variable "x1" for surgery status from time1 variable
rosmap_baseline_sub_long$tstart[rosmap_baseline_sub_long$tstart==0] <- rosmap_baseline_sub_long$time_bl[rosmap_baseline_sub_long$tstart==0]

form_death_cox_lhs_string <- "Surv(time=tstart,time2=tstop,event=death,type = 'counting')"
form_rhs_string <- paste0(Xmat_names,collapse = " + ")
form_death_cox <- as.formula(paste(form_death_cox_lhs_string, "~", form_rhs_string, "+ ad_tv_ind"))

freq_wb_death_lt <- 
  SemiCompRisksFreq::FreqSurv_HReg2(tstart | tstop + death ~ 
        race_fctwhite + msex + marital_bl_fctmarried + education_15plus + apoe4_any + ad_tv_ind, 
    data = rosmap_baseline_sub_long,hazard = "wb")

prior_mean_centering <- get_wb_AFT(freq_wb_death_lt)$est[1:2]
freq_beta <- get_wb_AFT(freq_wb_death_lt)$est[-(1:2)]
prior_prec_centering <- MASS::ginv(get_wb_AFT(freq_wb_death_lt)$vcov[1:2,1:2])/10
init_death_func <- function(){
  list(centering=prior_mean_centering, 
       beta=freq_beta[1:5], #excludes covariate with potential tv effect, ad_tv (HARDCODED!!)
       w = rep(1/5,5), #hardcoding that J will be 5
       dirichlet_alpha_param=as.array(1))
}

#TBP prior with weibull baseline
J_temp <- 5
for(effect_label in c("invariant", 
                      "piece",
                      "spline",
                      NULL)){
  for(baseline_temp in c("weibull",
                         NULL)){
    baseline_label <- switch(baseline_temp, "weibull"="wb","lognormal"="ln")
    effect_temp <- switch(effect_label, "invariant"="piecewise","piece"="piecewise","spline"="spline")
    nP_tv_temp <- switch(effect_label, "invariant"=1,"piece"=NULL,"spline"=4)
    knots_temp <- if(effect_label=="piece") knots_piece_death else NULL
    control_temp <- list(adapt_delta=0.9)
    n_chains_temp <- 3
    test_fit_tbp <- BayesSurv_AFTtvstan_tvcov(Y = Y_death,Xmat = Xmat,
                                        Xtv_time = Xtv_time_diag,
                                        baseline = baseline_temp,
                                        prior_centering = "mvn",
                                        prior_mean_centering = prior_mean_centering,
                                        prior_prec_centering = prior_prec_centering,
                                        prior_intercept = "none", m_intercept = 0, sd_intercept = 0,
                                        prior_scale = "gamma", a_scale=0.3, b_scale=0.05,
                                        prior_beta = "none", m_beta=0, sd_beta=10,
                                        prior_beta_tv = "none", m_beta_tv=0, sd_beta_tv=10,
                                        knots = knots_temp,
                                        tv_type = effect_temp, nP_tv = nP_tv_temp,
                                        tbp=TRUE, J = J_temp,
                                        dirichlet_alpha_fixed = FALSE,
                                        dirichlet_alpha_data = 2,
                                        a_dirichlet_alpha = 1,
                                        b_dirichlet_alpha = 1,
                                        n_chains= n_chains_temp,n_cores = n_chains_temp, 
                                        seed = if(effect_label=="piece") 2 else 1, 
                                        init = init_death_func,
                                        n_sample = 10000, n_warmup=2000)
    print(summary(test_fit_tbp$stan_fit,pars=test_fit_tbp$pars)$summary)
    check_hmc_diagnostics(test_fit_tbp$stan_fit)
    print(plot(test_fit_tbp$stan_fit,plotfun="trace",pars=c("w","intercept","scale")))
    saveRDS(test_fit_tbp, file = paste0(ROSMAPtemp,"ROSMAP_",baseline_label,"_tbp_deathtv_",
                                        effect_label,".RDS"))
  }
}
