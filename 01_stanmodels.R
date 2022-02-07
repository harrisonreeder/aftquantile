
#read in packages
sapply(c("survival","flexsurv","rstpm2","SemiCompRisks","SemiCompRisksPen","rstan","loo",
         "tidyverse","readxl","lubridate","tidylog","survminer","visdat","AFTTV"),
       function(x){if(!suppressWarnings(require(package=x,character.only = T)))
       {install.packages(pkgs=x);require(package=x,character.only = T)}})

source("00_import_rosmap_data.R")
source("BayesSurv_AFTtvstan_2022-01-20.R")
ROSMAPtemp <- "[path]"
####Data Prep####
rosmap_baseline_sub <- rosmap_baseline[complete.cases(rosmap_baseline[,c("race_fctwhite","msex","marital_bl_fctmarried", "education_15plus","apoe4_any")]),]

Y_diag <- as.matrix(cbind(rosmap_baseline_sub$time_ad_dx_mid, #Left "edge" of censoring
                          ifelse(rosmap_baseline_sub$ad_diagnosis==1, 
                                 rosmap_baseline_sub$time_ad_dx_mid,Inf), #Right "edge" of censoring
                          rosmap_baseline_sub$time_bl)) #enrollment time

Y_death <- as.matrix(cbind(rosmap_baseline_sub$time_last, #Left "edge" of censoring
                           ifelse(rosmap_baseline_sub$death==1, 
                                  rosmap_baseline_sub$time_last,Inf), #Right "edge" of censoring
                           rosmap_baseline_sub$time_bl))  #enrollment time

Y_diag_scale <- Y_diag/10
Y_death_scale <- Y_diag/10

Xmat <- as.matrix(rosmap_baseline_sub[,c("race_fctwhite","msex","marital_bl_fctmarried", "education_15plus","apoe4_any")])
Xmat_names <- colnames(Xmat)

Xmat_baselineonly <- as.matrix(rosmap_baseline_sub[,c("race_fctwhite","msex","marital_bl_fctmarried", "education_15plus")])
Xmat_baselineonly_names <- colnames(Xmat_baselineonly)

#for use when tv covariate is apoe, which is time-invariant
Xtv_vec <- rosmap_baseline_sub$apoe4_any

#for use when tv covariate is AD onset, which is time-varying
Xtv_time_diag <- rosmap_baseline_sub$time_ad_dx_mid


#### RUN STAN MODELS ####
stan_AFT_tvcov_piece_compiled <- stan_model("stan/AFT_all_tvcov_piece.stan")
stan_AFT_tvcov_spline_compiled <- stan_model("stan/AFT_all_tvcov_spline.stan")

stan_AFT_tbp_tvcov_piece_compiled <- stan_model("stan/AFT_tbp_tvcov_piece.stan")
stan_AFT_tbp_tvcov_spline_compiled <- stan_model("stan/AFT_tbp_tvcov_spline.stan")

#### AD ONSET ####
knots_piece <- c(0,7.5,15,22.5,30,Inf)
knots_piece_scale <- knots_piece / 10

#parametric models
for(effect_label in c("invariant", "piece","spline",
  NULL)){
  for(baseline_temp in c("weibull", "lognormal",
                       NULL)){
    baseline_label <- switch(baseline_temp, "weibull"="wb","lognormal"="ln")
    effect_temp <- switch(effect_label, "invariant"="piecewise","piece"="piecewise","spline"="spline")
    nP_tv_temp <- switch(effect_label, "invariant"=1,"piece"=NULL,"spline"=4)
    knots_temp <- if(effect_label=="piece") knots_piece else NULL
    control_temp <- if(effect_label=="spline") list(adapt_delta=0.999) else NULL
    test_fit <- BayesSurv_AFTtvstan(Y = Y_diag,Xmat = Xmat_baselineonly,
                                    Xtv = Xtv_vec,
                                    baseline = baseline_temp,
                                    prior_intercept = "none", m_intercept = 0, sd_intercept = 0,
                                    prior_scale = "gamma", a_scale=0.3, b_scale=0.05,
                                    prior_beta = "none", m_beta=0, sd_beta=10,
                                    prior_beta_tv = "none", m_beta_tv=0, sd_beta_tv=10,
                                    knots = knots_temp,
                                    tv_type = effect_temp, nP_tv = nP_tv_temp,
                                    n_chains=3, n_sample = 10000, n_warmup=2000,
                                    control = control_temp)
    saveRDS(test_fit, file = paste0(ROSMAPtemp,"ROSMAP_",baseline_label,"_apoe_",
                                              effect_label,"_2022-01-27.RDS"))
  }
}

#TBP prior with weibull baseline
J_temp <- 5
initfinvariant <- function() { list(intercept=3, w = rep(1/J_temp,J_temp),
                                scale_raw=2, beta=numeric(4),
                                beta_tv=as.array(0), #beta_tv=rep(0,5), 
                                dirichlet_alpha_param=as.array(1)) }
initfpiece <- function() { list(intercept=3, w = rep(1/J_temp,J_temp),
                            scale_raw=2, beta=numeric(4),
                            beta_tv=rep(0,5), dirichlet_alpha_param=as.array(1)) }
initfspline <- function() { list(intercept=3, w = rep(1/J_temp,J_temp),
                            scale_raw=2, beta=numeric(4),
                            beta_tv=rep(0,4), dirichlet_alpha_param=as.array(1)) }
for(effect_label in c("invariant", "piece", "spline",
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
                                              init = if(effect_label == "piece"){
                                                initfpiece 
                                              } else if(effect_label == "spline"){
                                                initfspline
                                              } else { initfinvariant },
                                        n_sample = 10000, n_warmup=2000)
    saveRDS(test_fit_tbp, file = paste0(ROSMAPtemp,"ROSMAP_",baseline_label,"_tbp_apoe_",
                                        effect_label,"_2022-01-20.RDS"))
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
                              prior_intercept = "none", m_intercept = 0, sd_intercept = 0,
                              prior_scale = "gamma", a_scale=0.3, b_scale=0.05,
                              prior_beta = "none", m_beta=0, sd_beta=10,
                              prior_beta_tv = "none", m_beta_tv=0, sd_beta_tv=10,
                              tv_type = effect_temp, nP_tv = nP_tv_temp, knots = knots_temp,
                              n_chains=3, n_sample = 10000, n_warmup=2000)
    saveRDS(test_fit, file = paste0(ROSMAPtemp,"ROSMAP_",baseline_label,"_deathtv_",
                                    effect_label,"_2022-01-27.RDS"))
  }
}

#TBP prior with weibull baseline
J_temp <- 5
initfinvariant <- function() { list(intercept=3.5, w = rep(1/J_temp,J_temp),
                                    scale_raw=3.15, beta=numeric(5),
                                    beta_tv=as.array(-0.75), #beta_tv=rep(0,5), 
                                    dirichlet_alpha_param=as.array(1)) }
initfpiece <- function() { list(intercept=3.5, w = rep(1/J_temp,J_temp),
                                scale_raw=3.15, beta=numeric(5),
                                # beta_tv=as.array(-1),
                                beta_tv=rep(-0.75,5), dirichlet_alpha_param=as.array(1)) }
initfspline <- function() { list(intercept=3.5, w = rep(1/J_temp,J_temp),
                                 scale_raw=3.15, beta=numeric(5),
                                 beta_tv=rep(-0.75,4), dirichlet_alpha_param=as.array(1)) }
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
                                        seed = 1, 
                                        init = if(effect_label == "piece"){
                                          initfpiece 
                                        } else if(effect_label == "spline"){
                                          initfspline
                                        } else { initfinvariant },
                                        n_sample = 10000, n_warmup=2000)
    saveRDS(test_fit_tbp, file = paste0(ROSMAPtemp,"ROSMAP_",baseline_label,"_tbp_deathtv_",
                                        effect_label,"_2022-01-20.RDS"))
  }
}
