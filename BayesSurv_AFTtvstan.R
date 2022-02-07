
source("/Users/reederh/Dropbox/Harrison Reeder/Paper2_AFT/Code/AFT_functions_spline_2022-01-20.R")
# source("/Users/reederh/Dropbox/Harrison Reeder/Paper2_AFT/Code/AFT_functions_tvcov_2021-11-03.R")


#' Bayesian Accelerated Failure Time Model with Percentile-varying Effects
#'
#' @param Y an n by 3 matrix with columns: left end of interval, right end of interval, left truncation time
#' @param Xmat an n by p matrix of covariates
#' @param Xvec_tv an n-length vector with the covariate having a time-varying effect
#' @param Xvec_tv_time an n-length vector with the covariate having a time-varying effect
#' @param prior_list a list of priors being used
#' @param tuning_vec a vec of tuning parameters
#' @param hyper_list a list of hyperparameters
#' @param knots a vector of knots (doesn't start with 0, does end with maximum time)
#' @param numReps numeric total number of iterations
#' @param thin numeric number of iterations to take before each sample
#' @param burninPerc numeric proportion of numReps to toss as burn-in
#' @param n_chains numeric number of chains to initialize, if start_list is null
#' @param start_list a list of lists of start values. If null, n_chains sets number of chains
#'
#' @return a list with outputs
#' @export
BayesSurv_AFTtvstan <- function(Y, Xmat=NULL, Xtv=NULL, baseline, 
                                prior_intercept, m_intercept, sd_intercept,
                                prior_scale, a_scale, b_scale,
                                prior_beta, m_beta, sd_beta,
                                prior_beta_tv, m_beta_tv, sd_beta_tv,
                                tv_type, nP_tv, knots=NULL, 
                                tbp=FALSE, J=10,
                                dirichlet_alpha_fixed=TRUE,
                                a_dirichlet_alpha=NULL, b_dirichlet_alpha=NULL, 
                                dirichlet_alpha_data=1,
                                n_sample=1000, n_warmup=1000, n_chains=4, n_cores=4,
                                init = "random", seed = 1, control=NULL){
  # browser()

  N <- nrow(Y)
  if(!is.null(Xmat)){
    nP <- NCOL(Xmat) #number of time-invariant covariates with time-invariant effects
  } else{
    nP <- 0
    Xmat <- matrix(data=0,nrow=N,ncol=0)
  }
  
  N_o <- sum(Y[,1]==Y[,2])
  N_m <- nrow(Y) - N_o
  N_l <- sum(Y[,3]>0)

  #take prior scale as either character and convert, or as numeric directly.
  if(tolower(prior_scale) %in% c("ig","invgamma","inv-gamma")){
    prior_scale_num <- 1
  } else if(tolower(prior_scale) %in% c("g","gamma")){
    prior_scale_num <- 2
  } else{
    prior_scale_num <- 0
  }
  prior_intercept_num <- if(tolower(prior_intercept) %in% c("n","normal")) 1 else 0
  prior_beta_num <- if(tolower(prior_beta) %in% c("n","normal")) 1 else 0
  prior_beta_tv_num <- if(tolower(prior_beta_tv) %in% c("n","normal")) 1 else 0
  
  d_list <- list(N_o = N_o, N_m = N_m, N_l = N_l,
                 nP = nP, # number of betas
                 nP_base=nP, # number of betas
                 # data for censored subjects
                 y_o=Y[Y[,1]==Y[,2],1],
                 X_o=Xmat[Y[,1]==Y[,2],],
                 Xtv_o=Xtv[Y[,1]==Y[,2]],
                 Xtv_time_o=numeric(N_o),
                 # data for uncensored subjects
                 y_m=Y[Y[,1]!=Y[,2],1],
                 X_m=Xmat[Y[,1]!=Y[,2],],
                 Xtv_m=Xtv[Y[,1]!=Y[,2]],
                 Xtv_time_m=numeric(N_m),
                 #placeholder data for left-truncated subjects
                 L_l=Y[Y[,3]>0,3],
                 X_l=Xmat[Y[,3]>0,],
                 Xtv_l=Xtv[Y[,3]>0],
                 Xtv_time_l=numeric(N_l),
                 L_ind_o=rep(-1,N_o),
                 L_ind_m=rep(-1,N_m),
                 #intercept hyperparameters
                 prior_intercept = prior_intercept_num,
                 m_intercept = if(prior_intercept_num == 0) 1 else m_intercept,
                 sd_intercept = if(prior_intercept_num == 0) 1 else sd_intercept,
                 #scale hyperparameters
                 prior_scale = prior_scale_num,
                 a_scale=if(prior_scale_num == 0) 1 else a_scale,
                 b_scale=if(prior_scale_num == 0) 1 else b_scale,
                 #beta hyperparameters
                 prior_beta = prior_beta_num,
                 m_beta = if(prior_beta_num == 0) 1 else m_beta,
                 sd_beta = if(prior_beta_num == 0) 1 else sd_beta,
                 #beta_tv hyperparameters
                 prior_beta_tv = prior_beta_tv_num,
                 m_beta_tv = if(prior_beta_tv_num == 0) 1 else m_beta_tv,
                 sd_beta_tv=if(prior_beta_tv_num == 0) 1 else sd_beta_tv,
                 #baseline specification
                 weib_ind = as.numeric(tolower(baseline) %in% c("wb","weibull")),
                 #for use with tbp baseline
                 J=J,
                 dirichlet_alpha_fixed=as.numeric(dirichlet_alpha_fixed),
                 dirichlet_alpha_data=if(dirichlet_alpha_fixed) 
                   array(dirichlet_alpha_data,1) else numeric(0),
                 a_dirichlet_alpha=if(dirichlet_alpha_fixed) 1 else a_dirichlet_alpha,
                 b_dirichlet_alpha=if(dirichlet_alpha_fixed) 1 else b_dirichlet_alpha,
                 NULL
  )
  
  
  #Now, fix up left truncated times
  if(N_l > 0){
    #make a lookup table linking which observations are 
    #censored/uncensored with which are truncated/not
    #col 1 is the index of all the original data
    #col 2 is the index among the observed data
    #col 3 is the index among the censored data
    #col 4 is the corresponding index for the truncated data
    Y_ind_mat <- matrix(data = -1,nrow = nrow(Y),ncol=4)
    Y_ind_mat[,1] <- 1:nrow(Y)
    Y_ind_mat[Y[,1]==Y[,2],2] <- 1:sum(Y[,1]==Y[,2])
    Y_ind_mat[Y[,1]!=Y[,2],3] <- 1:sum(Y[,1]!=Y[,2])
    Y_ind_mat[Y[,3]>0,4] <- 1:sum(Y[,3]>0)
    d_list[["L_ind_o"]] <- Y_ind_mat[Y[,1]==Y[,2],4]
    d_list[["L_ind_m"]] <- Y_ind_mat[Y[,1]!=Y[,2],4]
    
  }
  
  #now, deal with time-varying effects
  if(!(tolower(tv_type) %in% c("pw","piecewise","spline"))){
    nP_tv <- 0
  } else{

    #now, set up knots
    if(is.null(knots)){
      stopifnot(!is.null(nP_tv))
      if(tolower(tv_type) %in% c("pw","piecewise")){
        knots <- c(0, quantile(Y[Y[,1]==Y[,2],1], 
                               probs = seq(from=0,to=1,length.out = nP_tv+1))[-c(1,nP_tv+1)])
      } else if(tolower(tv_type) %in% c("spline")){
        #get knots according to the time after covariate jump
        knots<- get_default_knots_vec(y = Y[,1], delta = as.numeric(Y[,1]==Y[,2]), 
                                      p_tv = nP_tv, tv_type="rp")
      }
    }
    if(knots[1] != 0){ knots <- c(0,knots) }
    if(knots[length(knots)] < 1e8){ knots <- c(knots, Inf) }
    d_list[["knots_tv"]] <- knots
    
    #now, generate the bases here
    if(tolower(tv_type) %in% c("pw","piecewise")){
      basis <- pw_cum_mat(y=Y[,1],
                                 knots=knots[-length(knots)], #remove the Inf knot
                                 intercept = TRUE) #note we include intercept here!
      d_list[["basis_piece_o"]] <- basis[Y[,1]==Y[,2],,drop=FALSE]
      d_list[["basis_piece_m"]] <- basis[Y[,1]!=Y[,2],,drop=FALSE]
      
      d_list[["basis_piece_l"]] <- matrix(data=0, nrow=0, ncol=NCOL(basis))
      if(N_l > 0){
        d_list[["basis_piece_l"]] <- pw_cum_mat(y=Y[Y[,3]>0,3],
                                                       knots=knots[-length(knots)],
                                                       intercept = TRUE)
      }
      
    } else if(tolower(tv_type) %in% c("spline")){
      basis <- get_basis_tv(x=Y[,1], tv_type="rp",
                            knots=knots[-c(1,length(knots))], #remove the Inf knot
                            deriv = FALSE, intercept=TRUE, #note we include intercept now!
                            flexsurv_compatible = FALSE)
      dbasis <- get_basis_tv(x=Y[,1], tv_type="rp",
                             knots=knots[-c(1,length(knots))], #remove the Inf knot
                             deriv = TRUE, intercept=TRUE, #note we include intercept now!
                             flexsurv_compatible = FALSE)
      d_list[["basis_spline_o"]] <- basis[Y[,1]==Y[,2],,drop=FALSE]
      d_list[["dbasis_spline_o"]] <- dbasis[Y[,1]==Y[,2],,drop=FALSE]
      d_list[["basis_spline_m"]] <- basis[Y[,1]!=Y[,2],,drop=FALSE]
      
      d_list[["basis_spline_l"]] <- matrix(data=0, nrow=0, ncol=NCOL(basis))
      if(N_l > 0){
        d_list[["basis_spline_l"]] <- get_basis_tv(x=Y[Y[,3]>0,3], tv_type="rp",
                                                   knots=knots[-c(1,length(knots))], #remove the Inf knot
                                                   deriv = FALSE, intercept=TRUE, #note we include intercept now!
                                                   flexsurv_compatible = FALSE)
      }
      
      #pass in a complete basis spline to try QR decomposition.
      d_list[["basis_spline_full"]] <- rbind(d_list[["basis_spline_o"]],
                                             d_list[["basis_spline_m"]],
                                             d_list[["basis_spline_l"]])  
    }
    nP_tv <- NCOL(basis)
    d_list[["nP_tv_piece"]] <- d_list[["nP_tv_spline"]] <- nP_tv

  }
  
  
  par <- c("intercept","scale","beta")
  if(tolower(tv_type) %in% c("pw","piecewise")){
    par <- c(par,"beta_tv")
    if(tbp){
      par <- c(par,"w",if(!dirichlet_alpha_fixed) "dirichlet_alpha_param" else NULL)
      stan_fit <- rstan::sampling(stan_AFT_tbp_tvcov_piece_compiled,data=d_list,
                                  iter=n_warmup + n_sample,warmup=n_warmup,
                                  chains=n_chains, cores=n_cores, init = init, control = control,
                                  seed = seed, pars=c(par,"log_lik"), include=TRUE)
    } else{
      stan_fit <- rstan::sampling(stan_AFT_tvcov_piece_compiled,data=d_list,
                                  iter=n_warmup + n_sample,warmup=n_warmup,
                                  chains=n_chains,cores=n_cores, init = init, control = control,
                                  seed = seed, pars=c(par,"log_lik"), include=TRUE)
    }
  } else if(tolower(tv_type) %in% c("spline")){
    par <- c(par,#"R_beta_tv",
             "beta_tv")
    if(tbp){
      par <- c(par,"w",if(!dirichlet_alpha_fixed) "dirichlet_alpha_param" else NULL)
      stan_fit <- rstan::sampling(stan_AFT_tbp_tvcov_spline_compiled,data=d_list,
                                  iter=n_warmup + n_sample,warmup=n_warmup,
                                  chains=n_chains,cores=n_cores, init = init, control = control,
                                  seed = seed, pars=c(par,"log_lik"), include=TRUE)
    } else{
      stan_fit <- rstan::sampling(stan_AFT_tvcov_spline_compiled,data=d_list,
                                  iter=n_warmup + n_sample,warmup=n_warmup,
                                  chains=n_chains,cores=n_cores, init = init, control = control,
                                  seed = seed, pars=c(par,"log_lik"), include=TRUE)
    }
  } else{
    if(tbp){
      par <- c(par,"w",if(!dirichlet_alpha_fixed) "dirichlet_alpha_param" else NULL)
      stan_fit <- rstan::sampling(stan_AFT_tbp_invariant_compiled,data=d_list,
                                  iter=n_warmup + n_sample,warmup=n_warmup,
                                  chains=n_chains,cores=n_cores, init = init, control = control,
                                  seed = seed, pars=c(par,"log_lik"), include=TRUE)
    } else{
      stan_fit <- rstan::sampling(stan_AFT_invariant_compiled,data=d_list,
                                  iter=n_warmup + n_sample,warmup=n_warmup,
                                  chains=n_chains,cores=n_cores, init = init, control = control,
                                  seed = seed, pars=c(par,"log_lik"), include=TRUE)
    }
  }

  outlist <- list(
    stan_fit=stan_fit,
    pars=par,
    nP=nP,
    prior_intercept=prior_intercept,
    m_intercept=m_intercept,
    sd_intercept=sd_intercept,
    prior_scale=prior_scale,
    a_scale=a_scale,
    b_scale=b_scale,
    prior_beta=prior_beta,
    m_beta=m_beta,
    sd_beta=sd_beta,
    prior_beta_tv=prior_beta_tv,
    m_beta_tv=m_beta_tv,
    sd_beta_tv=sd_beta_tv,
    baseline=baseline,
    knots=knots,
    nP_tv=nP_tv,
    tv_type=tv_type,
    tbp=tbp, 
    J=if(tbp) J else NULL,
    dirichlet_alpha_fixed=dirichlet_alpha_fixed,
    a_dirichlet_alpha=a_dirichlet_alpha, 
    b_dirichlet_alpha=b_dirichlet_alpha,
    dirichlet_alpha_data=dirichlet_alpha_data,
    init=init, seed = seed, 
    NULL
  )
}

BayesSurv_AFTtvstan_tvcov <- function(Y, Xmat=NULL, Xtv=NULL, Xtv_time,baseline,
                                prior_intercept, m_intercept, sd_intercept,
                                prior_scale, a_scale, b_scale,
                                prior_beta, m_beta, sd_beta,
                                prior_beta_tv, m_beta_tv, sd_beta_tv,
                                tv_type, nP_tv, knots=NULL,
                                tbp=FALSE, J=10,
                                dirichlet_alpha_fixed=TRUE,
                                a_dirichlet_alpha=NULL, b_dirichlet_alpha=NULL, 
                                dirichlet_alpha_data=1,
                                n_sample=1000, n_warmup=1000, n_chains=4, n_cores=4,
                                init = "random", seed = 1, control=NULL){
  # browser()
  
  N <- nrow(Y)
  nP <- NCOL(Xmat) #number of time-invariant covariates with time-invariant effects
  
  #assume that default time-varying covariate is a binary jump
  if(is.null(Xtv)){
    Xtv <- rep(1,N)
  }
  
  N_o <- sum(Y[,1]==Y[,2])
  N_m <- nrow(Y) - N_o
  N_l <- sum(Y[,3]>0)
  # stopifnot(N_l==0)
  
  #if there are any Xtv_times that are less than the time of study entry, that's bad.
  stopifnot(all(Y[,3] <= Xtv_time))
  
  #take prior scale as either character and convert, or as numeric directly.
  if(tolower(prior_scale) %in% c("ig","invgamma","inv-gamma")){
    prior_scale_num <- 1
  } else if(tolower(prior_scale) %in% c("g","gamma")){
    prior_scale_num <- 2
  } else{
    prior_scale_num <- 0
    a_scale <- 1 #placeholder
    b_scale <- 1 #placeholder
  }
  prior_intercept_num <- if(tolower(prior_intercept) %in% c("n","normal")) 1 else 0
  prior_beta_num <- if(tolower(prior_beta) %in% c("n","normal")) 1 else 0
  prior_beta_tv_num <- if(tolower(prior_beta_tv) %in% c("n","normal")) 1 else 0
  
  d_list <- list(N_o = N_o, N_m = N_m, N_l = N_l,
                 nP = nP, # number of betas
                 nP_base=nP, # number of betas
                 # data for censored subjects
                 y_o=Y[Y[,1]==Y[,2],1],
                 X_o=Xmat[Y[,1]==Y[,2],],
                 Xtv_o=Xtv[Y[,1]==Y[,2]],
                 Xtv_time_o=Xtv_time[Y[,1]==Y[,2]],
                 # data for uncensored subjects
                 y_m=Y[Y[,1]!=Y[,2],1],
                 X_m=Xmat[Y[,1]!=Y[,2],],
                 Xtv_m=Xtv[Y[,1]!=Y[,2]],
                 Xtv_time_m=Xtv_time[Y[,1]!=Y[,2]],
                 #placeholder data for left-truncated subjects (not permitted currently for tv covariates)
                 L_l=Y[Y[,3]>0,3], #the left truncation times
                 X_l=Xmat[Y[,3]>0,], #baseline covariates for left truncated subjects
                 Xtv_l=numeric(N_l), #set time-varying covariate to 0 ensuring no "jump" happened before study entry
                 Xtv_time_l=Y[Y[,3]>0,1], #set jump times to be the same as event times so that there is no excess computation
                 L_ind_o=rep(-1,N_o),
                 L_ind_m=rep(-1,N_m),
                 #intercept hyperparameters
                 prior_intercept = prior_intercept_num,
                 m_intercept = if(prior_intercept_num == 0) 1 else m_intercept,
                 sd_intercept = if(prior_intercept_num == 0) 1 else sd_intercept,
                 #scale hyperparameters
                 prior_scale = prior_scale_num,
                 a_scale=if(prior_scale_num == 0) 1 else a_scale,
                 b_scale=if(prior_scale_num == 0) 1 else b_scale,
                 #beta hyperparameters
                 prior_beta = prior_beta_num,
                 m_beta = if(prior_beta_num == 0) 1 else m_beta,
                 sd_beta = if(prior_beta_num == 0) 1 else sd_beta,
                 #beta_tv hyperparameters
                 prior_beta_tv = prior_beta_tv_num,
                 m_beta_tv = if(prior_beta_tv_num == 0) 1 else m_beta_tv,
                 sd_beta_tv=if(prior_beta_tv_num == 0) 1 else sd_beta_tv,
                 #baseline specification
                 weib_ind = as.numeric(tolower(baseline) %in% c("wb","weibull")),
                 #for use with tbp baseline
                 J=J,
                 dirichlet_alpha_fixed=as.numeric(dirichlet_alpha_fixed),
                 dirichlet_alpha_data=if(dirichlet_alpha_fixed) 
                   array(dirichlet_alpha_data,1) else numeric(0),
                 a_dirichlet_alpha=if(dirichlet_alpha_fixed) 1 else a_dirichlet_alpha,
                 b_dirichlet_alpha=if(dirichlet_alpha_fixed) 1 else b_dirichlet_alpha,
                 NULL
  )
  
  Xtv_time_mat <- as.matrix(cbind(
    pmin(Y[,1],Xtv_time),
    pmax(0,Y[,1] - pmin(Y[,1],Xtv_time))
  ))

  #Now, fix up left truncated times
  if(N_l > 0){
    #make a lookup table linking which observations are 
    #censored/uncensored with which are truncated/not
    #col 1 is the index of all the original data
    #col 2 is the index among the observed data
    #col 3 is the index among the censored data
    #col 4 is the corresponding index for the truncated data
    Y_ind_mat <- matrix(data = -1,nrow = nrow(Y),ncol=4)
    Y_ind_mat[,1] <- 1:nrow(Y)
    Y_ind_mat[Y[,1]==Y[,2],2] <- 1:sum(Y[,1]==Y[,2])
    Y_ind_mat[Y[,1]!=Y[,2],3] <- 1:sum(Y[,1]!=Y[,2])
    Y_ind_mat[Y[,3]>0,4] <- 1:sum(Y[,3]>0)
    d_list[["L_ind_o"]] <- Y_ind_mat[Y[,1]==Y[,2],4]
    d_list[["L_ind_m"]] <- Y_ind_mat[Y[,1]!=Y[,2],4]
    
  }
  
  #now, deal with time-varying effects
  if(!(tolower(tv_type) %in% c("pw","piecewise","spline"))){
    nP_tv <- 0
  } else{
    
    #now, set up knots
    if(is.null(knots)){
      stopifnot(!is.null(nP_tv))
      if(tolower(tv_type) %in% c("pw","piecewise")){
        knots <- c(0, quantile(Xtv_time_mat[Y[,1]==Y[,2] & Xtv_time_mat[,2] > 0,2], 
                               probs = seq(from=0,to=1,length.out = nP_tv+1))[-c(1,nP_tv+1)])
      } else if(tolower(tv_type) %in% c("spline")){
        #get knots according to the time after covariate jump
        knots<- get_default_knots_vec(y = Xtv_time_mat[Xtv_time_mat[,2] > 0,2],
                                      delta = as.numeric(Y[,1]==Y[,2])[Xtv_time_mat[,2] > 0], 
                                      p_tv = nP_tv, tv_type="rp")
      }
    }
    if(knots[1] != 0){ knots <- c(0,knots) }
    if(knots[length(knots)] < 1e8){ knots <- c(knots, Inf) }
    d_list[["knots_tv"]] <- knots
    
    #now, generate the bases here
    if(tolower(tv_type) %in% c("pw","piecewise")){
      basis <- pw_cum_mat(y=Xtv_time_mat[,2],
                                 knots=knots[-length(knots)], #remove the Inf knot
                                 intercept = TRUE) #note we include intercept here!
      d_list[["basis_piece_o"]] <- basis[Y[,1]==Y[,2],,drop=FALSE]
      d_list[["basis_piece_m"]] <- basis[Y[,1]!=Y[,2],,drop=FALSE]
      
      d_list[["basis_piece_l"]] <- matrix(data=0, nrow=N_l, ncol=NCOL(basis))

    } else if(tolower(tv_type) %in% c("spline")){
      basis <- get_basis_tv(x=Xtv_time_mat[,2], tv_type="rp",
                            knots=knots[-c(1,length(knots))], #remove the 0 and Inf knot
                            deriv = FALSE, intercept=TRUE, #note we include intercept now!
                            flexsurv_compatible = FALSE)
      dbasis <- get_basis_tv(x=Xtv_time_mat[,2], tv_type="rp",
                             knots=knots[-c(1,length(knots))], #remove the 0 and Inf knot
                             deriv = TRUE, intercept=TRUE, #note we include intercept now!
                             flexsurv_compatible = FALSE)
      d_list[["basis_spline_o"]] <- basis[Y[,1]==Y[,2],,drop=FALSE]
      d_list[["dbasis_spline_o"]] <- dbasis[Y[,1]==Y[,2],,drop=FALSE]
      d_list[["basis_spline_m"]] <- basis[Y[,1]!=Y[,2],,drop=FALSE]
      
      d_list[["basis_spline_l"]] <- matrix(data=0, nrow=N_l, ncol=NCOL(basis))
    }
    nP_tv <- NCOL(basis)
    d_list[["nP_tv_piece"]] <- d_list[["nP_tv_spline"]] <- nP_tv
    
  }
  
  
  par <- c("intercept","scale","beta")
  if(tolower(tv_type) %in% c("pw","piecewise")){
    par <- c(par,"beta_tv")
    if(tbp){
      par <- c(par,"w",if(!dirichlet_alpha_fixed) "dirichlet_alpha_param" else NULL)
      stan_fit <- rstan::sampling(stan_AFT_tbp_tvcov_piece_compiled,data=d_list,
                                  iter=n_warmup + n_sample,warmup=n_warmup,
                                  chains=n_chains,cores=n_cores, init = init, control = control,
                                  seed = seed, pars=c(par,"log_lik"), include=TRUE)
    } else{
      stan_fit <- rstan::sampling(stan_AFT_tvcov_piece_compiled,data=d_list,
                                  iter=n_warmup + n_sample,warmup=n_warmup,
                                  chains=n_chains,cores=n_cores, init = init, control = control,
                                  seed = seed, pars=c(par,"log_lik"), include=TRUE)
    }
  } else if(tolower(tv_type) %in% c("spline")){
    par <- c(par,"beta_tv")
    if(tbp){
      par <- c(par,"w",if(!dirichlet_alpha_fixed) "dirichlet_alpha_param" else NULL)
      stan_fit <- rstan::sampling(stan_AFT_tbp_tvcov_spline_compiled,data=d_list,
                                  iter=n_warmup + n_sample,warmup=n_warmup,
                                  chains=n_chains,cores=n_cores, init = init, control = control,
                                  seed = seed, pars=c(par,"log_lik"), include=TRUE)
    } else{
      stan_fit <- rstan::sampling(stan_AFT_tvcov_spline_compiled,data=d_list,
                                  iter=n_warmup + n_sample,warmup=n_warmup,
                                  chains=n_chains,cores=n_cores, init = init, control = control,
                                  seed = seed, pars=c(par,"log_lik"), include=TRUE)
    }
  } else{
    if(tbp){
      par <- c(par,"w",if(!dirichlet_alpha_fixed) "dirichlet_alpha_param" else NULL)
      stan_fit <- rstan::sampling(stan_AFT_tbp_invariant_compiled,data=d_list,
                                  iter=n_warmup + n_sample,warmup=n_warmup,
                                  chains=n_chains,cores=n_cores, init = init, control = control,
                                  seed = seed, pars=c(par,"log_lik"), include=TRUE)
    } else{
      stan_fit <- rstan::sampling(stan_AFT_invariant_compiled,data=d_list,
                                  iter=n_warmup + n_sample,warmup=n_warmup,
                                  chains=n_chains,cores=n_cores, init = init, control = control,
                                  seed = seed, pars=c(par,"log_lik"), include=TRUE)
    }
  }
  
  outlist <- list(
    stan_fit=stan_fit,
    knots=knots,
    pars=par,
    nP=nP,
    nP_tv=nP_tv,
    prior_intercept=prior_intercept,
    m_intercept=m_intercept,
    sd_intercept=sd_intercept,
    prior_scale=prior_scale,
    a_scale=a_scale,
    b_scale=b_scale,
    prior_beta=prior_beta,
    m_beta=m_beta,
    sd_beta=sd_beta,
    prior_beta_tv=prior_beta_tv,
    m_beta_tv=m_beta_tv,
    sd_beta_tv=sd_beta_tv,
    baseline=baseline,
    tv_type=tv_type,
    tbp=tbp, 
    J=if(tbp) J else NULL,
    dirichlet_alpha_fixed=dirichlet_alpha_fixed,
    a_dirichlet_alpha=a_dirichlet_alpha, 
    b_dirichlet_alpha=b_dirichlet_alpha,
    dirichlet_alpha_data=dirichlet_alpha_data,
    init = init,
    NULL
  )
}


V0_func <- function(t=NULL, t_base_time=NULL, t_tv_time=NULL,
              newXtv, newXtv_time, 
              tv_type, beta_tv,
              knots,basis){
  # browser()
  if(is.null(t_base_time) || is.null(t_tv_time)){
    t_base_time <- pmin(t,newXtv_time) 
    t_tv_time <- pmax(0,t-t_base_time)
  } else if(is.null(t)){
      t <- t_base_time + t_tv_time
  } else{ stop("must provide either t or t_base_time and t_tv_time")}
  
  
  if(is.null(basis) && tolower(tv_type) %in% c("pw","piecewise","spline") ){
    stopifnot(!is.null(knots))
    knots <- if(knots[1] !=0) c(0,knots) else knots
    knots <- if(!(tail(knots,1) > 1e8)) c(knots,Inf) else knots #choose 1e8 as "arbitrarily large"
    if(tolower(tv_type) %in% c("pw","piecewise")){
      basis <- pw_cum_mat(y=t_tv_time,
                                 knots=knots[-length(knots)], #remove the Inf knot
                                 intercept = TRUE) #note we include intercept here!
    } else if (tv_type == "spline"){
      #if tv_type is spline but there is not time varying time, then skip basis formation.
      if(all(t_tv_time==0)){
        tv_type <- "none"
      } else{
        basis <- get_basis_tv(x=t_tv_time,
                              tv_type="rp",
                              knots=knots[-c(1,length(knots))],
                              deriv = FALSE, intercept=TRUE,
                              flexsurv_compatible = FALSE)
      }
    }
  }
  
  if(tolower(tv_type) %in% c("pw","piecewise")){
    if(all(newXtv==0)){
      V0_temp <- t
    } else if(all(newXtv==1)){
      V0_temp <- t_base_time + basis %*% exp(-beta_tv)
    } else{
      #if the new value isn't 0 or 1, then newXtv must be the same length as t_tv_time so that the
      #dimensionality of newXtv %*% t(beta_tv) matches dimensionality of basis matrix.
      stopifnot(length(newXtv)==NROW(t_tv_time))
      V0_temp <- t_base_time + rowSums(basis * exp(-newXtv %*% t(beta_tv)))
    }
  } else if (tolower(tv_type) == "spline"){
    V0_temp <- t_base_time + 
      t_tv_time * exp(-newXtv * basis %*% beta_tv)
  } else{
    V0_temp <- t
  }
  V0_temp
}

S_func <- function(t,
                   newXtv, newXtv_time, 
                   tv_type, beta_tv,
                   knots, basis,
                   oldXmat, newXmat, beta,
                   S0_func, baseline, intercept, scale,
                   tbp, w,
                   type){
  # browser()
  V0_temp <- V0_func(t = t,
                     newXtv=newXtv, newXtv_time=newXtv_time, 
                     tv_type=tv_type, beta_tv=beta_tv,
                     knots=knots,basis=basis)
  
  if(is.null(S0_func)){
    stopifnot(!(is.null(baseline) || is.null(intercept) || is.null(scale)|| is.null(tbp)))
    S0_func <- function(q,inter,sc,tbp,w){
      if(tolower(baseline=="weibull")){
        temp <- pweibull(q=q,scale=exp(inter),shape=sc,lower.tail = FALSE,log.p = FALSE)
      } else{
        temp <- plnorm(q=q,meanlog=inter,sdlog=sc,lower.tail = FALSE,log.p=FALSE)
      }
      if(tbp){
        J <- length(w)
        temp <- splines2::bernsteinPoly(x = temp,
                                          degree = J-1,
                                          intercept = TRUE,
                                          Boundary.knots = c(0,1),
                                          integral = TRUE) %*% w * J
      }
      temp
    }
  }
  
  if(type=="marginal"){
    stopifnot(!is.null(oldXmat))
    # oldXmatbeta <- oldXmat %*% beta
    if(NROW(V0_temp)==1){
      S_out <- mean(sapply(X = 1:NROW(oldXmat),
                           FUN = function(x){S0_func(q=V0_temp,
                                                     inter = intercept + oldXmat[x,] %*% beta,
                                                     sc = scale, tbp = tbp, w = w)}))
    } else{
      S_out <- rowMeans(sapply(X = 1:NROW(oldXmat),
                               FUN = function(x){S0_func(q=V0_temp,
                                                         inter = intercept + oldXmat[x,] %*% beta,
                                                         sc = scale, tbp = tbp, w = w)}))
    }
  } else {
    S_out <- S0_func(q=V0_temp, inter = intercept + newXmat %*% beta, 
                        sc = scale, tbp = tbp, w = w)
  }
  S_out
}

#not actually used for now.
Vinv_func <- function(t,
                      newXtv, newXtv_time, 
                      tv_type, beta_tv,
                      knots, basis,
                      newXmat, beta){
  
  if(tolower(tv_type)%in% c("pw","piece","piecewise")){
    ##something here
  } else if(tolower(tv_type)=="spline"){
    tryCatch(stats::uniroot(f = function (t){
      base_time <- pmin(t,newXtv_time)
      tv_time <- pmax(0,t-base_time)
      # print(paste("tv_time:",tv_time))
      if(tv_time==0){
        basis_spline <- t(rep(0,length(beta_tv)))
      } else{
        basis_spline <- get_basis_tv(x=tv_time,knots=knots,tv_type = "rp",
                                     intercept = TRUE,deriv = FALSE,flexsurv_compatible = FALSE)
      }
      t_obj - exp(- x_base %*% beta_base) * (base_time + tv_time * exp(-basis_spline %*% beta_tv * x_tv))},
      lower = lower, upper = upper)$root,
      error=function(e){return(NA)})
  } else{
    
  }
}

Sinv_func <- function(p,
                   newXtv, newXtv_time, 
                   tv_type, beta_tv,
                   knots,basis,
                   oldXmat, newXmat, beta,
                   S0_func, baseline, intercept, scale, 
                   tbp, w,
                   type, lower=0, upper=1000){
  
  # if type is marginal, or tv_type is spline, go the numerical route, 
  # otherwise maybe I can manually invert!
  # if(tolower(type)=="marginal"){
    # Mark Clements has a truly vectorized version of uniroot in the rstpm2 package
    # So, I can also try using it...
    #https://stat.ethz.ch/pipermail/r-help/2019-April/462477.html
  tryCatch(rstpm2::vuniroot(f = function (t){p -
      S_func(t=t,newXtv=newXtv, newXtv_time=newXtv_time, 
             tv_type=tv_type, beta_tv=beta_tv,
             knots=knots,basis=basis,
             oldXmat=oldXmat,newXmat=newXmat,beta=beta,
             S0_func=S0_func,baseline=baseline,intercept=intercept,scale=scale,
             tbp = tbp, w = w,
             type=type)},
      lower = rep(lower,length(p)), upper = rep(upper,length(p)))$root,
  error=function(e){return(rep(NA,length(t)))})
  
    # #For now, this is just the uniroot function put into sapply to "vectorize" it
    # #https://stat.ethz.ch/pipermail/r-help/2011-April/276457.html
    # sapply(p,function(x){tryCatch(stats::uniroot(f = function (t){x -
    #     S_func(t=t,newXtv=newXtv, newXtv_time=newXtv_time,
    #            tv_type=tv_type, beta_tv=beta_tv,
    #            knots=knots,basis=basis,
    #            oldXmat=oldXmat,newXmat=newXmat,beta=beta,
    #            S0_func=S0_func,baseline=baseline,intercept=intercept,scale=scale,type=type)},
    #     lower = lower, upper = upper)$root,
    #     error=function(e){return(NA)})})
  # } else{
  #    # to perhaps speed things up further, for conditional piecewise and time-invariant setting,
  #    # could directly implement Vinv and compute things that way.
  # }

}

#function to generate output matrix of conditional or marginal survivor functions
#default is to show 1 vs 0 of the time-varying covariate
predict.AFTtvstan <- function(stan_fit, 
                              t_seq,
                              oldXmat=NULL, 
                              newXmat=matrix(data=0,nrow=1,ncol=stan_fit$nP), 
                              newXtv=0, newXtv_time=0,
                              type=NULL, chain = NULL,
                              thin=1, verbose=FALSE,
                              summarize=TRUE){
  # browser()
  # n_ind <- length(newXtv)
  # if(length(newXtv_time)==1){newXtv_time <- rep(newXtv_time,n_ind)}
  # if(is.null(newXmat)){
  #   newXmat <- matrix(data=0,nrow=2,ncol=stan_fit$nP)
  # } else if(is.vector(newXmat)){
  #   newXmat <- matrix(data=newXmat,nrow=n_ind,ncol=nP,byrow=TRUE)
  # }
  # stopifnot(nrow(newXmat)==length(newXtv) && length(newXtv) == length(newXtv_time))
  
  stan_array <- as.array(stan_fit$stan_fit)
  
  if(!is.null(chain)){
    stan_array <- stan_array[,paste0("chain:",chain),,drop=FALSE]
  }
  
  
  int_temp <- as.vector(stan_array[,,"intercept"])
  scale_temp <- as.vector(stan_array[,,"scale"])
  if(stan_fit$nP > 1){
    beta_temp <- apply(X = stan_array[,,paste0("beta[",1:stan_fit$nP,"]"),drop=FALSE],
                       MARGIN=3,FUN = as.vector)
  } else if(stan_fit$nP == 1){
    beta_temp <- as.matrix(as.vector(stan_array[,,"beta[1]"]))
  } else{
    newXmat <- matrix(data=0,nrow=1,ncol=1)
    beta_temp <- as.matrix(rep(0,length(int_temp)))
  }
  if(stan_fit$tv_type != "none"){
    if(stan_fit$nP_tv > 1){
      beta_tv_temp <- apply(X = stan_array[,,paste0("beta_tv[",1:stan_fit$nP_tv,"]"),drop=FALSE],
                            MARGIN=3,FUN = as.vector)
    } else{
      beta_tv_temp <- as.matrix(as.vector(stan_array[,,"beta_tv[1]"]))
    }
  }

  if(stan_fit$tbp){
    w_temp <- apply(X = stan_array[,,paste0("w[",1:stan_fit$J,"]"),drop=FALSE],
                    MARGIN=3,FUN = as.vector)
  }
    
    
  # S0_func <- if(tolower(stan_fit$baseline=="weibull")){
  #     function(q,inter,sc){pweibull(q=q,scale=exp(inter),shape=sc,lower.tail = FALSE,log.p = FALSE)}
  #   } else{
  #     function(q,inter,sc){plnorm(q=q,meanlog=inter,sdlog=sc,lower.tail = FALSE,log.p=FALSE)}
  #   }
  
  tseq_base_time <- pmin(t_seq,newXtv_time) 
  tseq_tv_time <- pmax(0,t_seq-tseq_base_time)
  
  if(tolower(stan_fit$tv_type) %in% c("pw","piecewise")){
    tseq_tv_basis <- pw_cum_mat(y=tseq_tv_time,
                                       knots=stan_fit$knots[-length(stan_fit$knots)], #remove the Inf knot
                                       intercept = TRUE) #note we include intercept here!
  } else if (stan_fit$tv_type == "spline"){
    tseq_tv_basis <- get_basis_tv(x=tseq_tv_time,
                                  tv_type="rp",
                                  knots=stan_fit$knots[-c(1,length(stan_fit$knots))],
                                  deriv = FALSE, intercept=TRUE,
                                  flexsurv_compatible = FALSE)
  } else{
    tseq_tv_basis <- NULL
  }
  
  samp_ind <- (1:length(int_temp))[1:length(int_temp) %% thin == 0]
  
  out_mat <- matrix(data=NA,nrow=length(samp_ind),ncol=length(t_seq))
  # out_array <- array(data = NA,dim = c(length(int_temp),length(t_seq),n_ind))
  
  for(ind in 1:length(samp_ind)){
    if(verbose && ind %% 10 == 0){print(paste0(ind," of ",length(samp_ind)," after thinning by factor of ", thin))}
    j <- samp_ind[ind]
    
    out_mat[ind,] <- S_func(t=t_seq, 
                            newXtv = newXtv,newXtv_time = newXtv_time, 
                            oldXmat = oldXmat, newXmat = newXmat,
                            beta_tv = beta_tv_temp[j,], beta = beta_temp[j,],
                            basis = tseq_tv_basis, knots = NULL,
                            S0_func = NULL, baseline = stan_fit$baseline, 
                            intercept = int_temp[j],scale = scale_temp[j],
                            tbp = stan_fit$tbp, w = if(stan_fit$tbp) w_temp[j,] else NULL,
                            tv_type = stan_fit$tv_type, type=type)
  }
  if(summarize){
    return(t(apply(X = out_mat,MARGIN = 2,FUN = quantile,probs=c(0.5,0.025,0.975))))
  } else{
    return(out_mat)
  }
}




AF.AFTtvstan_tvcov <- function(stan_fit, p_seq, t_Xtv_seq=0, p_Xtv_seq=NULL, 
                               oldXmat=NULL, newXtv_num=1, newXtv_denom=0, 
                         newXmat_num=matrix(data=0,nrow=1,ncol=stan_fit$nP),
                         newXmat_denom=matrix(data=0,nrow=1,ncol=stan_fit$nP), 
                         type=NULL, thin=1, verbose=FALSE,long_out=FALSE,
                         summarize=TRUE){
  # browser()
  stan_array <- as.array(stan_fit$stan_fit)
  int_temp <- as.vector(stan_array[,,"intercept"])
  scale_temp <- as.vector(stan_array[,,"scale"])
  if(stan_fit$nP > 1){
    beta_temp <- apply(X = stan_array[,,paste0("beta[",1:stan_fit$nP,"]")],
                       MARGIN=3,FUN = as.vector)
  } else if(stan_fit$nP == 1){
    beta_temp <- as.matrix(as.vector(stan_array[,,"beta[1]"]))
  } else{
    newXmat_num <- newXmat_denom <- matrix(data=0,nrow=1,ncol=1)
    beta_temp <- as.matrix(rep(0,length(int_temp)))
  }
  if(stan_fit$tv_type != "none"){
    if(stan_fit$nP_tv > 1){
      beta_tv_temp <- apply(X = stan_array[,,paste0("beta_tv[",1:stan_fit$nP_tv,"]")],
                            MARGIN=3,FUN = as.vector)
    } else{
      beta_tv_temp <- as.matrix(as.vector(stan_array[,,"beta_tv[1]"]))
    }
  }
  
  if(stan_fit$tbp){
    w_temp <- apply(X = stan_array[,,paste0("w[",1:stan_fit$J,"]")],
                    MARGIN=3,FUN = as.vector)
  }
  
  samp_ind <- (1:length(int_temp))[1:length(int_temp) %% thin == 0]
  n_samp <- length(samp_ind)
  
  #if the denominator is null, we can substantially speed things up! so flag that
  Xtv_denom_nullflag <- all(newXtv_denom == 0)
  
  #if p sequence is given, override t sequence given, otherwise use t sequence
  p_Xtv_givenflag <- !is.null(p_Xtv_seq)
  t_Xtv_nonstandard <- p_Xtv_givenflag | any(t_Xtv_seq != 0) #are we dealing with a time-varying covariate or no?
  Xtv_seq_length <- if(p_Xtv_givenflag) length(p_Xtv_seq) else length(t_Xtv_seq)
  #if t is given, then we will estimate and store the corresponding quantiles p
  #if p is given, then we will estimate and store the corresponding survival times t
  out_Xtv_component <- array(data = NA,dim = c(Xtv_seq_length,n_samp))
  out_array <- array(data = 1,dim = c(length(p_seq),Xtv_seq_length,n_samp))
  

  for(ind in 1:length(samp_ind)){
    if(verbose && ind %% 10 == 0){print(paste0(ind," of ",n_samp," after thinning by factor of ", thin))}
    j <- samp_ind[ind]
    
    #if p quantiles are specified, fill in t_Xtv_seq, and vice versa
    if(p_Xtv_givenflag){
      #start with an empty vector of zeroes
      t_Xtv_seq <- numeric(length(p_Xtv_seq))
      #replace nonzero elements, but if p = 1 then leave Sinv(p) = 0
      t_Xtv_seq[!(p_Xtv_seq %in% c(0,1))] <- Sinv_func(p=p_Xtv_seq[!(p_Xtv_seq %in% c(0,1))],
                                                          newXtv = newXtv_denom, newXmat = newXmat_denom, #note denom used
                                                          oldXmat = oldXmat, newXtv_time = 0,
                                                          beta = beta_temp[j,], beta_tv = beta_tv_temp[j,],
                                                          tv_type = stan_fit$tv_type, knots = stan_fit$knots, basis=NULL,
                                                          S0_func = NULL, baseline = stan_fit$baseline,
                                                          intercept = int_temp[j], scale = scale_temp[j],
                                                          tbp = stan_fit$tbp, w = if(stan_fit$tbp) w_temp[j,] else NULL,
                                                          type = type)
      t_Xtv_seq[p_Xtv_seq == 0] <- Inf
      out_Xtv_component[,ind] <- t_Xtv_seq
    } else{ #if t times are specified, compute corresponding p's
      # p_Xtv_seq <- S_func(p=t_Xtv_seq,
      #                         newXtv = newXtv_denom, newXmat = newXmat_denom, #note denom used
      #                         oldXmat = oldXmat, newXtv_time = 0,
      #                         beta = beta_temp[j,], beta_tv = beta_tv_temp[j,],
      #                         tv_type = stan_fit$tv_type, knots = stan_fit$knots, basis=NULL,
      #                         S0_func = NULL, baseline = stan_fit$baseline,
      #                         intercept = int_temp[j], scale = scale_temp[j],
      #                         tbp = stan_fit$tbp, w = if(stan_fit$tbp) w_temp[j,] else NULL,
      #                         type = type)
      # out_Xtv_component[,ind] <- p_Xtv_seq
    }
    
    #if denominator does not have time-varying component, then we only need to compute it once instead of for every (possibly) changing value of Xtv 
    if(Xtv_denom_nullflag){
      S_inv_denom_temp <- Sinv_func(p=p_seq, 
                                    newXtv = newXtv_denom, newXtv_time = newXtv_denom, #defaults to 0
                                    oldXmat = oldXmat, newXmat = newXmat_denom, 
                                    beta = beta_temp[j,], beta_tv = beta_tv_temp[j,], 
                                    tv_type = stan_fit$tv_type, knots = stan_fit$knots, basis=NULL,
                                    S0_func = NULL, baseline = stan_fit$baseline,
                                    intercept = int_temp[j], scale = scale_temp[j],
                                    tbp = stan_fit$tbp, w = if(stan_fit$tbp) w_temp[j,] else NULL,
                                    type = type)
    }

    for(k in 1:length(t_Xtv_seq)){
      #if denominator DOES have time-varying component that's not 0, then recompute it for entire t_Xtv_seq
      if(!Xtv_denom_nullflag){
        S_inv_denom_temp <- Sinv_func(p=p_seq, 
                                      newXtv = newXtv_denom, newXtv_time = t_Xtv_seq[k],
                                      oldXmat = oldXmat, newXmat = newXmat_denom, 
                                      beta = beta_temp[j,], beta_tv = beta_tv_temp[j,], 
                                      tv_type = stan_fit$tv_type, knots = stan_fit$knots, basis=NULL,
                                      S0_func = NULL, baseline = stan_fit$baseline,
                                      intercept = int_temp[j], scale = scale_temp[j],
                                      tbp = stan_fit$tbp, w = if(stan_fit$tbp) w_temp[j,] else NULL,
                                      type = type)
      }
      #only need to explicitly compute for areas of the surface where the 
      #covariate jump quantile is before the survival quantile of interest
      #set everything else to 1.
      if(t_Xtv_nonstandard){
        relevant_AFs <- which(S_inv_denom_temp >= t_Xtv_seq[k])
        if(length(relevant_AFs) == 0){next}
        #change only the relevant values
        S_inv_num_temp <- Sinv_func(p=p_seq[relevant_AFs], newXtv = newXtv_num, newXtv_time = t_Xtv_seq[k],
                                    oldXmat = oldXmat, newXmat = newXmat_num, 
                                    beta = beta_temp[j,], beta_tv = beta_tv_temp[j,], 
                                    tv_type = stan_fit$tv_type, knots = stan_fit$knots, basis=NULL,
                                    S0_func = NULL, baseline = stan_fit$baseline,
                                    intercept = int_temp[j], scale = scale_temp[j],
                                    tbp = stan_fit$tbp, w = if(stan_fit$tbp) w_temp[j,] else NULL,
                                    type = type)
        out_array[relevant_AFs,k,ind] <- S_inv_num_temp/S_inv_denom_temp[relevant_AFs]
      } else{
        S_inv_num_temp <- Sinv_func(p=p_seq, newXtv = newXtv_num, newXtv_time = t_Xtv_seq[k],
                                    oldXmat = oldXmat, newXmat = newXmat_num, 
                                    beta = beta_temp[j,], beta_tv = beta_tv_temp[j,], 
                                    tv_type = stan_fit$tv_type, knots = stan_fit$knots, basis=NULL,
                                    S0_func = NULL, baseline = stan_fit$baseline,
                                    intercept = int_temp[j], scale = scale_temp[j],
                                    tbp = stan_fit$tbp, w = if(stan_fit$tbp) w_temp[j,] else NULL,
                                    type = type)
        out_array[,k,ind] <- S_inv_num_temp/S_inv_denom_temp
      }
    }
  }
  
  if(summarize){
    if(length(t_Xtv_seq)==1){
      out_mat <- t(apply(out_array,MARGIN = c(1,2),FUN = quantile, probs=c(0.5,0.025,0.975),na.rm=TRUE)[,,1])
      # t(apply(out_mat,MARGIN = c(1,2),FUN = quantile, probs=c(0.5,0.025,0.975), na.rm=TRUE))
    } else{
      out_mat <- apply(out_array,MARGIN = c(1,2),FUN = median,na.rm=TRUE)
    }
    out_Xtv_component_median <- apply(out_Xtv_component,MARGIN = 1,FUN = median,na.rm=TRUE)
    
    if(long_out){
      # if(p_Xtv_givenflag){
      #   as.data.frame(
      #     cbind(expand.grid(p=p_seq,Sinvxp=p_Xtv_seq),
      #           Sinvx=rep(out_Xtv_component_median,each=length(p_seq)),
      #           AF = round(as.vector(out_mat),4)))
      # } else{
      #   as.data.frame(
      #     cbind(expand.grid(p=p_seq,Sinvx=t_Xtv_seq),
      #           Sinvxp=rep(out_Xtv_component_median,each=length(p_seq)),
      #           AF = round(as.vector(out_mat),4)))
      # }
      
      #this assumes we've supplied a t_Xtv_seq for now.
      return(as.data.frame(cbind(expand.grid(p=p_seq,Sinvx=t_Xtv_seq),AF = round(as.vector(out_mat),4))))
    } else{
      return(out_mat)
    }
  } else{
    return(out_array)
  }

}

