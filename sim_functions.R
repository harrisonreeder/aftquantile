
#' Generate Matrix of Values of Inverse of Integrated Piecewise Constant Function
#'
#' @param y
#' @param knots
#' @param x_base
#' @param beta_base
#' @param x_tv
#' @param beta_tv
#'
#' @return
#' @export
pw_cum_mat_inv <- function(y, knots, x_base=x_base,beta_base=beta_base, x_tv, beta_tv){
  
  y_adj <- y #as.numeric(y * exp(-x_base %*% beta_base))
  #p0 is number of knots (including 0)
  
  #these are the knots on the original time scale
  if(knots[1] != 0){knots <- c(0,knots)}
  
  #n by p0 matrix with knots for each subject, shifted by baseline amount
  adj_knots_mat <- exp(-x_base %*% beta_base) %*% t(knots)
  #n by (p0-1) matrix with differences for each subject, shifted by baseline amount
  # tstar1*exp(-x_base trans beta_base), (tstar2 - tstar1)*exp(-x_base trans beta_base), etc.
  adj_knots_diff <- t(diff(t(adj_knots_mat)))
  
  #p0
  nP0 <- length(beta_tv)
  n <- length(y)
  
  stopifnot(nP0 > 1)
  
  #n by (p0-1) matrix with (exp(-x_i*beta1),exp(-x_i*beta2),...), omitting last beta
  temp_out <- exp(-x_tv %*% t(beta_tv[-nP0]))
  
  
  #n by (p0-1) matrix with each column giving corrected interval length between each knot
  #tstar1*exp(-x_base trans beta_base) (tstar2-tstar1)*exp(-x_base trans beta_base - x*beta1)
  adj_knots_diff_mat <- adj_knots_diff
  adj_knots_diff_mat[,-1] <- adj_knots_diff_mat[,-1] * temp_out
  
  #n by p0 matrix with each column being the jth knot for the ith subject
  adj_knots_mat <- cbind(0,t(apply(adj_knots_diff_mat,MARGIN = 1,cumsum)))
  
  #add column of 0's as placeholder, as this will be output matrix
  adj_knots_diff_mat <- cbind(adj_knots_diff_mat,0)
  
  #an n length vector of integers indicating which interval the observation is in
  #https://stackoverflow.com/questions/27472948/vectorization-of-findinterval
  cut_cats <- rowSums(y_adj >= adj_knots_mat) #mild question of whether this should be open or closed
  
  #an n length vector capturing the residual time accrued in the final interval of each observation
  #https://stackoverflow.com/questions/20617371/getting-elements-of-a-matrix-with-vectors-of-coordinates
  y_last <- y_adj - adj_knots_mat[cbind(1:n,cut_cats)]
  
  #a loop through each observation to finalize the matrix of time intervals
  for(i in 1:n){
    adj_knots_diff_mat[i,cut_cats[i]] <- y_last[i]
    adj_knots_diff_mat[i,-c(1:cut_cats[i])] <- 0
  }
  
  #adj_knots_diff_mat tallies up the total accrued time of individual i through all intervals
  #(after transforming to the inverse scale)
  
  return(adj_knots_diff_mat)
}




#' Compute Inverse of Integrated Covariate Process V
#'
#' @param t_obj
#' @param x_tv
#' @param beta_tv
#' @param xbeta_base
#' @param x_base
#' @param beta_base
#' @param tv_type
#' @param inv_basis
#' @param knots
#'
#' @return
#' @export
Vx_inv <- function(t_obj, x_tv=NULL, beta_tv=NULL,
                   xbeta_base=NULL, x_base=NULL, beta_base=NULL,
                   tv_type, inv_basis=NULL, knots=NULL){
  # browser()
  temp <- NULL
  
  if(is.null(xbeta_base)){
    if(is.null(x_base)){
      if(is.null(beta_base)){
        warning("we're assuming no baseline covariates, because xbeta_base, x_base, and beta_base are all NULL.")
        xbeta_base <- 0
        x_base <- numeric(length(t_obj))
        beta_base <- 0
      } else{
        warning("x_base and xbeta_base are both NULL, so we assume all baseline covariates are set to 1.")
        #just assume that we want to plot with x set to 1 for everyone
        xbeta_base <- sum(beta_base)
        x_base <- matrix(1, nrow=length(t_obj), ncol=length(beta_base))
      }
    } else{
      xbeta_base <- x_base %*% as.matrix(beta_base)
    }
  }
  
  if(tv_type=="piecewise"){
    if(is.null(x_tv)){
      warning("x_tv was null, so we're replacing it with 1's")
      x_tv <- rep(1,length(t_obj))
    }
    if(is.null(inv_basis)){
      inv_basis <- pw_cum_mat_inv(y=t_obj,knots=knots,
                                  x_base=x_base,beta_base=beta_base,
                                  x_tv=x_tv,beta_tv=beta_tv)
    }
    #matrix with columns
    #exp(x_base trans beta_base)  exp(x_base trans beta_base + x_tv * beta_tv1) exp(x_base trans beta_base + x_tv * beta_tv2)
    mult_temp <- cbind(1,exp(x_tv %*% t(beta_tv))) * as.numeric(exp(xbeta_base))
    
    #vector with row sums of the above quantities
    temp <- as.numeric(rowSums(inv_basis * mult_temp))
  } else{
    
    #now that we're in non-piecewise
    tvec_xbeta <- exp(xbeta_base) * t_obj
    
    if(is.null(x_tv)){
      if(is.null(beta_tv)){
        if(tv_type != "baseline"){
          stop("for everything except tv_type='baseline', must provide either xbeta_tv, or beta_tv (and optionally x_tv).")
        }
      } else{
        xbeta_tv <- rep(1,length(t_obj)) %*% t(beta_tv)
      }
    } else{
      xbeta_tv <- x_tv %*% t(beta_tv)
    }
    
    
    if(tv_type=="baseline"){temp <- tvec_xbeta}
    else if(tv_type=="constant"){temp <- tvec_xbeta*exp(xbeta_tv)}
    else if(tv_type=="step"){temp <- tvec_xbeta + (exp(xbeta_tv)-1)*(tvec_xbeta-tstar)*as.numeric(tvec_xbeta>tstar)}
    else if(tv_type=="logplusone"){
      #faster ifelse: https://github.com/ICJIA/r-user-group/issues/11
      flag <- xbeta_tv != 1
      
      temp <- rep(NA, length(flag))
      temp[flag] <- ( (1-xbeta_tv)*tvec_xbeta + 1 )^(1/(1-xbeta_tv)) - 1
      temp[!flag] <- exp(tvec_xbeta) - 1
    } else{stop("tv_type must be 'piecewise', 'baseline', 'constant', 'step', 'logplusone'")}
  }
  as.numeric(temp)
}



#' Simulate Data under Time-Varying AFT Model
#'
#' @param x_base
#' @param x_tv
#' @param beta_base_true
#' @param beta_tv_true
#' @param tv_type
#' @param S0_inv_func
#' @param dist
#' @param intercept
#' @param scale
#' @param cens
#' @param knots
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sim_AFTtv_better <- function (x_base, x_tv, beta_base_true, beta_tv_true,tv_type,
                              S0_inv_func=NULL,dist=NULL,intercept=NULL,scale=NULL,
                              cens, knots,...) {
  n <- dim(x_base)[1]
  p <- dim(x_base)[2]
  
  #this is just one way of sampling from the baseline distribution,
  #though you could also sample from it "directly" !!
  if(!is.null(S0_inv_func)){
    u <- stats::runif(n = n, min = 0, max = 1)
    T_baseline <- S0_inv_func(u, ...)
  } else if(is.null(dist) ){
    stop("must provide either S0_inv_func or specify dist as 'weibull' or 'lognormal'")
  } else{
    stopifnot(tolower(dist) %in% c("weibull","lognormal") & !is.null(intercept) & !is.null(scale))
    T_baseline <- switch(tolower(dist),
                         weibull=rweibull(n = n,scale=exp(intercept),shape=1/scale),
                         lognormal=rlnorm(n = n,meanlog = intercept,sdlog=scale))
  }
  
  T_temp <- Vx_inv(t_obj = T_baseline, x_base=x_base, x_tv=x_tv,
                   beta_base=beta_base_true, beta_tv=beta_tv_true,
                   tv_type=tv_type, knots=knots)
  
  delta <- rep(NA, n)
  y <- T_temp
  C_temp <- stats::runif(n = n, min = cens[1], max = cens[2])
  
  ind1 <- which(T_temp < C_temp)
  y[ind1] <- T_temp[ind1]
  delta[ind1] <- 1
  ind0 <- which(T_temp >= C_temp)
  y[ind0] <- C_temp[ind0]
  delta[ind0] <- 0
  
  data.frame(y=y, delta=delta)
}