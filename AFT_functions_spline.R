####HELPER KNOT AND BASIS FUNCTIONS####

#' Generate Matrix of Values of Integrated Piecewise Constant Function
#'
#' This helper function takes in a vector of event times, and
#'   generates a matrix having each row showing the time accrued within each
#'   interval between consecutive elements of a vector of knots.
#'
#' @param y vector of event times.
#' @param knots increasing vector of cutpoints. If it does not start with 0, one will be appended to the start.
#'   However, it should not include Inf at the end.
#' @param intercept if true, includes column corresponding with 'first' interval, and if false does not.
#'   It makes sense to include an intercept if the time-varying covariate is not also included in the "baseline",
#'   otherwise, there would be an identifiability issue.
#'
#' @return a numeric matrix with, with (number of nonzero knots + 1) columns, and with rows corresponding to elements of y.
#' @export
pw_cum_mat <- function(y, knots, intercept=TRUE){
  if(knots[1] != 0){knots <- c(0,knots)}
  #vector giving the length of the intervals between each knot
  knots_diff <- diff(knots)
  #count of the number of intervals in each list
  num_int <- c(length(knots))
  n <- length(y)
  #matrix with each row being an observation, and each column being an interval, capturing how much time is accrued in each interval
  knots_mat <- matrix(c(knots_diff,0),nrow=n,ncol=num_int,byrow = TRUE)
  #an n length vector of integers indicating which interval the observation is in
  cut_cats <- findInterval(x = y, vec = knots)
  #an n length vector capturing the residual time accrued in the final interval of each observation
  y_last <- y-knots[cut_cats]
  #a loop through each observation to finalize the matrix of time intervals
  for(i in 1:n){
    knots_mat[i,cut_cats[i]] <- y_last[i]
    knots_mat[i,-c(1:cut_cats[i])] <- 0
  }
  
  #removing the intercept changes the parameterization, so that every parameter is a change from the "reference"
  #group which is the first interval. (it basically subtracts the parameter from the first interval from every subsequent interval
  #could be worth changing to instead be the 'chained' parameterization,
  #in which we have each one be difference from the one previous?
  if(!intercept){
    if(ncol(knots_mat) <= 2){stop("include at least two nonzero knots.")}
    knots_mat <- knots_mat[,-1]
  }
  return(knots_mat)
}




ns_d <- function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x))
{
  
  #inspired by rstpm2 package nsxD function, extending the ns function from splines
  #generate derivative basis of natural cubic spline function
  #I purposefully write this function as the derivative of ns, for use with Royston parmar model
  
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  
  #temporarilty set aside all NA values
  if (nas <- any(nax))
    x <- x[!nax]
  
  #create a boolean flag for whether any of the given x values are outside of the boundaries
  outside <- if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    (ol <- x < Boundary.knots[1L]) | (or <- x > Boundary.knots[2L])
  } else {
    if (length(x) == 1L)
      Boundary.knots <- x * c(7, 9)/8
    FALSE
  }
  
  #if df is given instead of actual knot locations, fill in knots based on quantiles
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - 1L - intercept
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d",
                       1L + intercept), domain = NA)
    }
    #notice this funny nesting, 'knots' inside is sequence of percentiles, and then reassigned as the quantiles of internal data
    knots <- if (nIknots > 0L) {
      knots <- seq.int(from = 0, to = 1, length.out = nIknots +
                         2L)[-c(1L, nIknots + 2L)]
      stats::quantile(x[!outside], knots)
    }
  } else{ nIknots <- length(knots) }
  Aknots <- sort(c(rep(Boundary.knots, 4L), knots)) #this code essentially puts two copies of each boundary on the ends, forming intervals
  if (any(outside)) {
    basis <- array(0, c(length(x), nIknots + 4L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      tt <- splines::splineDesign(Aknots, rep(k.pivot, sum(ol)), 4,
                                  1)
      basis[ol, ] <- tt
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      tt <- splines::splineDesign(Aknots, rep(k.pivot, sum(or)), 4,
                                  1)
      basis[or, ] <- tt
    }
    if (any(inside <- !outside))
      basis[inside, ] <- splines::splineDesign(Aknots, x[inside],
                                               4,derivs = 1)
  } else{ basis <- splines::splineDesign(Aknots, x, ord = 4L,derivs = 1) }
  const <- splines::splineDesign(Aknots, Boundary.knots, ord = 4L, derivs = c(2L,
                                                                              2L))
  if (!intercept) {
    const <- const[, -1, drop = FALSE]
    basis <- basis[, -1, drop = FALSE]
  }
  qr.const <- qr(t(const))
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L),
                                                     drop = FALSE])#uses QR decomp to orthogonalize the basis
  n.col <- ncol(basis)
  
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = 3L, knots = if (is.null(knots)) numeric() else knots,
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("ns_d", "basis", "matrix")
  basis
}


#' Generate List of Knot Location Vectors from Event Times
#'
#' This function creates a list containing three numeric vectors. Each numeric vector
#'   is a sequence of increasing integers giving the location of knots used for
#'   spline and piecewise baseline hazard specifications. This function
#'   generates this list according to the conventions and recommended locations
#'   of these knots, which depends on the choice of hazard specification, number
#'   of knots requested, and distribution of observed event times.
#'
#' @param y vector of event times (for now, set aside interval censoring)
#' @param delta vector of event indicators (for now, set aside interval censoring)
#' @param p_tv Integer indicating how many parameters should be associated with the time-varying effect (total, including intercept)
#'
#' @return A list of three increasing sequences of integers, each corresponding to
#'   the knots for the flexible model on the corresponding transition baseline hazard.
#' @export
get_default_knots_vec <- function(y,delta,p_tv,tv_type){

  if(tolower(tv_type) %in% c("bspline","bs")){
    # Recall, cubic bsplines are built by columns that are basically an orthogonalized transformation of:
    # (1, y, y^2, y^3, (y-internalknot1)_+^3,(y-internalknot2)_+^3,(y-internalknot3)_+^3),...
    # Therefore, for cubic spline the minimum number of parameters in a hazard must be 4 (aka no internal knots)
    if(p_tv < 4){stop("Cubic B-Spline Specification must have at least 4 parameters in each hazard.")}
    #vector of quantiles at which we will set knots
    quantile_seq <- seq(from = 0,to = 1, length.out = p_tv-2)[-c(1,p_tv-2)]
    #vector of knots corresponding to chosen quantiles (including boundary knots of 0 and maximum observed time.)
    knots_tv <- c(0,stats::quantile(y[delta==1],quantile_seq),max(y))
  } else if(tolower(tv_type) %in% c("piecewise","pw")){
    #vector of quantiles at which we will set knots
    quantile_seq <- seq(from = 0,to = 1, length.out = p_tv+1)[-c(1,p_tv+1)]
    #vector of knots corresponding to chosen quantiles (including lower boundary knot at 0, but no upper bound knot.)
    knots_tv <- c(0,stats::quantile(y[delta==1],quantile_seq))
  } else if(tolower(tv_type) %in% c("royston-parmar","rp")){
    #in rp model, boundary knots are set directly at endpoints, so no need to fix at 0
    #vector of quantiles at which we will set knots
    quantile_seq <- seq(from = 0,to = 1, length.out = p_tv)
    #rp model puts splines on log scale, but let's specify knots on original scale and
    #transform them in the get_basis function (makes the knots themselves more interpretable and consistent)
    knots_tv <- stats::quantile(y[delta==1],quantile_seq)
  } else {return(NULL)}
  return(knots_tv)
}

#' Get Basis Function Values for Flexible Hazard Specifications
#'
#' @param x Numeric vector of event times (e.g., \code{y1} or \code{y2}) at which
#'   to generate basis function values.
#' @param knots Increasing vector of integers corresponding to the knots used
#'   in the desired spline or piecewise specification. Often an element of
#'   list generated from \code{\link{get_default_knots_list}}.
#' @param deriv Boolean for whether returned values should be from derivatives of
#'   basis functions if \code{TRUE}, or original basis functions if \code{FALSE}.
#' @param flexsurv_compatible For Royston-Parmar basis, boolean for whether returned values should be untransformed,
#'   as in \code{flexsurv} package, if \code{TRUE}, or whether they should be transformed to remove correlation,
#'   as in ns package otherwise.
#'
#' @return A matrix with each row corresponding to an element of the input, and
#'   each column giving the corresponding basis function value.
#' @export
get_basis_tv <- function(x,knots,tv_type,
                         intercept=TRUE,
                         deriv=FALSE,
                         flexsurv_compatible=FALSE){

  #the exact form of the knots passed into this function comes from the above function

  if(tolower(tv_type) %in% c("bspline","bs")){
    if(deriv){return(NULL)} #we don't compute derivatives of bspline basis functions
    #compute spline basis, with boundary knots at 0 and maximum observed time, including an intercept.
    basis_out <- splines::bs(x = x,knots = knots[-c(1,length(knots))],Boundary.knots = knots[c(1,length(knots))],intercept = TRUE)
  } else if(tolower(tv_type) %in% c("piecewise","pw")){
    if(deriv){return(NULL)} #we don't compute derivatives of piecewise basis functions
    basis_out <- pw_cum_mat(y = x,knots = knots)
  } else if(tolower(tv_type) %in% c("royston-parmar","rp")){
    #lower "boundary" knot is set at minimum observed event times, not necessarily strictly 0
    knots_log <- log(knots) #for rp we generate basis on log scale, which means transforming knots and data
    if(any(is.infinite(knots_log) | is.na(knots_log))){stop("For royston-parmar model, all knots must be positive.")}
    temp_log <- log(x)
    temp_log[is.infinite(temp_log)] <- NA
    if(deriv){
      if(flexsurv_compatible){
        basis_out <- flexsurv::dbasis(x=temp_log,knots=knots_log)
        basis_out[is.na(basis_out)] <- 1 #can't set this to 0, because it is then logged and that causes a mess even when it multiplies with delta1delta2 and would otherwise be 0
        if(!intercept){basis_out <- basis_out[,-1]}
      } else{
        # basis_out <- splines2::naturalSpline(x = temp_log,
        #                                      knots = knots_log[-c(1,length(knots_log))],
        #                                      Boundary.knots = knots_log[c(1,length(knots_log))],
        #                                      intercept = intercept, derivs = 1)
        basis_out <- ns_d(x = temp_log,
                          knots = knots_log[-c(1,length(knots_log))],
                          Boundary.knots = knots_log[c(1,length(knots_log))],
                          intercept = intercept)
        basis_out[is.na(basis_out)] <- 1 #can't set this to 0, because it is then logged and that causes a mess even when it multiplies with delta1delta2 and would otherwise be 0
      }
    } else{
      if(flexsurv_compatible){
        # basis_out <- flexsurv::basis(x = temp_log,knots = knots_log)
        basis_out[is.na(basis_out)] <- 0 #This doesn't cause a problem for illness-death model because the changed rows are zeroed out of the likelihood by deltas
        if(!intercept){basis_out <- basis_out[,-1]}
      } else{
        # basis_out <- splines2::naturalSpline(x = temp_log,
        #                                      knots = knots_log[-c(1,length(knots_log))],
        #                                      Boundary.knots = knots_log[c(1,length(knots_log))],
        #                                      intercept = intercept, deriv = 0)
        basis_out <- splines::ns(x = temp_log,
                                 knots = knots_log[-c(1,length(knots_log))],
                                 Boundary.knots = knots_log[c(1,length(knots_log))],
                                 intercept = intercept)
        basis_out[is.na(basis_out)] <- 0 #This doesn't cause a problem for illness-death model because the changed rows are zeroed out of the likelihood by deltas
      }
    }
  } else {return(NULL)}

  return(basis_out)
}

#' 
#' #' Compute Covariate Process V under spline specification
#' #'
#' #' This function computes the integrated covariate process V used for time-varying AFT models.
#' #'
#' #' @param t_obj
#' #' @param x_base
#' #' @param beta_base
#' #' @param x_tv
#' #' @param beta_tv
#' #' @param xbeta_base
#' #' @param xbeta_tv
#' #' @param tv_type
#' #' @param basis
#' #' @param knots
#' #'
#' #' @return
#' #' @export
#' AFspline <- function(x_base=NULL, beta_base=NULL, xbeta_base=NULL,
#'                      x_tv=NULL, beta_tv=NULL, basis=NULL, knots=NULL){
#'   #browser()
#' 
#'   #I'll reimplement the need to manually compute the basis if needed
#'   # if(is.null(basis)){
#'   #   if(is.null(knots)){stop("Need to supply a basis, or a vector of knots.")}
#'   #   basis <- get_basis_tv(x=t_obj, knots=knots,flexsurv_compatible = TRUE)
#'   # }
#' 
#'   #compute vector of -x * spline
#'   temp <- as.numeric(-basis %*% beta_tv) * x_tv
#' 
#'   #Lastly, work out what to do with the basesline covariates
#'   if(is.null(xbeta_base)){
#'     if(is.null(x_base)){
#'       if(is.null(beta_base)){
#'         warning("we're assuming no baseline covariates, because xbeta_base, x_base, and beta_base are all NULL.")
#'         xbeta_base <- 0
#'       } else{
#'         warning("x_base and xbeta_base are both NULL, so we assume all baseline covariates are set to 1.")
#'         xbeta_base <- sum(beta_base)
#'       }
#'     } else{
#'       #calculate xbeta_base using x_base and beta_base
#'       xbeta_base <- x_base %*% as.matrix(beta_base)
#'     }
#'   }
#' 
#'   #returning the multiplicative factor on the cumulative scale
#'   as.numeric(exp(-xbeta_base - temp))
#' }
#' 
#' 
#' #' Title
#' #'
#' #' @param para
#' #' @param y
#' #' @param delta
#' #' @param x_base
#' #' @param x_tv
#' #' @param baseline
#' #' @param basis
#' #' @param dbasis
#' #'
#' #' @return
#' #' @export
#' nll_AFTtv_spline <- function (para, y, delta, x_base, x_tv,
#'                               baseline = "weibull", basis=NULL, dbasis=NULL){
#'   # browser()
#' 
#'   #number of baseline EFFECTS (including the baseline effect for a covariate that will also have a time-varying effect)
#'   nP_base <- ncol(x_base)
#'   #number of time-varying EFFECTS (excluding the baseline effect.)
#'   nP_tv <- if(is.null(basis)) 1 else ncol(basis)
#' 
#'   nP0 <- if(tolower(baseline) %in% c("weibull","lognormal")) 2 else 0
#'   nP_tot <- nP_base + nP_tv + nP0
#'   stopifnot(length(para) == nP_tot)
#' 
#'   intercept_temp <- para[1]
#'   logsigma_temp <- para[2]
#'   beta_base_temp <- if(nP_base > 0) para[(1+nP0):(nP0+nP_base)] else 0
#'   beta_tv_temp <- if(nP_tv > 0) para[(1+nP0+nP_base):(nP0+nP_base+nP_tv)] else 0
#' 
#'   xbeta_base_temp <- x_base %*% as.matrix(beta_base_temp)
#' 
#'   if(tolower(baseline)=="weibull"){
#'     logS0 <- function(x){stats::pweibull(q=x,scale = exp(intercept_temp),
#'                                          shape = exp(logsigma_temp),
#'                                          lower.tail = FALSE, log.p = TRUE)}
#'     logh0 <- function(x){flexsurv::hweibull(x=x,scale = exp(intercept_temp),
#'                                             shape = exp(logsigma_temp), log = TRUE)}
#'   } else if(tolower(baseline)=="lognormal"){
#'     logS0 <- function(x){stats::plnorm(q=x, meanlog = intercept_temp,
#'                                        sdlog = exp(logsigma_temp),
#'                                        lower.tail = FALSE, log.p = TRUE)}
#'     logh0 <- function(x){log(flexsurv::hlnorm(x=x, meanlog = intercept_temp,
#'                                               sdlog = exp(logsigma_temp)))} #logging manually bc of a bug in flexsurv for now
#'   }
#' 
#'   #now, to compute the relevant functions
#'   xspline_temp <- as.numeric(basis %*% beta_tv_temp) * x_tv
#'   dxspline_temp <- as.numeric(dbasis %*% beta_tv_temp) * x_tv
#'   V_temp <- y * exp(-xbeta_base_temp - xspline_temp)
#'   #derived via the chain rule from V
#'   logv <- - xbeta_base_temp - xspline_temp +
#'     log1p(-dxspline_temp)
#' 
#'   #This was what I had previously and I'm certain it was wrong
#'   # logv <- -xbeta_base_temp - xspline_temp - dxspline_temp# - y
#' 
#'   ll <- delta * (logh0(V_temp) + logv) + logS0(V_temp)
#'   mean(-ll)
#' }
#' 
#' 
#' 
#' ##**************************************##
#' ####FUNCTIONS FOR MARGINALIZED RESULTS####
#' ##**************************************##
#' 
#' 
#' 
#' #' Compute inverse of regression-standardized survival function
#' #'
#' #' @param p
#' #' @param x_base_aug
#' #' @param x_tv_aug
#' #' @param beta_base
#' #' @param beta_tv
#' #' @param int
#' #' @param shp
#' #' @param tv_type
#' #' @param knots
#' #' @param lower
#' #' @param upper
#' #'
#' #' @return
#' #' @export
#' Vx_spline_inv = function(t_obj, beta_base, beta_tv,knots=NULL,lower=0.01,upper=100){
#'   #for now, I'm assuming we're only having a single beta_base for an x_base = 1,
#'   #a single time t, and a spline specification for beta_tv with a single x_tv = 1
#'   #also, a lognormal baseline hazard by default for now.
#'   #browser()
#'   tryCatch(stats::uniroot(f = function (t){
#'     basis_spline <- get_basis_tv(x=t,knots=knots,tv_type = "rp",
#'                                  intercept = FALSE,deriv = FALSE,flexsurv_compatible = FALSE)
#'     t_obj - t * exp(-beta_base - basis_spline %*% beta_tv)},
#'     lower = lower, upper = upper)$root,
#'     error=function(e){return(NA)})
#' }
#' 
#' #' Compute regression-standardized survival function
#' #'
#' #' @param t
#' #' @param x_base_aug
#' #' @param x_tv_aug
#' #' @param beta_base
#' #' @param beta_tv
#' #' @param int
#' #' @param shp
#' #' @param tv_type
#' #' @param knots
#' #' @param baseline
#' #' @param ...
#' #'
#' #' @return
#' #' @export
#' Smarg_spline <- function(t, x_base_aug, x_tv_aug, beta_base, beta_tv,
#'                   int, shp, tv_type, basis,baseline="weibull",...){
#'   #aug prefix is a reminder that the intention is that this is a matrix with
#'   #the (time-varying) exposure 'set to' a value, and everything else marginalized out
#'   # browser()
#'   V_temp <- t * exp(- x_base_aug %*% beta_base - x_tv_aug * (basis %*% beta_tv))
#'   if(tolower(baseline) %in% c("wb","weibull")){
#'     S_temp_vec <- stats::pweibull(q=V_temp,
#'                                   scale = exp(int), shape = shp,lower.tail = FALSE)
#'   } else{
#'     S_temp_vec <- stats::plnorm(q=V_temp,
#'                                 meanlog = int, sdlog = shp,lower.tail = FALSE)
#'   }
#'   mean(S_temp_vec)
#' }
#' 
#' Smarg_spline_inv = function(p,x_base_aug, x_tv_aug, beta_base, beta_tv,
#'                      int, shp, tv_type, knots=NULL,baseline="weibull",lower=0,upper=100){
#'   # browser()
#'   tryCatch(stats::uniroot(f = function (t){
#'     basis_spline <- get_basis_tv(x=t,knots=knots,tv_type = "rp",
#'                                  intercept = FALSE,deriv = FALSE,flexsurv_compatible = FALSE)
#' 
#'     p - Smarg_spline(t=t, x_base_aug = x_base_aug, x_tv_aug=x_tv_aug,
#'             beta_base = beta_base, beta_tv = beta_tv,basis = basis_spline,
#'             int = int, shp = shp, baseline=baseline)},
#'       lower = lower, upper = upper)$root,
#'       error=function(e){return(NA)})
#' }
#' 
#' 
#' #### FUNCTIONS FOR TIME-VARYING COVARIATE ####
#' 
#' 
#' 
#' #' Title
#' #'
#' #' @param para
#' #' @param y
#' #' @param delta
#' #' @param x_base
#' #' @param x_tv
#' #' @param baseline
#' #' @param basis
#' #' @param dbasis
#' #'
#' #' @return
#' #' @export
#' nll_AFTtv_tvcov_spline <- function (para, y, delta, x_base, x_tv_time, x_tv,
#'                               baseline = "weibull", basis=NULL, dbasis=NULL){
#'   # browser()
#'   
#'   #number of baseline EFFECTS (EXCLUDING the baseline effect for a covariate that will also have a time-varying effect)
#'   #notice how this differs from above, because now we're putting all effects for the tvcov together.
#'   nP_base <- ncol(x_base)
#'   #number of time-varying EFFECTS (excluding the baseline effect.)
#'   nP_tv <- if(is.null(basis)) 1 else ncol(basis)
#'   
#'   nP0 <- if(tolower(baseline) %in% c("weibull","lognormal")) 2 else 0
#'   nP_tot <- nP_base + nP_tv + nP0
#'   stopifnot(length(para) == nP_tot)
#'   
#'   intercept_temp <- para[1]
#'   logsigma_temp <- para[2]
#'   beta_base_temp <- if(nP_base > 0) para[(1+nP0):(nP0+nP_base)] else 0
#'   beta_tv_temp <- if(nP_tv > 0) para[(1+nP0+nP_base):(nP0+nP_base+nP_tv)] else 0
#'   
#'   xbeta_base_temp <- x_base %*% as.matrix(beta_base_temp)
#'   
#'   if(tolower(baseline)=="weibull"){
#'     logS0 <- function(x){stats::pweibull(q=x,scale = exp(intercept_temp),
#'                                          shape = exp(logsigma_temp),
#'                                          lower.tail = FALSE, log.p = TRUE)}
#'     logh0 <- function(x){flexsurv::hweibull(x=x,scale = exp(intercept_temp),
#'                                             shape = exp(logsigma_temp), log = TRUE)}
#'   } else if(tolower(baseline)=="lognormal"){
#'     logS0 <- function(x){stats::plnorm(q=x, meanlog = intercept_temp,
#'                                        sdlog = exp(logsigma_temp),
#'                                        lower.tail = FALSE, log.p = TRUE)}
#'     logh0 <- function(x){log(flexsurv::hlnorm(x=x, meanlog = intercept_temp,
#'                                               sdlog = exp(logsigma_temp)))} #logging manually bc of a bug in flexsurv for now
#'   }
#'   
#'   #how much time is the person in the x_tv=0 state?
#'   base_time <- pmin(y,x_tv_time)
#'   #how much time is the person in the x_tv=1 state?
#'   tv_time <- pmax(0,y-base_time)
#'   
#'   if(is.null(x_tv)){
#'     # warning("x_tv was null, so we're assuming the time-varying covariate is a 0/1 jump for everyone.")
#'     x_tv <- 1
#'   }
#'   
#'   #now, to compute the relevant functions
#'   xspline_temp <- as.numeric(basis %*% beta_tv_temp) * x_tv
#'   xspline_temp[tv_time==0] <- 0 #just a placeholder value, because this doesn't get used in the nll anyways
#'   dxspline_temp <- as.numeric(dbasis %*% beta_tv_temp) * x_tv
#'   dxspline_temp[tv_time==0] <- 0
#'   V_temp <- exp(-xbeta_base_temp) * 
#'     ( base_time + tv_time * exp(-xspline_temp))
#'   #derived via the chain rule from V
#'   logv <- - xbeta_base_temp - xspline_temp +
#'     log1p(-dxspline_temp)
#'   
#'   #This was what I had previously and I'm certain it was wrong
#'   # logv <- -xbeta_base_temp - xspline_temp - dxspline_temp# - y
#'   
#'   ll <- delta * (logh0(V_temp) + logv) + logS0(V_temp)
#'   mean(-ll)
#' }
#' 
#' #' Compute inverse of regression-standardized survival function
#' #'
#' #' @param p
#' #' @param x_base_aug
#' #' @param x_tv_aug
#' #' @param beta_base
#' #' @param beta_tv
#' #' @param int
#' #' @param shp
#' #' @param tv_type
#' #' @param knots
#' #' @param lower
#' #' @param upper
#' #'
#' #' @return
#' #' @export
#' #' 
#' Vx_tvcov_spline_inv = function(t_obj, x_tv_time, x_tv, x_base=as.matrix(0), beta_base=0, beta_tv,knots=NULL,lower=0.01,upper=100){
#'   #for now, I'm assuming we're only having a single beta_base for an x_base = 1,
#'   #a single time t, and a spline specification for beta_tv with a single x_tv = 1
#'   #also, a lognormal baseline hazard by default for now.
#'   # browser()
#'   tryCatch(stats::uniroot(f = function (t){
#'     base_time <- pmin(t,x_tv_time)
#'     tv_time <- pmax(0,t-base_time)
#'     # print(paste("tv_time:",tv_time))
#'     if(tv_time==0){
#'       basis_spline <- t(rep(0,length(beta_tv)))
#'     } else{
#'       basis_spline <- get_basis_tv(x=tv_time,knots=knots,tv_type = "rp",
#'                                    intercept = TRUE,deriv = FALSE,flexsurv_compatible = FALSE)
#'     }
#'     t_obj - exp(- x_base %*% beta_base) * (base_time + tv_time * exp(-basis_spline %*% beta_tv * x_tv))},
#'     lower = lower, upper = upper)$root,
#'     error=function(e){return(NA)})
#' }
#' 

