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
  # browser()
  
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
  
  # #somehow, for computing marginal af this is actually slower lol...
  # knots_mat <- splines2::ibs(y,knots=knots[-1], degree = 0,
  #                            Boundary.knots = c(knots[1],1e305), intercept = TRUE)
  
  
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
