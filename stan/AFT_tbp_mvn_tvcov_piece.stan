functions {

  //this simple code finds the index of the first knot larger than the value
  //we use this to 'evaluate' the time-varying effect at time y
  //knots must be in increasing order, and must start with 0 and end with infinity.
  //BUT FOR NOW THERE'S NO CHECK OF THIS!
  //it returns 0 if the outcome is in the first interval,
  //otherwise returns the index of the basis column corresponding to the interval its in
  int firstgreater(real y, vector knots){
    int ind = 1;
    while(y > knots[ind]){ //if y=0, then returns 0
      ind = ind + 1;
    }
    return ind-1;
  }


  //function giving the survival under TBP prior
  //Now, I'm gonna try and implement it in as numerically "safe" a way as I possibly can, working with logged things as much as possible.
  real tbp_lccdf(vector y, real scale, vector xbeta_inter, vector logw, vector log_beta_surv_const, int weib_ind){
    
    int n = num_elements(y);
    int J = num_elements(logw);

    real log_surv_base;
    real log_F_base;
    real log_running_DeltaS0;    
    vector[n] log_surv_tbp;

    // print("scale = ", scale);
    // print("log_beta_surv_const = ", log_beta_surv_const);
    // print("logw = ", logw);
    for(i in 1:n){
      //here, I think we actually can keep everything logged
      if(weib_ind){
        log_surv_base = weibull_lccdf(y[i] | scale, exp(xbeta_inter[i]));      
//        log_F_base = weibull_lcdf(y[i] | scale, exp(xbeta_inter[i]));
      } else{
        log_surv_base = lognormal_lccdf(y[i] | xbeta_inter[i], scale);
//        log_F_base = lognormal_lcdf(y[i] | xbeta_inter[i], scale);
      }
      log_F_base = log1m_exp(log_surv_base);

      log_running_DeltaS0 = log1m_exp(J * log_F_base);
      log_surv_tbp[i] = logw[1] + log_running_DeltaS0;
      // print("person = ", i);
      // print("y = ", y[i]);
      // print("xbeta_inter = ", xbeta_inter[i]);
      // print("log_surv_base = ", log_surv_base);
      // print("log_F_base = ", log_F_base);
      // print("log_running_DeltaS0 after iteration 1 = ", log_running_DeltaS0);
      // print("log_surv_tbp after iteration 1 = ", log_surv_tbp[i]);
      for(j in 2:J){
        log_running_DeltaS0 = log_diff_exp(log_running_DeltaS0, log_beta_surv_const[j-1] + (j-1) * log_surv_base + (J - j + 1) * log_F_base);
        // print("log_running_DeltaS0 after next iteration = ", log_running_DeltaS0);
        if(is_nan(log_running_DeltaS0)) break;
        log_surv_tbp[i] = log_sum_exp(log_surv_tbp[i], logw[j] + log_running_DeltaS0);
        // print("log_surv_tbp after next iteration = ", log_surv_tbp[i]);
      }
    }

    return sum(log_surv_tbp);
  }

  //function giving the density under TBP prior
  real tbp_lpdf(vector y, real scale, vector xbeta_inter, vector logwconst, int weib_ind){

    int n = num_elements(y);
    int J = num_elements(logwconst);

    real log_surv_base;
    real log_F_base;
    real log_dens_base;    
    vector[n] log_dens_tbp;
    
    for(i in 1:n){

      //here, I think we actually can keep everything logged
      if(weib_ind){
        log_surv_base = weibull_lccdf(y[i] | scale, exp(xbeta_inter[i]));      
        //log_F_base = weibull_lcdf(y[i] | scale, exp(xbeta_inter[i]));
      } else{
        log_surv_base = lognormal_lccdf(y[i] | xbeta_inter[i], scale);
        //log_F_base = lognormal_lcdf(y[i] | xbeta_inter[i], scale);
      }
      log_F_base = log1m_exp(log_surv_base);

      log_dens_tbp[i] = logwconst[1] + (J-1) * log_F_base;
      // log f_tbp = log( sum_j( exp( log_w_j + beta_lpdf(S_0, j) ) ) ) 
      // for numerical stability and speed we're instead retaining the underlying log S_0     
      for(j in 2:J){
        log_dens_tbp[i] = log_sum_exp(log_dens_tbp[i], logwconst[j] + (j - 1) * log_surv_base + (J-j) * log_F_base);
      }
    }

    //additionally compute the log density
    if(weib_ind){
      log_dens_base = weibull_lpdf(y | scale, exp(xbeta_inter));
    } else{
      log_dens_base = lognormal_lpdf(y | xbeta_inter, scale);
    }

    return log_dens_base + sum(log_dens_tbp);  
  }

  //NOW, VERSIONS FOR INDIVIDUAL LOG LIKELIHOOD CONTRIBUTIONS

  real tbp_lccdf_i(real y_i, real scale, real xbeta_inter_i, vector logw, vector log_beta_surv_const, int weib_ind){

    int J = num_elements(logw);
    real log_surv_base;
    real log_F_base;
    real log_running_DeltaS0;    
    real log_surv_tbp;

    if(weib_ind){
      log_surv_base = weibull_lccdf(y_i | scale, exp(xbeta_inter_i));      
      //log_F_base = weibull_lcdf(y_i | scale, exp(xbeta_inter_i));
    } else{
      log_surv_base = lognormal_lccdf(y_i | xbeta_inter_i, scale);
      //log_F_base = lognormal_lcdf(y_i | xbeta_inter_i, scale);
    }
    log_F_base = log1m_exp(log_surv_base);

    log_running_DeltaS0 = log1m_exp(J * log_F_base);
    log_surv_tbp = logw[1] + log_running_DeltaS0;
    for(j in 2:J){
      log_running_DeltaS0 = log_diff_exp(log_running_DeltaS0, log_beta_surv_const[j-1] + (j-1) * log_surv_base + (J - j + 1) * log_F_base);
      if(is_nan(log_running_DeltaS0)) break;
      log_surv_tbp = log_sum_exp(log_surv_tbp, logw[j] + log_running_DeltaS0);
    }
    return log_surv_tbp;
  }


  //function giving the density under TBP prior
  real tbp_lpdf_i(real y_i, real scale, real xbeta_inter_i, vector logwconst, int weib_ind){

    int J = num_elements(logwconst);
    real log_surv_base;
    real log_F_base;
    real log_dens;
    
    //here, I think we actually can keep everything logged
    if(weib_ind){
      log_surv_base = weibull_lccdf(y_i | scale, exp(xbeta_inter_i));      
      //log_F_base = weibull_lcdf(y_i | scale, exp(xbeta_inter_i));
    } else{
      log_surv_base = lognormal_lccdf(y_i | xbeta_inter_i, scale);
      //log_F_base = lognormal_lcdf(y_i | xbeta_inter_i, scale);
    }
    log_F_base = log1m_exp(log_surv_base);

    log_dens = logwconst[1] + (J-1) * log_F_base;
    // log f_tbp = log( sum_j( exp( log_w_j + beta_lpdf(S_0, j) ) ) ) 
    // for numerical stability and speed we're instead retaining the underlying log S_0     
    for(j in 2:J){
      log_dens = log_sum_exp(log_dens, logwconst[j] + (j - 1) * log_surv_base + (J-j) * log_F_base);
    }

    //additionally compute the log density
    if(weib_ind){
      log_dens += weibull_lpdf(y_i | scale, exp(xbeta_inter_i));
    } else{
      log_dens += lognormal_lpdf(y_i | xbeta_inter_i, scale);
    }

    return log_dens;
  }


}

data {
  //baseline specification
  int<lower=0> weib_ind; //indicator for weibull (1) vs log-normal (0)

  //covariate specification
  int<lower=0> nP_base; // number of baseline beta parameters
  int<lower=0> nP_tv_piece; // number of time-varying beta parameters

  vector[nP_tv_piece + 1] knots_tv; //vector of knots, plus a 0 at the start and an infinity at the end.
  
  // data for censored subjects
  int<lower=0> N_m; //number of censored obs
  matrix[N_m,nP_base] X_m; //matrix of covariates for censored obs
  vector[N_m] Xtv_m;  //vector of time_varying covariate
  vector[N_m] Xtv_time_m;  //vector of times at which Xtv jumped from 0 to the level Xtv
  vector[N_m] y_m; // vector of censoring times.
  matrix[N_m,nP_tv_piece] basis_piece_m; //matrix of nP_tv basis functions evaluated at y
  int L_ind_m[N_m];
  
  // data for observed subjects
  int<lower=0> N_o;
  matrix[N_o,nP_base] X_o;
  vector[N_o] Xtv_o;  
  vector[N_o] Xtv_time_o;  //vector of times at which Xtv jumped from 0 to the level Xtv
  vector[N_o] y_o; // vector of observed event times.
  matrix[N_o,nP_tv_piece] basis_piece_o; //matrix of nP_tv basis functions evaluated at y
  int L_ind_o[N_o];

  // data for left truncation observations
  int<lower=0> N_l;
  matrix[N_l,nP_base] X_l;
  vector[N_l] Xtv_l;  
  vector[N_l] Xtv_time_l;  //vector of times at which Xtv jumped from 0 to the level Xtv
  vector[N_l] L_l; // vector of observed event times.
  matrix[N_l,nP_tv_piece] basis_piece_l; //matrix of nP_tv basis functions evaluated at y
  
  //inputs for model
  int<lower=1> J; //number of basis functions for bernstein polynomial

  //hyperparameters (assuming that we put an mvn prior on the centering parameters)
  vector[2] prior_mean_centering;
  matrix[2,2] prior_prec_centering;
  int<lower=0> prior_beta; //indicator for whether the prior should be applied for beta
  real<lower=0> m_beta;
  real<lower=0> sd_beta;
  int<lower=0> prior_beta_tv; //indicator for whether the prior should be applied for beta_tv
  real<lower=0> m_beta_tv;
  real<lower=0> sd_beta_tv;

  int<lower = 0, upper = 1> dirichlet_alpha_fixed; //indicator for whether the prior should be applied for dirichlet alpha param
  real<lower=0> dirichlet_alpha_data[dirichlet_alpha_fixed ? 1 : 0]; //dirichlet_alpha_data is size 0 if dirichlet_alpha_fixed is FALSE
  real<lower=0> a_dirichlet_alpha;
  real<lower=0> b_dirichlet_alpha;

}

transformed data {
  vector[J] log_beta_surv_const;
  vector[J] log_beta_dens_const;

  //basetime is the amount of time accrued before the observed jump in x_tv
  //it is 0 unless x_tv is a time-varying covariate
  vector[N_m] base_time_m;
  vector[N_o] base_time_o;
  vector[N_l] base_time_l;

  //tvtime is the amount of time accrued after the observed jump in x_tv
  vector[N_m] tv_time_m;
  vector[N_o] tv_time_o;
  vector[N_l] tv_time_l;

  //cut_cat gives which interval of the piecewise effect the observed subject's time falls.
  //set it to 0 if x doesn't jump at all.
  int cut_cat_m[N_m];
  int cut_cat_o[N_o];
  int cut_cat_l[N_l];

  for (i in 1:N_o) {
    base_time_o[i] = fmin(y_o[i],Xtv_time_o[i]);
    tv_time_o[i] = fmax(0,y_o[i]-base_time_o[i]);
    cut_cat_o[i] = firstgreater(tv_time_o[i], knots_tv);
  }

  for (j in 1:N_m) {
    base_time_m[j] = fmin(y_m[j],Xtv_time_m[j]);
    tv_time_m[j] = fmax(0,y_m[j]-base_time_m[j]);
    cut_cat_m[j] = firstgreater(tv_time_m[j], knots_tv);
  }

  //how to handle time-varying covariates with left truncation??
  //two options:
  // 1. assume that people JUST jumped at the time of study entry
  // 2. assume that people jumped at the time origin
  //we get 1 if Xtv_time_l is set to L_l, and 2 if Xtv_time_l is set to 0.
    //BUT BEWARE THIS SHOULD MATCH WHAT WAS DONE IN THE BASIS CONSTRUCTION
  for (k in 1:N_l) {
    base_time_l[k] = fmin(L_l[k],Xtv_time_l[k]);
    tv_time_l[k] = fmax(0,L_l[k]-base_time_l[k]);
    cut_cat_l[k] = firstgreater(tv_time_l[k], knots_tv);
  }  


  for(j in 1:J){
    log_beta_surv_const[j] = -lbeta(j * 1.0, J - j + 1.0) - log(j * 1.0);
    log_beta_dens_const[j] = -lbeta(j * 1.0, J - j + 1.0);
  }

  // print("log beta survival constants = ", log_beta_surv_const);
  // print("log beta density constants = ", log_beta_dens_const);

}

parameters {
  vector[2] centering; //intercept and scale
  
  vector[nP_base] beta; //baseline covariates
  vector[nP_tv_piece] beta_tv; //time-varying effect

  simplex[J] w; // coefficients for basis functions

  //dirichlet_alpha_param is size 0 if dirichlet_alpha_fixed is TRUE
  real<lower=0> dirichlet_alpha_param[dirichlet_alpha_fixed ? 0 : 1];

}

transformed parameters{
  //mvn prior parameters for intercept and logscale parameters
  real intercept = centering[1];
  real<lower=0> scale = exp(centering[2]);

  vector[N_m] xbeta_inter_m;
  vector[N_o] xbeta_inter_o;
  vector[N_l] xbeta_inter_l;

  vector[N_o] Vout_o;
  vector[N_m] Vout_m;
  vector[N_l] Vout_l;

  vector[N_o] logveval_o; //this is the lower case v function evaluated at each y_o


  vector[J] logw;
  real<lower=0> dirichlet_alpha;

  if(N_o > 0){
    if(nP_base > 0){
      xbeta_inter_o = intercept + X_o*beta;
    } else{
      xbeta_inter_o = to_vector(rep_array(intercept,N_o));    
    }
  }
  if(N_m > 0){
    if(nP_base > 0){
      xbeta_inter_m = intercept + X_m*beta;
    } else{
      xbeta_inter_m = to_vector(rep_array(intercept,N_m));    
    }
  }
  if(N_l > 0){
    if(nP_base > 0){
      xbeta_inter_l = intercept + X_l*beta;
    } else{
      xbeta_inter_l = to_vector(rep_array(intercept,N_l));    
    }
  }


  //this is the definition of the V function corresponding to the time AFTER the covariate jumps:
  //basis_rowvec * exp(-x_tv*beta_tv);

  //the basis formula parameterizes relative to the first interval
  //so, first beta_tv is effect in second interval relative to first, second is in third interval relative to first, etc.
  
  //generate the integrated value used in the likelihood for observed data
  for (i in 1:N_o) {
    if(cut_cat_o[i] == 0){
      logveval_o[i] = 0;
      Vout_o[i] = base_time_o[i];
    } else {
      logveval_o[i] = -Xtv_o[i] * beta_tv[ cut_cat_o[i] ];
      Vout_o[i] = base_time_o[i] + basis_piece_o[i] * exp(-Xtv_o[i]*beta_tv);
    }
  }
  
  //generate the integrated value used in the likelihood for censored obs
  for (j in 1:N_m) {
    if(cut_cat_m[j] == 0){
      Vout_m[j] = base_time_m[j];
    } else {
      Vout_m[j] = base_time_m[j] + basis_piece_m[j] * exp(-Xtv_m[j]*beta_tv);
    }
  }

  //generate the integrated value used in the likelihood for left-truncated obs
  for (k in 1:N_l) {
    if(cut_cat_l[k] == 0){
      Vout_l[k] = base_time_l[k];
    } else {
      Vout_l[k] = base_time_l[k] + basis_piece_l[k] * exp(-Xtv_l[k]*beta_tv);
    }
  }

  //choose which dirichlet hyperparameter should be used accordingly
  if (dirichlet_alpha_fixed) {
    dirichlet_alpha = dirichlet_alpha_data[1];
  } else {
    dirichlet_alpha = dirichlet_alpha_param[1];
  }

  logw = log(w);

  // print("dirichlet_alpha = ", dirichlet_alpha);
  // print("logw = ", logw);
  // print("scale = ", scale);
  // print("beta = ", beta);
  // print("beta_tv = ", beta_tv);

}

model {

  //mvn prior on centering distribution parameters
  centering ~ multi_normal_prec(prior_mean_centering,prior_prec_centering);

  if(prior_beta==1){
    beta ~ normal(m_beta, sd_beta^2);    
  }
  if(prior_beta_tv==1){
    beta_tv ~ normal(m_beta_tv, sd_beta_tv^2);    
  }

  //gamma hyperprior on dirichlet parameter, following Hanson paper
  if(!dirichlet_alpha_fixed){
    dirichlet_alpha_param ~ gamma(a_dirichlet_alpha,b_dirichlet_alpha); //hyperprior for dirichlet
  }
  w ~ dirichlet(rep_vector(dirichlet_alpha,J));

  // print("log posterior before log likelihood = ", target());
  
  // evaluate likelihood for censored and uncensored subjects
  if(N_o > 0){
    target += tbp_lpdf(Vout_o| scale, xbeta_inter_o, logw + log_beta_dens_const, weib_ind);
    // print("log posterior after observed tbp pdf = ", target());
    target += logveval_o;    
    // print("log posterior after tv correction for observed events = ", target());
  }
  if(N_m > 0){
    target += tbp_lccdf(Vout_m| scale, xbeta_inter_m, logw, log_beta_surv_const, weib_ind);
    // print("log posterior after censored events = ", target());
  }
  if(N_l > 0){
    target += -tbp_lccdf(Vout_l| scale, xbeta_inter_l, logw, log_beta_surv_const, weib_ind);    
    // print("log posterior after left truncation events = ", target());
  }

}

generated quantities {

  vector[N_o + N_m] log_lik;

  for (i in 1:N_o) {
    log_lik[i] = 0;
    log_lik[i] += tbp_lpdf_i(Vout_o[i], scale, xbeta_inter_o[i], logw + log_beta_dens_const, weib_ind);

    if(N_l > 0){
      if(L_ind_o[i] != -1){
        log_lik[i] += -tbp_lccdf_i(Vout_l[L_ind_o[i]], scale, xbeta_inter_l[L_ind_o[i]], logw, log_beta_surv_const, weib_ind);
      }
    }
    log_lik[i] += logveval_o[i];
  }

  for (j in 1:N_m) {
    log_lik[N_o + j] = 0;
    log_lik[N_o + j] += tbp_lccdf_i(Vout_m[j], scale, xbeta_inter_m[j], logw, log_beta_surv_const, weib_ind);      

    if(N_l > 0){
      if(L_ind_m[j] != -1){
        log_lik[N_o + j] += -tbp_lccdf_i(Vout_l[L_ind_m[j]], scale, xbeta_inter_l[L_ind_m[j]], logw, log_beta_surv_const, weib_ind);
      }
    }
  }

}


