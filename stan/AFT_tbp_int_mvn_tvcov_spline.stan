functions {

  //function giving the survival under TBP prior
  //Now, I'm gonna try and implement it in as numerically "safe" a way as I possibly can, working with logged things as much as possible.
  real tbp_lccdf(vector y, real scale, vector xbeta_inter, vector logw, vector log_beta_surv_const, int weib_ind){
    
    int n = num_elements(y);
    int J = num_elements(logw);

    real log_surv_base;
    real log_F_base;
    real log_running_DeltaS0;    
    vector[n] log_surv_tbp;

    for(i in 1:n){
      //here, I think we actually can keep everything logged
      if(weib_ind){
        log_surv_base = weibull_lccdf(y[i] | scale, exp(xbeta_inter[i]));      
        // log_F_base = weibull_lcdf(y[i] | scale, exp(xbeta_inter[i]));
      } else{
        log_surv_base = lognormal_lccdf(y[i] | xbeta_inter[i], scale);
        // log_F_base = lognormal_lcdf(y[i] | xbeta_inter[i], scale);
      }
      log_F_base = log1m_exp(log_surv_base);

      log_running_DeltaS0 = log1m_exp(J * log_F_base);
      log_surv_tbp[i] = logw[1] + log_running_DeltaS0;
      for(j in 2:J){
        log_running_DeltaS0 = log_diff_exp(log_running_DeltaS0, log_beta_surv_const[j-1] + (j-1) * log_surv_base + (J - j + 1) * log_F_base);
        if(is_nan(log_running_DeltaS0)) break;
        log_surv_tbp[i] = log_sum_exp(log_surv_tbp[i], logw[j] + log_running_DeltaS0);
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
        // log_F_base = weibull_lcdf(y[i] | scale, exp(xbeta_inter[i]));
      } else{
        log_surv_base = lognormal_lccdf(y[i] | xbeta_inter[i], scale);
        // log_F_base = lognormal_lcdf(y[i] | xbeta_inter[i], scale);
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
      // log_F_base = weibull_lcdf(y_i | scale, exp(xbeta_inter_i));
    } else{
      log_surv_base = lognormal_lccdf(y_i | xbeta_inter_i, scale);
      // log_F_base = lognormal_lcdf(y_i | xbeta_inter_i, scale);
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
      // log_F_base = weibull_lcdf(y_i | scale, exp(xbeta_inter_i));
    } else{
      log_surv_base = lognormal_lccdf(y_i | xbeta_inter_i, scale);
      // log_F_base = lognormal_lcdf(y_i | xbeta_inter_i, scale);
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

  int<lower=0> nP_base; // number of baseline beta parameters (including the baseline component of a time-varying effect)
  int<lower=0> nP_tv_spline; // number of time-varying beta parameters (excluding the baseline component of a time-varing effect)
  
  // data for censored subjects
  int<lower=0> N_m; //number of censored obs
  matrix[N_m,nP_base] X_m; //matrix of covariates for censored obs
  vector[N_m] Xtv_m;  //vector of covariate with time-varying effect
  vector[N_m] Xtv_time_m;  //vector of times at which Xtv jumped from 0 to the level Xtv
  vector[N_m] y_m; // vector of censoring times.
  matrix[N_m,nP_tv_spline] basis_spline_m; //matrix of nP_tv basis functions evaluated at y
  int L_ind_m[N_m];
  
  // data for observed subjects
  int<lower=0> N_o;
  matrix[N_o,nP_base] X_o;
  vector[N_o] Xtv_o;  
  vector[N_o] Xtv_time_o;  //vector of times at which Xtv jumped from 0 to the level Xtv
  vector[N_o] y_o; // vector of observed event times.
  matrix[N_o,nP_tv_spline] basis_spline_o; //matrix of nP_tv basis functions evaluated at y
  matrix[N_o,nP_tv_spline] dbasis_spline_o; //matrix of nP_tv dbasis functions evaluated at y
  int L_ind_o[N_o];

  // data for left truncation observations
  int<lower=0> N_l;
  matrix[N_l,nP_base] X_l;
  vector[N_l] Xtv_l;  
  vector[N_l] Xtv_time_l;  //vector of times at which Xtv jumped from 0 to the level Xtv
  vector[N_l] L_l; // vector of observed event times.
  matrix[N_l,nP_tv_spline] basis_spline_l; //matrix of nP_tv basis functions evaluated at y  
  
  // data for interval-censored subjects
  int<lower=0> N_i; //number of censored obs
  matrix[N_i,nP_base] X_i; //matrix of covariates for censored obs
  vector[N_i] Xtv_i;  //vector of time_varying covariate
  vector[N_i] Xtv_time_i;  //vector of times at which Xtv jumped from 0 to the level Xtv
  vector[N_i] y_iL; // vector of left interval-censoring times.
  vector[N_i] y_iU; // vector of right interval-censoring times.
  matrix[N_i,nP_tv_spline] basis_spline_iL; //matrix of nP_tv basis functions evaluated at yL
  matrix[N_i,nP_tv_spline] basis_spline_iU; //matrix of nP_tv basis functions evaluated at yU
  int L_ind_i[N_i];

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
  vector[N_m] base_time_m;
  vector[N_o] base_time_o;
  vector[N_i] base_time_i;
  vector[N_l] base_time_l;

  //tvtime is the amount of time accrued after the observed jump in x_tv
  vector[N_m] tv_time_m;
  vector[N_o] tv_time_o;
  vector[N_i] tv_time_iL;
  vector[N_i] tv_time_iU;
  vector[N_l] tv_time_l;

  //a new version of Xtv with value set to 0 if y <= x_time.
  vector[N_o] Xtv_nonzero_o;
  vector[N_m] Xtv_nonzero_m;
  vector[N_i] Xtv_nonzero_i;
  vector[N_l] Xtv_nonzero_l;


  //Here I'm using a few helpful expressions:
  //https://mc-stan.org/docs/2_27/reference-manual/conditional-operator-section.html
  //https://mc-stan.org/docs/2_28/functions-reference/absolute-functions.html

  for (i in 1:N_o) {
    base_time_o[i] = fmin(y_o[i],Xtv_time_o[i]);
    tv_time_o[i] = fmax(0,y_o[i]-base_time_o[i]);
    Xtv_nonzero_o[i] = int_step(tv_time_o[i]) ? Xtv_o[i] : 0;
  }

  for (j in 1:N_m) {
    base_time_m[j] = fmin(y_m[j],Xtv_time_m[j]);
    tv_time_m[j] = fmax(0,y_m[j]-base_time_m[j]);
    Xtv_nonzero_m[j] = int_step(tv_time_m[j]) ? Xtv_m[j] : 0;
  }

  for (l in 1:N_i) {
    base_time_i[l] = fmin(y_iL[l],Xtv_time_i[l]);
    tv_time_iL[l] = fmax(0,y_iL[l]-base_time_i[l]);
    tv_time_iU[l] = fmax(0,y_iU[l]-base_time_i[l]);
    Xtv_nonzero_i[l] = int_step(tv_time_iL[l]) ? Xtv_i[l] : 0;
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
    Xtv_nonzero_l[k] = int_step(tv_time_l[k]) ? Xtv_l[k] : 0;
  }

  for(j in 1:J){
    log_beta_surv_const[j] = -lbeta(j * 1.0, J - j + 1.0) - log(j * 1.0);
    log_beta_dens_const[j] = -lbeta(j * 1.0, J - j + 1.0);
  }


}

parameters {
  vector[2] centering; //intercept and scale

  vector[nP_base] beta; //baseline covariates
  vector[nP_tv_spline] beta_tv; //time-varying effect

  simplex[J] w; // coefficients for basis functions

  //dirichlet_alpha_param is size 0 if dirichlet_alpha_fixed is TRUE
  real<lower=0> dirichlet_alpha_param[dirichlet_alpha_fixed ? 0 : 1];

}

transformed parameters{
  //mvn prior parameters for intercept and logscale parameters
  real intercept = centering[1];
  real<lower=0> scale = exp(centering[2]);

  // transform the vectors of covariates
  vector[N_m] xbeta_inter_m;
  vector[N_o] xbeta_inter_o;
  vector[N_i] xbeta_inter_i;
  vector[N_l] xbeta_inter_l;

  vector[N_o] Vout_o;
  vector[N_m] Vout_m;
  vector[N_i] Vout_iL;
  vector[N_i] Vout_iU;
  vector[N_l] Vout_l;

  vector[N_o] xspline_o;
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
  if(N_i > 0){
    if(nP_base > 0){
      xbeta_inter_i = intercept + X_i*beta;
    } else{
      xbeta_inter_i = to_vector(rep_array(intercept,N_i));    
    }
  }  
  if(N_l > 0){
    if(nP_base > 0){
      xbeta_inter_l = intercept + X_l*beta;
    } else{
      xbeta_inter_l = to_vector(rep_array(intercept,N_l));    
    }
  }

  //we use this twice when computing things so let's store it
  //we multiply terms that only arise after the covariate jump by tv_nonzero to set them to 0 otherwise.
  if(N_o > 0){
    xspline_o = (basis_spline_o * beta_tv) .* Xtv_nonzero_o; //assume that this is 0 when there is no time accrued post-jump
    Vout_o = base_time_o + tv_time_o .* exp(-xspline_o);    
    logveval_o = -xspline_o + log1p(-(dbasis_spline_o * beta_tv) .* Xtv_nonzero_o);
  }
  if(N_m > 0){
    Vout_m = base_time_m + tv_time_m .* exp(-(basis_spline_m * beta_tv) .* Xtv_nonzero_m);
  }
  if(N_i > 0){
    Vout_iL = base_time_i + tv_time_iL .* exp(-(basis_spline_iL * beta_tv) .* Xtv_nonzero_i);
    Vout_iU = base_time_i + tv_time_iU .* exp(-(basis_spline_iU * beta_tv) .* Xtv_nonzero_i);
  }
  if(N_l > 0){
    Vout_l = base_time_l + tv_time_l .* exp(-(basis_spline_l * beta_tv) .* Xtv_nonzero_l);
  }

  //choose which dirichlet hyperparameter should be used accordingly
  if (dirichlet_alpha_fixed) {
    dirichlet_alpha = dirichlet_alpha_data[1];
  } else {
    dirichlet_alpha = dirichlet_alpha_param[1];
  }

  logw = log(w);

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

  // evaluate likelihood for censored and uncensored subjects
  if(N_o > 0){
    target += tbp_lpdf(Vout_o| scale, xbeta_inter_o, logw + log_beta_dens_const, weib_ind);
    target += logveval_o;
  }
  if(N_m > 0){
    target += tbp_lccdf(Vout_m| scale, xbeta_inter_m, logw, log_beta_surv_const, weib_ind);
  }
  if(N_l > 0){
    target += -tbp_lccdf(Vout_l| scale, xbeta_inter_l, logw, log_beta_surv_const, weib_ind);    
  }
  for (l in 1:N_i) {
    target += log_diff_exp( tbp_lccdf_i(Vout_iL[l], scale, xbeta_inter_i[l], logw, log_beta_surv_const, weib_ind),
                            tbp_lccdf_i(Vout_iU[l], scale, xbeta_inter_i[l], logw, log_beta_surv_const, weib_ind));
  }

}

generated quantities {

  vector[N_o + N_m + N_i] log_lik;

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

  for (l in 1:N_i) {
    log_lik[N_o + N_m + l] = 0;
    log_lik[N_o + N_m + l] += log_diff_exp( tbp_lccdf_i(Vout_iL[l], scale, xbeta_inter_i[l], logw, log_beta_surv_const, weib_ind),
                                            tbp_lccdf_i(Vout_iU[l], scale, xbeta_inter_i[l], logw, log_beta_surv_const, weib_ind));

    if(N_l > 0){
      if(L_ind_i[l] != -1){
        log_lik[N_o + N_m + l] += -tbp_lccdf_i(Vout_l[L_ind_i[l]], scale, xbeta_inter_l[L_ind_i[l]], logw, log_beta_surv_const, weib_ind);
      }
    }
  }

}


