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

  //hyperparameters
  int<lower=0> prior_scale; //indicator for whether the prior should be applied for scale
  real<lower=0> a_scale;
  real<lower=0> b_scale;
  int<lower=0> prior_intercept; //indicator for whether the prior should be applied for mu
  real<lower=0> m_intercept;
  real<lower=0> sd_intercept;
  int<lower=0> prior_beta; //indicator for whether the prior should be applied for beta
  real<lower=0> m_beta;
  real<lower=0> sd_beta;
  int<lower=0> prior_beta_tv; //indicator for whether the prior should be applied for beta_tv
  real<lower=0> m_beta_tv;
  real<lower=0> sd_beta_tv;

}

transformed data{
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

}
parameters {
  real intercept;  
  real<lower=0> scale_raw; // Weibull Shape    

  vector[nP_base] beta; //baseline covariates
  vector[nP_tv_spline] beta_tv; //time-varying effect

}
 
transformed parameters{

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

  real<lower=0> scale;

  scale = (prior_scale == 1) ? sqrt(scale_raw) : scale_raw ;  

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
  
}

model {

  //match the prior used in the MCMC sampler of SemiCompRisks
  if(prior_scale==1){ //inverse gamma
    scale_raw ~ inv_gamma(a_scale, b_scale);    
  } else if(prior_scale==2){ //gamma
    scale_raw ~ gamma(a_scale,b_scale);
  }
  if(prior_intercept==1){
    intercept ~ normal(m_intercept, sd_intercept^2);    
  }
  if(prior_beta==1){
    beta ~ normal(m_beta, sd_beta^2);    
  }
  if(prior_beta_tv==1){
    beta_tv ~ normal(m_beta_tv, sd_beta_tv^2);    
  }
  
  if(weib_ind == 1){ //weibull baseline
    if(N_o > 0){
      target += weibull_lpdf(Vout_o | scale, exp(xbeta_inter_o) );
    }
    if(N_m > 0){
      target += weibull_lccdf(Vout_m | scale, exp(xbeta_inter_m) );      
    }
    if(N_l > 0){
      target += -weibull_lccdf(Vout_l | scale, exp(xbeta_inter_l) );
    }
    //calculate interval contributions one by one on the log scale
    for (l in 1:N_i) {
      target += log_diff_exp( weibull_lccdf(Vout_iL[l] | scale, exp(xbeta_inter_i[l]) ) ,
                              weibull_lccdf(Vout_iU[l] | scale, exp(xbeta_inter_i[l]) ));
    }    
  } else{ //lognormal baseline, with "scale" square rooted (to naturally comport with inv-gamma
    if(N_o > 0){
      target += lognormal_lpdf(Vout_o | xbeta_inter_o, scale);
    }
    if(N_m > 0){
      target += lognormal_lccdf(Vout_m | xbeta_inter_m, scale);
    }
    if(N_l > 0){
      target += -lognormal_lccdf(Vout_l | xbeta_inter_l, scale);
    }
    //calculate interval contributions one by one on the log scale
    for (l in 1:N_i) {
      target += log_diff_exp( lognormal_lccdf(Vout_iL[l] | xbeta_inter_i[l], scale ) ,
                              lognormal_lccdf(Vout_iU[l] | xbeta_inter_i[l], scale ));
    }
  }

  if(N_o > 0){
    target += logveval_o;
  }

}

generated quantities {

  vector[N_o + N_m + N_i] log_lik;

  for (i in 1:N_o) {
    log_lik[i] = 0;
    if(weib_ind == 1){
      log_lik[i] += weibull_lpdf(Vout_o[i] | scale, exp(xbeta_inter_o[i]) );
    } else{
      log_lik[i] += lognormal_lpdf(Vout_o[i] | xbeta_inter_o[i], scale);
    }

    if(N_l > 0){
      if(L_ind_o[i] != -1){
        if(weib_ind == 1){
          log_lik[i] += -weibull_lccdf(Vout_l[L_ind_o[i]] | scale, exp(xbeta_inter_l[L_ind_o[i]]));
        } else{
          log_lik[i] += -lognormal_lccdf(Vout_l[L_ind_o[i]] | xbeta_inter_l[L_ind_o[i]], scale);          
        }
      }      
    }
    log_lik[i] += logveval_o[i];
  }

  for (j in 1:N_m) {
    log_lik[N_o + j] = 0;
    if(weib_ind == 1){
      log_lik[N_o + j] += weibull_lccdf(Vout_m[j] | scale, exp(xbeta_inter_m[j]));      
    } else{
      log_lik[N_o + j] += lognormal_lccdf(Vout_m[j] | xbeta_inter_m[j], scale);
    }

    if(N_l > 0){
      if(L_ind_m[j] != -1){
        if(weib_ind == 1){
          log_lik[N_o + j] += -weibull_lccdf(Vout_l[L_ind_m[j]] | scale, exp(xbeta_inter_l[L_ind_m[j]]));
        } else{
          log_lik[N_o + j] += -lognormal_lccdf(Vout_l[L_ind_m[j]] | xbeta_inter_l[L_ind_m[j]], scale);          
        }
      }
    }
  }

  for (l in 1:N_i) {
    log_lik[N_o + N_m + l] = 0;
    if(weib_ind == 1){
      log_lik[N_o + N_m + l] = log_diff_exp( weibull_lccdf(Vout_iL[l] | scale, exp(xbeta_inter_i[l]) ) ,
                        weibull_lccdf(Vout_iU[l] | scale, exp(xbeta_inter_i[l]) ));

    } else{
      log_lik[N_o + N_m + l] = log_diff_exp( lognormal_lccdf(Vout_iL[l] | xbeta_inter_i[l], scale ) ,
                              lognormal_lccdf(Vout_iU[l] | xbeta_inter_i[l], scale ));

    }

    if(N_l > 0){
      if(L_ind_m[l] != -1){
        if(weib_ind == 1){
          log_lik[N_o + N_m + l] += -weibull_lccdf(Vout_l[L_ind_i[l]] | scale, exp(xbeta_inter_l[L_ind_i[l]]));
        } else{
          log_lik[N_o + N_m + l] += -lognormal_lccdf(Vout_l[L_ind_i[l]] | xbeta_inter_l[L_ind_i[l]], scale);          
        }
      }
    }
  }


}
