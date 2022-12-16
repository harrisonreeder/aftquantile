#script to test on the cluster
#read in packages
library(survival)
library(rstan)
# library(loo)
library(splines2)

options(mc.cores=1)
rstan_options(auto_write=TRUE)

# #Locally
# iter <- 1
# N <- 2000
# stanseed <- 1234

#on the cluster
#these are 'environment variables' defined in the sbatch command
iter <- commandArgs(trailingOnly=TRUE) #iteration number (basically, "j"), passed via bash script execute line as "${SLURM_ARRAY_TASK_ID}"
N <- eval(parse(text=Sys.getenv('N')))
stanseed <- eval(parse(text=Sys.getenv('stanseed')))

R <- 300 #total number of sims
J <- 150 #total number of jobs parallelized (should match array=1-J in the command line)
P <- R/J #number of simulations that THIS job will do.

scripts <- "scripts/"
output <- "output/"
stanpath <- "stan/"

iter <- as.numeric(iter)

print("number of cores available in this session is:")
print(parallel::detectCores())

source("sim_functions.R")
source("AFT_functions_spline.R")
source("BayesSurv_AFTtvstan.R")

## Compile/load models ##
##*********************##

begin <- Sys.time()
stan_AFT_tbp_mvn_tvcov_piece_compiled <- stan_model("stan/AFT_tbp_mvn_tvcov_piece.stan",verbose = TRUE)
end <- Sys.time()
print("tbp pw aft compilation time:")
print(end-begin)

begin <- Sys.time()
stan_AFT_tbp_mvn_tvcov_spline_compiled <- stan_model("stan/AFT_tbp_mvn_tvcov_spline.stan",verbose = TRUE)
end <- Sys.time()
print("tbp spline aft compilation time:")
print(end-begin)

##*******************************##
####Specify simulation settings####
##*******************************##

dgm_base <- "lognormal"
nP <- 2 #two continuous covariates (and one binary)

#generate data to reflect the data application
int_ln <- 3.2
scl_ln <- 0.55

#would be good to come up with values for these that reflect the data application
beta_base_temp <- c(-0.5,0.5,-0.2)

#for piecewise tv models
knots_piece <- c(0,7.5,15,22.5,30)
alpha_tv_temp_piece <- c(0,0.3,0.45,0.5) #these are "changes" from beta3 above

# beta_tv_temp_piece <- c(0.5,-0.4,0.4,0)
# beta_tv_temp_piece <- c(0.2,0.2,0,-0.2) #these are absolute numbers, not "differences"

#true lognormal baseline hazard
t_seq_temp <- seq(from=0, to=40, length=101)[-1]
true_basehaz <- h0_func(t=t_seq_temp, inter=int_ln, sc = scl_ln,
                        tbp = FALSE, baseline=dgm_base)
true_basesurv <- plnorm(q=t_seq_temp, meanlog=int_ln,sdlog=scl_ln,
                        lower.tail = FALSE,log.p=FALSE)

#sampling points for evaluating acceleration factors
#marginal model takes longer to compute AF so we choose a few
p_seq_marg_temp <- c(0.9,0.75,0.5,0.25,0.1)
#conditional model is quick to compute AF so we can do many
p_seq_cond_temp <- seq(from=0.01,to=0.99,by=0.01)

##*************************##
####begin simulation loop####
##*************************##

begin_time_1 <- Sys.time()
for(r in 1:P){
  begin_time_r <- Sys.time()

  #counter that represents which simulation out of R this is. Also used to give each sim a unique seed.
  iter_seed<-(as.numeric(iter)-1)*P+r
  print(paste0("r: ",r," of ",P,", overall sim iter: ", iter_seed," of ",R))

  ## generate covariate data ##
  ##*************************##

  set.seed(iter_seed)
  #generate three covariates, one of which is binary
  #binary covariate will have time-varying effect
  x_bin_temp <- sample(x = c(0,1),size = N,replace = TRUE)
  x_base_only_temp <- matrix(data=rnorm(N*nP), ncol = nP)
  x_base_temp <- cbind(x_base_only_temp, x_bin_temp)
  colnames(x_base_temp) <- c("x1","x2","x3")

  # dgm_type <- "piecewise"
  for(dgm_type in c("baseline",
                    "piecewise",NULL)){
    begin_time_dgm <- Sys.time()
    ## generate outcome data ##
    ##***********************##

    set.seed(iter_seed)
    temp_dat <- sim_AFTtv_better(x_base = x_base_temp, x_tv = x_bin_temp,
                                 beta_base_true = beta_base_temp,
                                 tv_type = dgm_type,
                                 beta_tv_true = alpha_tv_temp_piece,
                                 knots = if(dgm_type=="piecewise") knots_piece else NULL,
                                 dist = dgm_base, intercept=int_ln, scale = scl_ln,
                                 cens = c(15,40))
    Y_temp <- as.matrix(cbind(temp_dat$y, #Left "edge" of censoring
                              ifelse(temp_dat$delta==1,
                                     temp_dat$y,Inf), #Right "edge" of censoring
                              0))

    ## fit models and save results ##
    ##*****************************##

    for(model_base in c(#"lognormal",
                        "weibull",
                        NULL)){
      begin_time_base <- Sys.time()

      ## fit frequentist parametric model for initialization ##
      ##*****************************************************##

      freq_model_temp <- survreg(Surv(time = temp_dat$y,event = temp_dat$delta) ~ x_base_temp,
                                 dist = model_base)
      #grab the means and vcov for intercept and log(scale) from frequentist fit
      prior_mean_centering <- summary(freq_model_temp)$table[c("(Intercept)","Log(scale)"),1]
      freq_vcov_centering <- vcov(freq_model_temp)[c("(Intercept)","Log(scale)"),c("(Intercept)","Log(scale)")]
      freq_beta <- freq_model_temp$coefficients[-1]

      #flip sign of second parameter to make consistent with how we've previously been specifying it
      prior_mean_centering[2] <- -prior_mean_centering[2]
      #flip sign of covariance to be consistent with how we've previously been specifying scale
      #this is basically a "delta method" type adjustment for the slippery naming of log(shape) = -log(scale)
      freq_vcov_centering[1,2] <- -freq_vcov_centering[1,2]
      freq_vcov_centering[2,1] <- -freq_vcov_centering[2,1]
      prior_prec_centering <- MASS::ginv(freq_vcov_centering)/10

      init_func_temp <- function(){
        list(centering=prior_mean_centering,
             beta=freq_beta[1:nP], #excludes covariate with potential tv effect
             # beta_tv=as.array(rep(tail(freq_beta,n = 1),nP_tv_temp)), #start at "no quantile-varying effect"
             # beta_tv=as.array(c(tail(freq_beta,n = 1),numeric(nP_tv_temp-1))),
             w = rep(1/5,5), #hardcoding that J will be 5
             dirichlet_alpha_param=as.array(1))
      }

      ## fit bayesian model ##
      ##********************##

      # model_type <- "piecewise"
      for(model_type in c("invariant","piecewise","spline")){
        print(paste0("dgm_type: ",dgm_type,", dgm_base: ",dgm_base))
        print(paste0("model_type: ",model_type,", model_base: ",model_base))
        begin_time_eff <- Sys.time()
        print(begin_time_eff)

        #labels for saving objects at the end
        dgm_base_label <- switch(dgm_base, "lognormal"="ln", "weibull"="wb")
        dgm_type_label <- switch(dgm_type, "baseline"="inv", "piecewise"="pw")
        model_base_label <- switch(model_base, "lognormal"="ln", "weibull"="wb")
        model_type_label <- switch(model_type, "invariant"="inv","piecewise"="pw","spline"="sp")

        #bookkeeping variables for calling the right models
        effect_temp <- switch(model_type, "invariant"="piecewise","piecewise"="piecewise","spline"="spline")
        nP_tv_temp <- switch(model_type, "invariant"=1,"piecewise"=5,"spline"=4)
        knots_temp <- if(model_type=="piecewise") knots_piece else NULL

        begin_time_fit <- Sys.time()
        control_temp <- list(adapt_delta=0.95)
        test_fit <- BayesSurv_AFTtvstan(Y = Y_temp,
                                        Xmat = x_base_only_temp,
                                        Xtv = x_bin_temp,
                                        baseline = model_base, #we're going for a weibull model this time
                                        tbp = TRUE, #apply tbp
                                        J = 5, dirichlet_alpha_fixed = FALSE, dirichlet_alpha_data = 1,
                                        a_dirichlet_alpha = 1, b_dirichlet_alpha = 1,
                                        prior_centering = "mvn", #for tbp models, use centering prior
                                        prior_mean_centering = prior_mean_centering,
                                        prior_prec_centering = prior_prec_centering,
                                        prior_intercept = "none", m_intercept = 0, sd_intercept = 0,
                                        prior_scale = "gamma", a_scale=0.3, b_scale=0.05,
                                        prior_beta = "none", m_beta=0, sd_beta=10,
                                        prior_beta_tv = "none", m_beta_tv=0, sd_beta_tv=10,
                                        knots = knots_temp,
                                        tv_type = effect_temp, nP_tv = nP_tv_temp,
                                        n_chains=1, n_cores = 1,
                                        # n_sample = 400, n_warmup=400,
                                        n_sample = 2000, n_warmup=2000,
                                        # n_chains=4, n_sample = 5000, n_warmup=3000,
                                        init=init_func_temp,
                                        seed=stanseed, control = control_temp)
        end_time_fit <- Sys.time()

        #create the output list and begin filling
        out_list_temp <- list()
        out_list_temp[["summary"]] <- summary(test_fit$stan_fit,
                                              pars=test_fit$pars)$summary
        print(out_list_temp[["summary"]])

        #catch if there are no samples...
        if(is.null(out_list_temp[["summary"]])){
          print("initialization failed, keep going!")
          next
        }

        out_list_temp[["diagnostics"]] <-
          c(n_divergences=get_num_divergent(test_fit$stan_fit),
            n_max_treedepth=get_num_max_treedepth(test_fit$stan_fit),
            n_low_bfmi_chains=sum(get_bfmi(test_fit$stan_fit)<0.2))
        check_hmc_diagnostics(test_fit$stan_fit)

        #save the array of samples EXCEPT for the log-likelihoods, which take
        #up disproportionately too much space...
        out_list_temp[["samp_array"]] <-
          as.array(test_fit$stan_fit)[,,!(dimnames(as.array(test_fit$stan_fit))[[3]] %in% paste0("log_lik[",1:N,"]")),drop=FALSE]
        out_list_temp[["x_base"]] <- x_base_temp
        out_list_temp[["dat"]] <- temp_dat
        if(model_type != "invariant"){
          out_list_temp[["knots"]] <- test_fit$knots
        }

        ## Compute statistics for performance metrics ##
        ##********************************************##

        #1. loocv metric
        out_list_temp[["loo"]] <- rstan::loo(test_fit$stan_fit)
        print(out_list_temp[["loo"]])

        #2. hazard from 0 to 40
        #this, I think I'll estimate by like a sort of riemann sum I guess
        cond_basehaz_temp <- predict.basehaz(stan_fit = test_fit,
                              t_seq = t_seq_temp,verbose = FALSE)
        out_list_temp[["cond_basehaz"]] <- cond_basehaz_temp

        #2. survivor curve from 0 to 40
        cond_basesurv_temp <-
          predict.AFTtvstan(stan_fit = test_fit,t_seq = t_seq_temp,
                            type ="conditional",newXtv = 0, newXtv_time = 0,
                            thin = 1,verbose=FALSE)
        out_list_temp[["cond_basesurv"]] <- cond_basesurv_temp

        #3. conditional af at various quantiles
        begin_time_condaf <- Sys.time()
        cond_af_temp <-
          AF.AFTtvstan_tvcov(stan_fit = test_fit, p_seq = p_seq_cond_temp,
              type = "conditional", thin=1,verbose=TRUE)
        out_list_temp[["cond_af"]] <- cond_af_temp
        end_time_condaf <- Sys.time()
        print(paste(round(difftime(end_time_condaf, begin_time_condaf,units="mins"),3),"minutes: conditional AF"))

        #4. marginal af at various quantiles
        begin_time_margaf <- Sys.time()
        #turns out that this is the same across quantiles in the invariant setting
        #so this will hopefully speed things up just a tad
        if(model_type == "invariant"){
          marg_af_temp <- AF.AFTtvstan_tvcov(stan_fit = test_fit, p_seq = 0.5,
                                             oldXmat = x_base_only_temp,
                                             type = "marginal", thin=1,
                                             verbose=TRUE)[rep(1,length(p_seq_marg_temp)),]
        } else{
          marg_af_temp <- AF.AFTtvstan_tvcov(stan_fit = test_fit, p_seq = p_seq_marg_temp,
                                              oldXmat = x_base_only_temp,
                                              type = "marginal", thin=1,
                                              verbose=TRUE)
        }
        out_list_temp[["marg_af"]] <- marg_af_temp
        end_time_margaf <- Sys.time()
        print(paste("Marg AF took",
          round(difftime(end_time_margaf, begin_time_margaf,units="mins"),3),"minutes."))

        out_list_temp[["timing"]] <- list(fit=end_time_fit-begin_time_fit,
                                      cond_af=end_time_condaf-begin_time_condaf#,
                                      # marg_af=end_time_margaf-begin_time_margaf
                                      )

        saveRDS(out_list_temp, file = paste0(output,"outlist_sim_",iter_seed,
                                        "_n_",N,
                                        "_dgm_",dgm_base_label,"_",dgm_type_label,
                                        "_model_tbp",model_base_label,"_",model_type_label,
                                        "_stanseed_",stanseed,
                                        ".RDS"))
        end_time_eff <- Sys.time()
        print(paste(round(difftime(end_time_eff, begin_time_eff,units="mins"),3),"minutes: model effect"))
      } #end loop through model effect types
      end_time_base <- Sys.time()
      print(paste(round(difftime(end_time_base, begin_time_base,units="mins"),3),"minutes: model baseline"))
    } #end loop through model baseline types
    end_time_dgm <- Sys.time()
    print(paste(round(difftime(end_time_dgm, begin_time_dgm,units="mins"),3),"minutes: dgm"))
  } #end loop through dgm effect types
  end_time_r <- Sys.time()
  print(paste(round(difftime(end_time_r, begin_time_r,units="mins"),3),"minutes: whole iteration"))

  avg_runtime <- round(difftime(end_time_r, begin_time_1,units="mins")/r,3)
  print(paste(avg_runtime,"minutes: average iteration"))
  print(paste(round(avg_runtime*(P-r)/60,3)," hours: expected remaining runtime."))
} #end loop through iterations
end_time_1 <- Sys.time()
print(paste(round(difftime(end_time_1, begin_time_1,units="mins"),3),"minutes: job"))

