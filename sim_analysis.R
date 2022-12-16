#script to test on the cluster
#read in packages
library(survival)
library(rstan)
library(loo)
library(tidyverse)
#install.packages("rstpm2")

# N <- 2500
N <- 2000
R <- 300

output <- "[path]" #place to store simulation results
source("sim_functions.R")
source("AFT_functions_spline.R")
source("BayesSurv_AFTtvstan.R")

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

#true lognormal baseline hazard
t_seq_temp <- seq(from=0, to=40, length=101)[-1]
true_basehaz <- h0_func(t=t_seq_temp, inter=int_ln, sc = scl_ln,
                        tbp = FALSE, baseline=dgm_base)
true_basesurv <- plnorm(q=t_seq_temp, meanlog=int_ln,sdlog=scl_ln,
                        lower.tail = FALSE,log.p=FALSE)

#sampling points for evaluating acceleration factors
#marginal model takes longer to compute AF so we choose a few
p_seq_marg_temp <- c(0.9,0.75,0.5,0.25,0.10)
#conditional model is quick to compute AF so we can do many
#this is insane, but I found that if I don't round to 2 decimals it yields 0.1 != 0.1
p_seq_cond_temp <- round(seq(from=0.01,to=0.99,by=0.01),2)

##************************************##
####Compute true effects of interest####
##************************************##
#true conditional acceleration factors
true_condaf_list <- list()

true_condaf_list[["baseline"]] <- rep(exp(beta_base_temp[3]),length(p_seq_cond_temp))
true_condaf_list[["piecewise"]] <-
  AF.AFTtv_tvcov(int_temp = int_ln, scale_temp = scl_ln,
                 beta_temp = as.matrix(0),
                 tv_type = "piecewise", knots_temp = knots_piece,
                 beta_tv_temp = t(c(beta_base_temp[3],beta_base_temp[3] + alpha_tv_temp_piece)),
                 baseline = dgm_base,
                 tbp = FALSE, w_temp = rep(1/5,5),
                 p_seq = p_seq_cond_temp,
                 type = "conditional", #defaults are for comparison setting others to 0
                 thin=1,verbose=TRUE)[,1]

## Compute true marginal acceleration factors by numerical integration ##
##*********************************************************************##
#true marginal acceleration factors
true_margaf_list <- list()

surv_ln_pw_func <- function(x, t, int, scl, beta_base, z, beta_tv, knots){
  # browser()
  if(z==1){
    basis <- pw_cum_mat(y=t,
                        knots=knots[-length(knots)], #remove the Inf knot
                        intercept = TRUE) #note we include intercept here!
    V0_temp <- basis %*% exp(-beta_tv)
  } else {V0_temp <- t}
  eta <- as.vector(beta_base %*% x)
  plnorm(q=V0_temp, meanlog=int + eta, sdlog=scl,
         lower.tail = FALSE,log.p=FALSE)
}
margsurv_ln_pw_integrand_func <- function(x, t, int, scl, beta_base, z, beta_tv, knots){
  surv_ln_pw_func(x=x,t=t,int=int,scl=scl,beta_base=beta_base,
                  z=z,beta_tv=beta_tv,knots=knots) * prod(dnorm(x))
}
margsurv_ln_pw_func <- function(t, int, scl, beta_base, z, beta_tv, knots){
  tryCatch(
  cubature::cubintegrate(margsurv_ln_pw_integrand_func,
                         lower=rep(-Inf,2),upper=rep(Inf,2),fDim=length(t),method="pcubature",
                         t=t,int=int,scl=scl,beta_base=beta_base,z=z,beta_tv=beta_tv,knots=knots)$integral,
  error=function(e){return(rep(NA,length(t)))})
}
marginvsurv_ln_pw_func <- function(p, int, scl, beta_base, z, beta_tv, knots){
    rstpm2::vuniroot(f = function (t){p -
      margsurv_ln_pw_func(t=t,int=int,scl=scl,beta_base=beta_base,
                          z=z,beta_tv=beta_tv,knots=knots)},
      lower = rep(1e-10,length(p)), upper = rep(1000,length(p)))$root#,
}

true_margaf_list[["baseline"]] <-
  marginvsurv_ln_pw_func(p = p_seq_marg_temp,
                       int = int_ln, scl = scl_ln,
                       beta_base = t(beta_base_temp[1:2]), z = 1,
                       knots = c(0,Inf), beta_tv=beta_base_temp[3])/
  marginvsurv_ln_pw_func(p = p_seq_marg_temp,
                         int = int_ln, scl = scl_ln,
                         beta_base = t(beta_base_temp[1:2]), z = 0,
                         knots = c(0,Inf), beta_tv=beta_base_temp[3])

true_margaf_list[["piecewise"]] <-
  marginvsurv_ln_pw_func(p = p_seq_marg_temp,
                         int = int_ln, scl = scl_ln,
                         beta_base = t(beta_base_temp[1:2]), z = 1,
                         knots = c(knots_piece,Inf),
                         beta_tv=c(beta_base_temp[3],beta_base_temp[3] + alpha_tv_temp_piece))/
  marginvsurv_ln_pw_func(p = p_seq_marg_temp,
                         int = int_ln, scl = scl_ln,
                         beta_base = t(beta_base_temp[1:2]), z = 0,
                         knots = c(knots_piece,Inf),
                         beta_tv=c(beta_base_temp[3],beta_base_temp[3] + alpha_tv_temp_piece))

# ## Check true marginal acceleration factors by monte carlo ##
# ##*********************************************************##
# set.seed(1234)
# temp_oldmat <- matrix(data=rnorm(1e6*nP), ncol = nP)
# true_marg_af_list[["baseline_mc"]] <-
#   AF.AFTtv_tvcov(int_temp = int_ln, scale_temp = scl_ln,
#                  beta_temp = t(beta_base_temp[1:2]), oldXmat = temp_oldmat,
#                  tv_type = "piecewise", knots_temp = c(0,Inf),
#                  beta_tv_temp = as.matrix(beta_base_temp[3]),
#                  baseline = dgm_base,
#                  tbp = FALSE, w_temp = rep(1/5,5),
#                  p_seq = p_seq_marg_temp,
#                  type = "marginal", #defaults are for comparison setting others to 0
#                  thin=1,verbose=TRUE)[,1]
#
# true_marg_af_list[["piecewise_mc"]] <-
#   AF.AFTtv_tvcov(int_temp = int_ln, scale_temp = scl_ln,
#                  beta_temp = t(beta_base_temp[1:2]), oldXmat = temp_oldmat,
#                  tv_type = "piecewise", knots_temp = knots_piece,
#                  beta_tv_temp = t(c(beta_base_temp[3],beta_base_temp[3] + alpha_tv_temp_piece)),
#                  baseline = dgm_base,
#                  tbp = FALSE, w_temp = rep(1/5,5),
#                  p_seq = p_seq_marg_temp,
#                  type = "marginal", #defaults are for comparison setting others to 0
#                  thin=1,verbose=TRUE)[,1]
# rbind(p_seq_marg_temp, do.call(rbind,true_marg_af_list))
# p_seq_marg_temp 0.9000000 0.7500000 0.5000000 0.2500000 0.1000000
# baseline        0.8187320 0.8187307 0.8187300 0.8187308 0.8187308
# baseline_mc     0.8187320 0.8187307 0.8187300 0.8187308 0.8187308
# piecewise       0.8187322 0.8187308 0.8912553 1.0875332 1.1975723
# piecewise_mc    0.8187323 0.8187308 0.8911339 1.0874597 1.1975424


##**************************************##
#### COMPUTE FINAL SIMULATION RESULTS ####
##**************************************##

#storage data frames for the simulation results
frame_af <- NULL
frame_base <- NULL
condaf_fullframe <- NULL
temp_condaf_list <- temp_margaf_list <- temp_basehaz_list <- list()

for(dgm_type in c("baseline",
                  "piecewise",NULL)){

  # model_type <- "piecewise"
  for(model_type in c("invariant","piecewise","spline")){

    # model_base <- "lognormal"
    for(model_base in c("lognormal",
                        "weibull",
                        NULL)){

      temp_condaf_list[[paste(dgm_type,model_type,model_base,sep = "_")]] <-
        matrix(ncol=length(p_seq_cond_temp),nrow=R)
      temp_margaf_list[[paste(dgm_type,model_type,model_base,sep = "_")]] <-
        matrix(ncol=length(p_seq_marg_temp),nrow=R)
      temp_basehaz_list[[paste(dgm_type,model_type,model_base,sep = "_")]] <-
        matrix(ncol=length(t_seq_temp),nrow=R)

      temp_margaf_mse_mat <- temp_margaf_bias_mat <- temp_margaf_sd_mat <- temp_margaf_cp_mat <-
        matrix(ncol=length(true_margaf_list[[dgm_type]]),nrow=R)
      temp_condaf_mse_mat <- temp_condaf_bias_mat <- temp_condaf_sd_mat <- temp_condaf_cp_mat <-
        matrix(ncol=length(true_condaf_list[[dgm_type]]),nrow=R)
      temp_isq_basehaz <- temp_looic <- numeric(R)

      print(paste0("dgm_type: ",dgm_type,", dgm_base: ",dgm_base))
      print(paste0("model_type: ",model_type,", model_base: ",model_base))
      for(r in 1:R){ #just 2 for now
        if(r %% 10 ==0){print(paste0("iteration: ", r, ", ",Sys.time()))}

        #labels for saving objects at the end
        dgm_base_label <- switch(dgm_base, "lognormal"="ln", "weibull"="wb")
        dgm_type_label <- switch(dgm_type, "baseline"="inv", "piecewise"="pw")
        model_base_label <- switch(model_base, "lognormal"="ln", "weibull"="wb")
        model_type_label <- switch(model_type, "invariant"="inv","piecewise"="pw","spline"="sp")

        if(!file.exists(paste0(output,"outlist_sim_",r,
                               "_n_",N,
                               "_dgm_",dgm_base_label,"_",dgm_type_label,
                               "_model_",model_base_label,"_",model_type_label,
                               "_2022-10-27.RDS"))){
          print("skipping file")
          next
        }

        out_list <- readRDS(file = paste0(output,"outlist_sim_",r,
                             "_n_",N,
                             "_dgm_",dgm_base_label,"_",dgm_type_label,
                             "_model_",model_base_label,"_",model_type_label,
                             "_2022-10-27.RDS"))

        #results on conditional acceleration factors
        temp_condaf_bias_mat[r,] <- out_list$cond_af[,"median"] - true_condaf_list[[dgm_type]]
        temp_condaf_sd_mat[r,] <- out_list$cond_af[,"sd"]
        temp_condaf_mse_mat[r,] <- temp_condaf_bias_mat[r,]^2 + temp_condaf_sd_mat[r,]^2
        temp_condaf_cp_mat[r,] <-
          as.numeric(out_list$cond_af[,"lb"] <= true_condaf_list[[dgm_type]] &
                       out_list$cond_af[,"ub"] >= true_condaf_list[[dgm_type]])

        #results on marginalized acceleration factors
        temp_margaf_bias_mat[r,] <- out_list$marg_af[,"median"] - true_margaf_list[[dgm_type]]
        temp_margaf_sd_mat[r,] <- out_list$marg_af[,"sd"]
        temp_margaf_mse_mat[r,] <- temp_margaf_bias_mat[r,]^2 + temp_margaf_sd_mat[r,]^2
        temp_margaf_cp_mat[r,] <-
          as.numeric(out_list$marg_af[,"lb"] <= true_margaf_list[[dgm_type]] &
                      out_list$marg_af[,"ub"] >= true_margaf_list[[dgm_type]])

        #results on baseline hazard and loo
        #do a sort of "riemann" integral using the 80 sample points for the baseline hazard
        #by definition diff(t_seq_temp) widths are all the same except for numerical issues
        #so just use the first one as the "rectangle width"
        temp_isq_basehaz[r] <-
          sum((out_list$cond_basehaz[,"median"] - true_basehaz)^2) * diff(t_seq_temp)[1]
        temp_looic[r] <- out_list$loo$estimates[3,1]


        temp_condaf_fullframe <- as.data.frame(t(out_list$cond_af[,"median"]))
        colnames(temp_condaf_fullframe) <- paste0("p",p_seq_cond_temp)
        temp_condaf_fullframe$r <- r
        temp_condaf_fullframe$dgm_type <- dgm_type
        temp_condaf_fullframe$model_type <- model_type
        temp_condaf_fullframe$model_base <- model_base

        condaf_fullframe <- rbind(condaf_fullframe,temp_condaf_fullframe)

        #fill in detailed matrices
        temp_condaf_list[[paste(dgm_type,model_type,model_base,sep = "_")]][r,] <-
          out_list$cond_af[,"median"]
        temp_margaf_list[[paste(dgm_type,model_type,model_base,sep = "_")]][r,] <-
          out_list$marg_af[,"median"]
        temp_basehaz_list[[paste(dgm_type,model_type,model_base,sep = "_")]][r,] <-
          out_list$cond_basehaz[,"median"]

      } #end loop through iterations

      condaf_length <- length(p_seq_cond_temp)
      frame_condaf <- data.frame(dgm_base=dgm_base,dgm_type=dgm_type,
                 model_type=model_type,model_base=model_base,
                 effect="cond",p=p_seq_cond_temp,
                 bias=apply(temp_condaf_bias_mat,2,median,na.rm=TRUE),
                 stddev=apply(temp_condaf_sd_mat,2,median,na.rm=TRUE),
                 emp_se=apply(temp_condaf_bias_mat,2,sd,na.rm=TRUE), #sd of bias is same as sd of original estimates
                 mse=apply(temp_condaf_mse_mat,2,sd,na.rm=TRUE),
                 cp=apply(temp_condaf_cp_mat,2,mean,na.rm=TRUE))
      frame_margaf <- data.frame(dgm_base=dgm_base,dgm_type=dgm_type,
                 model_type=model_type,model_base=model_base,
                 effect="marg",p=p_seq_marg_temp,
                 bias=apply(temp_margaf_bias_mat,2,median,na.rm=TRUE),
                 stddev=apply(temp_margaf_sd_mat,2,median,na.rm=TRUE),
                 emp_se=apply(temp_margaf_bias_mat,2,sd,na.rm=TRUE), #sd of bias is same as sd of original estimates
                 mse=apply(temp_margaf_mse_mat,2,sd,na.rm=TRUE),
                 cp=apply(temp_margaf_cp_mat,2,mean,na.rm=TRUE))
      frame_af <- rbind(frame_af,frame_condaf,frame_margaf)

      frame_base <- rbind(frame_base,
        data.frame(dgm_base=dgm_base,dgm_type=dgm_type,
          model_type=model_type,model_base=model_base,
          looic=median(temp_looic,na.rm=TRUE),
          isqbasehaz=median(temp_isq_basehaz,na.rm=TRUE)))
    } #end loop through model baseline types

  } #end loop through model effect types

} #end loop through dgm effect types

#table results
tab_condaf <- frame_af %>% select(-emp_se,-mse) %>%
  filter(effect=="cond", p %in% c(0.75,0.5,0.25)) %>%
  arrange(dgm_base,dgm_type,model_type,model_base,-p) %>% #make it go left to right...
  pivot_wider(names_from=p,values_from=c(bias,stddev, cp)) %>%
  mutate(id=paste0("dgm: ",dgm_base,"_",dgm_type,", model: ", model_type,"_",model_base)) %>%
  relocate(id) %>%
  select(-effect,-dgm_base,-dgm_type,-model_type,-model_base)
tab_condaf
#add spacer columns between the sets of results
print(xtable::xtable(
  cbind(tab_condaf[,1:4],NA,tab_condaf[,5:7],NA,tab_condaf[,8:10]),
  digits=3),include.rownames = FALSE)

tab_margaf <- frame_af %>% select(-emp_se,-mse) %>%
  filter(effect=="marg", p %in% c(0.75,0.5,0.25)) %>%
  pivot_wider(names_from=p,values_from=c(bias,stddev,
                                         cp)) %>%
  mutate(id=paste0("dgm: ",dgm_base,"_",dgm_type,", model: ", model_type,"_",model_base)) %>%
  relocate(id) %>%
  select(-effect,-dgm_base,-dgm_type,-model_type,-model_base)
tab_margaf
#add spacer columns between the sets of results
print(xtable::xtable(
  cbind(tab_margaf[,1:4],NA,tab_margaf[,5:7],NA,tab_margaf[,8:10]),
                           digits=3),include.rownames = FALSE)

tab_base <- frame_base %>%
  mutate(id=paste0("dgm: ",dgm_base,"_",dgm_type,", model: ", model_type,"_",model_base)) %>%
  relocate(id) %>% select(-dgm_base,-dgm_type,-model_type,-model_base)
tab_base %>% arrange(looic)
tab_base %>% filter(str_detect(id,"lognormal_baseline")) %>% arrange(looic)
tab_base %>% filter(str_detect(id,"lognormal_piecewise")) %>% arrange(looic)
tab_base %>% filter(str_detect(id,"lognormal_baseline")) %>% arrange(isqbasehaz)
tab_base %>% filter(str_detect(id,"lognormal_piecewise")) %>% arrange(isqbasehaz)

print(xtable::xtable(tab_base,digits=c(0,0,1,6)),include.rownames=FALSE)



#### TBP ####

stanseed <- 1234

#storage data frames for the simulation results
frame_af <- NULL
frame_base <- NULL

status_frame <- NULL
fittime_vec <- NULL

for(dgm_type in c("baseline",
                  "piecewise",NULL)){
  
  # model_type <- "piecewise"
  for(model_type in c("invariant","piecewise","spline")){
    
    model_base <- "tbpweibull"
    # for(model_base in c("lognormal",
    #                     "weibull",
    #                     NULL)){
    
    temp_margaf_mse_mat <- temp_margaf_bias_mat <- temp_margaf_sd_mat <- temp_margaf_cp_mat <-
      matrix(ncol=length(true_margaf_list[[dgm_type]]),nrow=R)
    temp_condaf_mse_mat <- temp_condaf_bias_mat <- temp_condaf_sd_mat <- temp_condaf_cp_mat <-
      matrix(ncol=length(true_condaf_list[[dgm_type]]),nrow=R)
    temp_isq_basehaz <- temp_looic <- rep(NA,R)
    
    print(paste0("dgm_type: ",dgm_type,", dgm_base: ",dgm_base))
    print(paste0("model_type: ",model_type,", model_base: ",model_base))
    for(r in 1:R){ #just 2 for now
      if(r %% 10 ==0){print(paste0("iteration: ", r, ", ",Sys.time()))}
      
      #labels for saving objects at the end
      dgm_base_label <- switch(dgm_base, "lognormal"="ln", "weibull"="tbpwb")
      # model_base_label <- switch(model_base, "lognormal"="ln", "weibull"="tbpwb")
      model_base_label <- "tbpwb"
      dgm_type_label <- switch(dgm_type, "baseline"="inv", "piecewise"="pw")
      model_type_label <- switch(model_type, "invariant"="inv","piecewise"="pw","spline"="sp")
      
      if(!file.exists(paste0(output,"outlist_sim_",r,
                             "_n_",N,
                             "_dgm_",dgm_base_label,"_",dgm_type_label,
                             "_model_",model_base_label,"_",model_type_label,
                             "_stanseed_",stanseed,
                             "_2022-10-27.RDS"))){
        print("skipping file")
        
        status_frame <- rbind(status_frame,
                              c(r=r, dgm_base=dgm_base_label, dgm_type=dgm_type_label,
                                model_base=model_base_label, model_type=model_type_label,
                                have_fit=0,actually_ran=NA, time=NA
                              ))
        
        next
      }
      
      out_list <- readRDS(file = paste0(output,"outlist_sim_",r,
                                        "_n_",N,
                                        "_dgm_",dgm_base_label,"_",dgm_type_label,
                                        "_model_",model_base_label,"_",model_type_label,
                                        "_stanseed_",stanseed,
                                        "_2022-10-27.RDS"))
      
      actually_worked <- as.numeric(min(out_list$summary[,"n_eff"],na.rm = TRUE) > 100)
      
      status_frame <- rbind(status_frame,
                            c(r=r, dgm_base=dgm_base_label, dgm_type=dgm_type_label,
                              model_base=model_base_label, model_type=model_type_label,
                              have_fit=1,actually_ran=actually_worked,
                              time=as.numeric(out_list$timing$fit, unit = "mins")
                            ))
      
      if(!actually_worked){
        print("ran but didn't actually work...")
        next
      }
      
      #results on conditional acceleration factors
      temp_condaf_bias_mat[r,] <- out_list$cond_af[,"median"] - true_condaf_list[[dgm_type]]
      temp_condaf_sd_mat[r,] <- out_list$cond_af[,"sd"]
      temp_condaf_mse_mat[r,] <- temp_condaf_bias_mat[r,]^2 + temp_condaf_sd_mat[r,]^2
      temp_condaf_cp_mat[r,] <-
        as.numeric(out_list$cond_af[,"lb"] <= true_condaf_list[[dgm_type]] &
                     out_list$cond_af[,"ub"] >= true_condaf_list[[dgm_type]])
      
      marg_af <- readRDS(file = paste0(output_marg,"margafmat_sim_",r,
                                       "_n_",N,
                                       "_dgm_",dgm_base_label,"_",dgm_type_label,
                                       "_model_",model_base_label,"_",model_type_label,
                                       "_stanseed_",stanseed,
                                       "_2022-10-27.RDS"))
      
      #results on marginalized acceleration factors
      temp_margaf_bias_mat[r,] <- marg_af[,"median"] - true_margaf_list[[dgm_type]]
      temp_margaf_sd_mat[r,] <- marg_af[,"sd"]
      temp_margaf_mse_mat[r,] <- temp_margaf_bias_mat[r,]^2 + temp_margaf_sd_mat[r,]^2
      temp_margaf_cp_mat[r,] <-
        as.numeric(marg_af[,"lb"] <= true_margaf_list[[dgm_type]] &
                     marg_af[,"ub"] >= true_margaf_list[[dgm_type]])
      
      #results on baseline hazard and loo
      #do a sort of "riemann" integral using the 80 sample points for the baseline hazard
      #by definition diff(t_seq_temp) widths are all the same except for numerical issues
      #so just use the first one as the "rectangle width"
      temp_isq_basehaz[r] <-
        sum((out_list$cond_basehaz[,"median"] - true_basehaz)^2) * diff(t_seq_temp)[1]
      temp_looic[r] <- out_list$loo$estimates[3,1]
      
    } #end loop through iterations
    
    condaf_length <- length(p_seq_cond_temp)
    frame_condaf <- data.frame(dgm_base=dgm_base,dgm_type=dgm_type,
                               model_type=model_type,model_base=model_base,
                               effect="cond",p=p_seq_cond_temp,
                               bias=apply(temp_condaf_bias_mat,2,median,na.rm=TRUE),
                               stddev=apply(temp_condaf_sd_mat,2,median,na.rm=TRUE),
                               emp_se=apply(temp_condaf_bias_mat,2,sd,na.rm=TRUE), #sd of bias is same as sd of original estimates
                               mse=apply(temp_condaf_mse_mat,2,sd,na.rm=TRUE),
                               cp=apply(temp_condaf_cp_mat,2,mean,na.rm=TRUE))
    frame_margaf <- data.frame(dgm_base=dgm_base,dgm_type=dgm_type,
                               model_type=model_type,model_base=model_base,
                               effect="marg",p=p_seq_marg_temp,
                               bias=apply(temp_margaf_bias_mat,2,median,na.rm=TRUE),
                               stddev=apply(temp_margaf_sd_mat,2,median,na.rm=TRUE),
                               emp_se=apply(temp_margaf_bias_mat,2,sd,na.rm=TRUE), #sd of bias is same as sd of original estimates
                               mse=apply(temp_margaf_mse_mat,2,sd,na.rm=TRUE),
                               cp=apply(temp_margaf_cp_mat,2,mean,na.rm=TRUE))
    frame_af <- rbind(frame_af,
                      frame_condaf,
                      frame_margaf
    )
    
    frame_base <- rbind(frame_base,
                        data.frame(dgm_base=dgm_base,dgm_type=dgm_type,
                                   model_type=model_type,model_base=model_base,
                                   looic=median(temp_looic,na.rm=TRUE),
                                   isqbasehaz=median(temp_isq_basehaz,na.rm=TRUE)))
    # } #end loop through model baseline types
    
  } #end loop through model effect types
  
} #end loop through dgm effect types

status_frame <- as.data.frame(status_frame)
status_frame$r <- as.numeric(status_frame$r)
status_frame$have_fit <- as.numeric(status_frame$have_fit)
status_frame$actually_ran <- as.numeric(status_frame$actually_ran)
status_frame$time <- as.numeric(status_frame$time)
sum(status_frame$have_fit)
sum(status_frame$actually_ran,na.rm=TRUE)
status_frame %>% group_by(dgm_base,dgm_type,model_type) %>%
  summarize(have_fit=sum(have_fit),actually_ran=sum(actually_ran,na.rm=TRUE)) %>%
  mutate(prop = actually_ran / have_fit)

#which models had some kind of problem? about 30% did..
status_frame %>% filter(!have_fit | actually_ran == 0) %>% View
#how many models actually "ran" at a bare minimum? 1233, so there are another 567 that didn't...
status_frame %>% summarize(sum(actually_ran,na.rm=TRUE))
#look by dataset, and see how many of the three models fit correctly
status_frame %>% group_by(r, dgm_base,dgm_type) %>%
  summarize(have_fit=sum(have_fit),actually_ran=sum(actually_ran,na.rm=TRUE)) %>%
  mutate(prop = actually_ran / have_fit) %>% ggplot(data=., aes(x=actually_ran)) + geom_histogram()
#which datasets had one or two run? Just try rerunning those from a different seed
status_frame %>% filter(!have_fit) %>% View
status_frame %>% group_by(r, dgm_base,dgm_type) %>%
  summarize(have_fit=sum(have_fit),actually_ran=sum(actually_ran,na.rm=TRUE)) %>%
  filter(have_fit==0)
df <- status_frame %>% group_by(r, dgm_base,dgm_type) %>%
  summarize(have_fit=sum(have_fit),actually_ran=sum(actually_ran,na.rm=TRUE)) %>%
  filter(have_fit==3 & actually_ran==1)
df$r[df$dgm_type=="inv"]
df$r[df$dgm_type=="pw"]
status_frame %>% filter(dgm_type=="inv" & r %in% df$r[df$dgm_type=="inv"])
#table results
tab_condaf <- frame_af %>% select(-emp_se,-mse) %>%
  filter(effect=="cond", p %in% c(0.75,0.5,0.25)) %>%
  arrange(dgm_base,dgm_type,model_type,model_base,-p) %>% #make it go left to right...
  pivot_wider(names_from=p,values_from=c(bias,stddev, cp)) %>%
  mutate(id=paste0("dgm: ",dgm_base,"_",dgm_type,", model: ", model_type,"_",model_base)) %>%
  relocate(id) %>%
  select(-effect,-dgm_base,-dgm_type,-model_type,-model_base)
tab_condaf
#add spacer columns between the sets of results
print(xtable::xtable(
  cbind(tab_condaf[,1:4],NA,tab_condaf[,5:7],NA,tab_condaf[,8:10]),
  digits=3),include.rownames = FALSE)

tab_margaf <- frame_af %>% select(-emp_se,-mse) %>%
  filter(effect=="marg", p %in% c(0.75,0.5,0.25)) %>%
  pivot_wider(names_from=p,values_from=c(bias,stddev,
                                         cp)) %>%
  mutate(id=paste0("dgm: ",dgm_base,"_",dgm_type,", model: ", model_type,"_",model_base)) %>%
  relocate(id) %>%
  select(-effect,-dgm_base,-dgm_type,-model_type,-model_base)
tab_margaf
#add spacer columns between the sets of results
print(xtable::xtable(
  cbind(tab_margaf[,1:4],NA,tab_margaf[,5:7],NA,tab_margaf[,8:10]),
  digits=3),include.rownames = FALSE)

tab_base <- frame_base %>%
  mutate(id=paste0("dgm: ",dgm_base,"_",dgm_type,", model: ", model_type,"_",model_base)) %>%
  relocate(id) %>% select(-dgm_base,-dgm_type,-model_type,-model_base)
tab_base %>% arrange(looic)
tab_base %>% filter(str_detect(id,"lognormal_baseline")) %>% arrange(looic)
tab_base %>% filter(str_detect(id,"lognormal_piecewise")) %>% arrange(looic)
tab_base %>% filter(str_detect(id,"lognormal_baseline")) %>% arrange(isqbasehaz)
tab_base %>% filter(str_detect(id,"lognormal_piecewise")) %>% arrange(isqbasehaz)

print(xtable::xtable(tab_base,digits=c(0,0,1,6)),include.rownames=FALSE)
