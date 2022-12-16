####Present analysis for paper####

####READ THINGS IN AND SET UP DATA
#read in packages
sapply(c("survival","flexsurv","rstpm2","SemiCompRisks","SemiCompRisksPen","rstan","loo",
         "tidyverse","readxl","lubridate","tidylog","survminer","visdat","AFTTV"),
       function(x){if(!suppressWarnings(require(package=x,character.only = T)))
       {install.packages(pkgs=x);require(package=x,character.only = T)}})

figurepath <- "[path]"
ROSMAPtemp <- "[path]"
three_color_qual <- RColorBrewer::brewer.pal(n=3,name="Set2")
four_color_qual <- rev(RColorBrewer::brewer.pal(n=4,name="Dark2"))

source("00_import_rosmap_data.R")
source("BayesSurv_AFTtvstan.R")

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

#### Baseline Table ####
totcount <- length(rosmap_baseline_sub$studyid)
catcounts <- tapply(rosmap_baseline_sub$studyid, rosmap_baseline_sub$outcome_cat,length)
catvars <- c("race_fctwhite","msex","marital_bl_fctmarried","education_15plus","apoe4_any")
basetable <- NULL
for(var in catvars){
  cutoffval <- 1 #all vars are binary 1/0 in our case
  varcounts <- tapply(rosmap_baseline_sub[,var][[1]], 
                      rosmap_baseline_sub$outcome_cat,
                      function(x){sum(x==cutoffval,na.rm=T)})
  varna <- tapply(rosmap_baseline_sub[,var][[1]],
                  rosmap_baseline_sub$outcome_cat,
                  function(x){sum(is.na(x))})
  testset <- rosmap_baseline_sub[,var][[1]]
  temprow <- c( 
    paste0(varcounts," (",round(varcounts/(catcounts-varna),3)*100,"%)"),
    NULL
  )
  namesvec <- c(rownames(basetable),var)
  
  basetable <-rbind(basetable,temprow)
  rownames(basetable) <- namesvec
}

colnames(basetable) <- levels(rosmap_baseline_sub$outcome_cat)

totcolumn <- NULL
for(var in catvars){
  varcountstot <- sum(rosmap_baseline_sub[,var][[1]]==cutoffval,na.rm=T)
  totcolumn <- c(totcolumn,
                 paste0(varcountstot," (",round(varcountstot/(totcount),3)*100,"%)"))
}

basetable <- cbind(total=totcolumn,basetable)
basetable <- rbind(total=paste0(c(totcount,catcounts), " (100%)"),basetable)
basetable
basetable <- basetable[,c("total","neither","ad_diagnosis_only","death_only","both")]

recode_varnames <- function(x){
  recode(x,
         "total"="Total",
         "race_fctwhite"="White Race/Ethnicity",
         "msex"="Male Sex",
         "marital_bl_fctmarried"="Married at Study Entry",
         "education_15plus"="15+ Years of Education",
         "apoe4_any"="APOE4 Genetic Variant")
}
rownames(basetable) <- recode_varnames(rownames(basetable))
colnames(basetable) <- recode(colnames(basetable),"total"="Total","neither"= "Censored prior to AD/dementia or death","ad_diagnosis_only"="AD/dementia and censored prior to death","death_only"="Death without AD/dementia","both"="AD/dementia diagnosis and death")
xtable::xtable(basetable)

#### Cox Modeling ####

# AD
form_AD_cox_lhs_string <- "Surv(time=time_bl,time2=time_ad_dx,event=ad_diagnosis,type = 'counting')"
form_rhs_string <- paste0(Xmat_names,collapse = " + ")
form_AD_cox <- as.formula(paste(form_AD_cox_lhs_string,"~",form_rhs_string))
form_AD_cox_tv <- update(form_AD_cox, . ~ . + tt(apoe4_any))

cox_AD_fit <- coxph(form_AD_cox, data=rosmap_baseline_sub)
summary(cox_AD_fit)
cox_AD_fit_zph <- cox.zph(cox_AD_fit,global = TRUE)
cox_AD_fit_zph
plot(cox_AD_fit_zph)

S0_cox_AD_m <- survexp(~apoe4_any, data=rosmap_baseline_sub, ratetable=cox_AD_fit)

#Death
rosmap_baseline_sub_long <- 
  tmerge(data1=rosmap_baseline_sub,
         data2=rosmap_baseline_sub, #set data as both "data1", "data2"
         id=studyid, #"subject" is our identifier variable for each person
         death = event(time_last, death), #define "death" by time2 and event2 variables
         ad_tv_ind=tdc(time_ad_dx)) #create a variable "x1" for surgery status from time1 variable
rosmap_baseline_sub_long$tstart[rosmap_baseline_sub_long$tstart==0] <- rosmap_baseline_sub_long$time_bl[rosmap_baseline_sub_long$tstart==0]

form_death_cox_lhs_string <- "Surv(time=tstart,time2=tstop,event=death,type = 'counting')"
form_death_cox <- as.formula(paste(form_death_cox_lhs_string, "~", form_rhs_string, "+ ad_tv_ind"))
form_death_cox_tv <- update(form_death_cox, . ~ . + tt(ad_tv_ind))

cox_death_fit <- coxph(form_death_cox, data=rosmap_baseline_sub_long)
summary(cox_death_fit)

cox_death_fit_zph <- cox.zph(cox_death_fit,global = TRUE)
cox_death_fit_zph
plot(cox_death_fit_zph)

S0_cox_death_m <- survexp(~ad_tv_ind, data=rosmap_baseline_sub_long,
                          ratetable=cox_death_fit)

#### Read In Stan Models ####
fit_list <- list()
fit_list[["apoe"]] <- list()
fit_list[["deathtv"]] <- list()
for(outcome in c("apoe","deathtv",
                 NULL)){
  for(eff in c("invariant","piece","spline",NULL)){
    print(paste(outcome,base,eff))
    temp_fit <- readRDS(file = paste0(ROSMAPtemp,"ROSMAP_",base,"_",
                                      outcome,"_",eff,".RDS"))
    fit_list[[outcome]][[paste(base,eff,sep="_")]] <- temp_fit
  }
}

for(outcome in c("apoe","deathtv",
                 NULL)){
  for(eff in c("invariant","piece","spline",NULL)){
    for(base in c("wb","ln")){
      print(paste(outcome,base,eff))
      print(fit_list[[outcome]][[paste(base,eff,sep="_")]]$knots)
    }
    print(paste(outcome,base,eff))
    print(fit_list[[outcome]][[paste("tbp",eff,sep="_")]]$knots)
  }
}

#### Model Evaluation ####
library("loo")
loo_list <- list()
loo_list[["deathtv"]] <- loo_list[["apoe"]] <- list()
loo_tab_list <- list()
#AD/dementia Onset
for(outcome in c(#"apoe",
                 "deathtv",
                 NULL)){
  for(eff in c("invariant","piece","spline",NULL)){
    for(base in c(#"wb","ln",
                  "tbp",
                  NULL)){
      print(paste(outcome,base,eff))
      temp_fit <- fit_list[[outcome]][[paste(base,eff,sep="_")]]
      log_lik_temp <- extract_log_lik(temp_fit$stan_fit, merge_chains = FALSE)
      rel_eff_temp <- relative_eff(exp(log_lik_temp), cores = 2)
      loo_list[[outcome]][[paste(base,eff,sep="_")]] <- loo(log_lik_temp, r_eff = rel_eff_temp, cores = 2)
      print(loo_list[[outcome]][[paste(base,eff,sep="_")]])
    }
  }
  loo_tab_list[[paste0(outcome,"_long")]] <- as.data.frame(loo_compare(x=loo_list[[outcome]])[,"looic",drop=FALSE]) %>% 
    rownames_to_column(var="model") %>%
    separate(model,sep="_",into=c("baseline","effect")) %>%
    arrange(factor(baseline,levels = c("ln","wb","tbp"))) %>%
    mutate(baseline = recode_factor(baseline,"tbp"="TBP (Weibull Centered)", "wb"="Weibull","ln"="log-Normal"),
           effect = recode_factor(effect,"invariant"="Constant", "piece"="Piecewise Linear","spline"="Restricted Cubic Spline")) %>% arrange(looic)
  loo_tab_list[[paste0(outcome,"_wide")]] <- loo_tab_list[[paste0(outcome,"_long")]] %>%
    pivot_wider(names_from=baseline,values_from=looic) %>%
    relocate(effect, `log-Normal`, `Weibull`, `TBP (Weibull Centered)`) %>% arrange(effect)
}
#copy and paste the tables from this output
invisible(lapply(loo_tab_list[str_detect(names(loo_tab_list),"wide")],function(x){
  # print(x)
  print(xtable::xtable(x,digits=1),include.rownames = FALSE)
  }))

#### Covariate Estimates ####
spline_list <- piece_list <- cov_list <- list()
spline_list[["deathtv"]] <- spline_list[["apoe"]] <- 
  piece_list[["deathtv"]] <- piece_list[["apoe"]] <- 
  cov_list[["deathtv"]] <- cov_list[["apoe"]] <- list()

for(outcome in c("apoe","deathtv",NULL)){
  for(eff in c(#"invariant",
               #"piece",
               "spline",
               NULL)){
    for(base in c("wb","ln","tbp")){
      print(paste(outcome,base,eff))
      temp_fit <- fit_list[[outcome]][[paste(base,eff,sep="_")]]
      print(check_hmc_diagnostics(temp_fit$stan_fit))
      summary_temp <- round(summary(temp_fit$stan_fit,pars=c("beta","beta_tv"))$summary,2)
      
      #now, to transform the parameters to the "beta, alpha" parameterization
      #which, in the sampler, they are not to reduce parameter correlation
      if(eff == "piece"){
        betatv_samp_mat <- apply(X = as.array(temp_fit$stan_fit)[,,paste0("beta_tv[",1:temp_fit$nP_tv,"]"),drop=FALSE],
                                 MARGIN=3,FUN = as.vector)
        piece_list[[outcome]][[paste(base,eff,sep="_")]] <- 
          t(apply(X = cbind(betatv_samp_mat[,1], betatv_samp_mat[,2:5] - betatv_samp_mat[,1]), 
                MARGIN=2, FUN=quantile, probs=c(0.5,0.025,0.975)))
      } else if(eff == "spline"){
        knots_temp <- temp_fit$knots
        if(outcome == "apoe"){
          Q_temp <- get_basis_tv(x=seq(from=0.01, to=40, by=0.01), tv_type="rp",
                                 knots=knots_temp[-c(1,length(knots_temp))], #remove the Inf knot
                                 deriv = FALSE, intercept=TRUE, #include the intercept
                                 flexsurv_compatible = FALSE)
          B_temp <- cbind(1, #this is the intercept
                          get_basis_tv(x=seq(from=0.01, to=40, by=0.01), tv_type="rp",
                                       knots=knots_temp[-c(1,length(knots_temp))], #remove the Inf knot
                                       deriv = FALSE, intercept=FALSE, #exclude the intercept
                                       flexsurv_compatible = FALSE))
        } else{
          Q_temp <- get_basis_tv(x=seq(from=0.01, to=20, by=0.01), tv_type="rp",
                                 knots=knots_temp[-c(1,length(knots_temp))], #remove the Inf knot
                                 deriv = FALSE, intercept=TRUE, #include the intercept
                                 flexsurv_compatible = FALSE)
          B_temp <- cbind(1, #this is the intercept
                          get_basis_tv(x=seq(from=0.01, to=20, by=0.01), tv_type="rp",
                                       knots=knots_temp[-c(1,length(knots_temp))], #remove the Inf knot
                                       deriv = FALSE, intercept=FALSE, #exclude the intercept
                                       flexsurv_compatible = FALSE))
        }
        betatv_samp_mat <- apply(X = as.array(temp_fit$stan_fit)[,,paste0("beta_tv[",1:temp_fit$nP_tv,"]"),drop=FALSE],
                                 MARGIN=3,FUN = as.vector)
        #transformation matrix
        transform_mat <- solve(t(B_temp) %*% B_temp) %*% t(B_temp) %*% Q_temp
        spline_list[[outcome]][[paste(base,eff,sep="_")]] <- 
          t(apply(transform_mat %*% t(betatv_samp_mat), MARGIN=1, FUN=quantile, probs=c(0.5,0.025,0.975)))
        
        y_temp <- Q_temp %*% betatv_samp_mat[1,]
        print(summary(lm( y_temp ~ B_temp-1)))
        
      }
      cov_list[[outcome]][[paste(base,eff,sep="_")]] <- paste0(summary_temp[,6]," (",summary_temp[,4],", ",summary_temp[,8],")")
    }
  }
}
cov_list[["apoe"]][[paste("blank","blank",sep="_")]] <- NA
summary_temp <- round(cbind(cox_AD_fit$coefficients,confint(cox_AD_fit)),2)
cov_list[["apoe"]][[paste("ph","invariant",sep="_")]] <- paste0(summary_temp[,1]," (",summary_temp[,2],", ",summary_temp[,3],")")
cov_list[["deathtv"]][[paste("blank","blank",sep="_")]] <- NA
summary_temp <- round(cbind(cox_death_fit$coefficients,confint(cox_death_fit)),2)
cov_list[["deathtv"]][[paste("ph","invariant",sep="_")]] <- paste0(summary_temp[,1]," (",summary_temp[,2],", ",summary_temp[,3],")")

cov_mat_list <- list()
for(outcome in c("apoe","deathtv")){
  cov_mat <- do.call(what=cbind,args=cov_list[[outcome]]) #not all vectors are same length that's ok!
  Xmat_names_temp <- if(outcome=="apoe") Xmat_baselineonly_names else Xmat_names
  rownames(cov_mat) <- c(Xmat_names_temp,
                         paste0("tv",1:(NROW(cov_mat)-length(Xmat_names_temp))))
  
  cov_outtable_temp <- as.data.frame(cov_mat) %>% rownames_to_column(var="variable") %>% 
    pivot_longer(cols=-variable,names_to="model",values_to="value") %>%
    separate(model,into = c("baseline","effect"),sep="_") %>% 
    arrange(factor(variable,levels = c(Xmat_names,paste0("tv",1:9))),
            factor(baseline,levels = c("blank","ph","ln","wb","tbp"))) %>%
    mutate(value = ifelse(effect=="spline" & str_detect(string = variable,pattern = "tv[5-9]"), "",value),
           value = ifelse(effect=="invariant" & str_detect(string = variable,pattern = "tv[2-9]"), "",value),
           baseline = recode(baseline,"ph"="Cox Proportional Hazards",
                             "tbp"="TBP (Weibull Centered)", "wb"="Weibull","ln"="log-Normal"),
           effect = recode(effect,"invariant"="Constant", "piece"="Piecewise Linear","spline"="Restricted Cubic Spline"))
  
  cov_outtable_wide <- cov_outtable_temp %>%
    pivot_wider(names_from="effect",values_from="value") %>%
    mutate(Constant = ifelse(str_detect(string = variable,pattern = "tv[2-9]"),"",Constant),
           `Restricted Cubic Spline` = ifelse(str_detect(string = variable,
                                                         pattern = paste0("tv[",fit_list[[outcome]][["wb_spline"]][["nP_tv"]] + 1,"-9]")),
                                              "",`Restricted Cubic Spline`),
           variable = ifelse(baseline=="blank",variable,as.character(baseline)),
           variable = recode_varnames(variable)) %>% select(-baseline,-blank)  
  
  #remove extraneous values for constant and restricted cubic spline, where necessary.
  cov_outtable_long <- cov_outtable_temp %>%
    pivot_wider(names_from="baseline",values_from="value") %>%
    mutate(variable = ifelse(effect=="blank",variable,as.character(effect)),
           variable = recode_varnames(variable)) %>% select(-effect,-blank)
  cov_mat_list[[paste(outcome,"wide",sep="_")]] <- cov_outtable_wide
  cov_mat_list[[paste(outcome,"long",sep="_")]] <- cov_outtable_long
}
#copy and paste the tables from this output
invisible(lapply(cov_mat_list[str_detect(names(cov_mat_list),"long")],function(x){
  print(xtable::xtable(x),include.rownames = FALSE)
}))

#to build final regression tables, replace time-varying coefficients with the
#corresponding values from below
View(cov_mat_list$apoe_long)
View(cov_mat_list$deathtv_long)
round(piece_list$apoe$ln_piece,2)
round(piece_list$apoe$wb_piece,2)
round(piece_list$apoe$tbp_piece,2)
round(piece_list$deathtv$ln_piece,2)
round(piece_list$deathtv$wb_piece,2)
round(piece_list$deathtv$tbp_piece,2)

round(spline_list$apoe$ln_spline,2)
round(spline_list$apoe$wb_spline,2)
round(spline_list$apoe$tbp_spline,2)
round(spline_list$deathtv$ln_spline,2)
round(spline_list$deathtv$wb_spline,2)
round(spline_list$deathtv$tbp_spline,2)

#### AD: Survivor Curves ####
t_seq_temp <- seq(from=0, to=40, length=80)[-1]

##generate plot data
exposed_survs_marg_AD_list <- baseline_survs_marg_AD_list <- list()
for(base in c("wb","ln", "tbp",
              NULL)){
  thin_amt <- 1
  print(base)
  baseline_survs_marg_AD_list <- c(baseline_survs_marg_AD_list,
                                   lapply(X = fit_list[["apoe"]][paste0(base,"_",c("invariant","piece","spline"))],
                                                function(x){
                                                  predict.AFTtvstan(stan_fit = x, t_seq = t_seq_temp,
                                                    oldXmat = Xmat_baselineonly,type ="marginal",
                                                    newXtv = 0, newXtv_time = 0,
                                                    thin = thin_amt,verbose=TRUE)
                                                }))
  saveRDS(baseline_survs_marg_AD_list,file = paste0(ROSMAPtemp,"ROSMAP_AD_baseline_survs_marg.RDS"))  
  exposed_survs_marg_AD_list <- c(exposed_survs_marg_AD_list,
                                  lapply(X = fit_list[["apoe"]][paste0(base,"_",c("invariant","piece","spline"))],
                                    function(x){
                                      predict.AFTtvstan(stan_fit = x, t_seq = t_seq_temp,
                                        oldXmat = Xmat_baselineonly, type ="marginal",
                                        newXtv = 1, newXtv_time = 0,
                                        thin = thin_amt,verbose=TRUE)
                                }))
  saveRDS(exposed_survs_marg_AD_list,file = paste0(ROSMAPtemp,"ROSMAP_AD_exposed_survs_marg.RDS"))  
}

#### AD: Acceleration Factors ####
p_seq_temp <- seq(from=0.01,to=0.99,by=0.02)

#generate and save marginal AFs
begin <- Sys.time()
AF_marg_AD_list <- list()
for(base in c("wb","ln", "tbp",
              NULL)){
  thin_amt <- 1
  print(base)
  AF_marg_AD_list <- c(AF_marg_AD_list,
                       lapply(X = fit_list[["apoe"]][paste0(base,"_",c("invariant","piece","spline"))],
                          function(x){
                            AF.AFTtvstan_tvcov(stan_fit = x, p_seq = p_seq_temp,
                               oldXmat = Xmat_baselineonly, type = "marginal", thin=thin_amt,verbose=TRUE)
                          }))
  saveRDS(AF_marg_AD_list,file = paste0(ROSMAPtemp,"ROSMAP_AD_AF_marg.RDS"))  
}
end <- Sys.time()
end-begin

#### Death: Survivor Curves ####
t_seq_temp <- seq(from=0, to=40, length=80)[-1]

#generate and save marginal curves for t = 5 and t = 20
exposed_survs_marg_death_list <- baseline_survs_marg_death_list <- list()
for(base in c("wb","ln", "tbp",
              NULL)){
  thin_amt <- 1
  print(base)
  baseline_survs_marg_death_list <- c(baseline_survs_marg_death_list,
                                   lapply(X = fit_list[["deathtv"]][paste0(base,"_",c("invariant","piece","spline"))],
                                          function(x){
                                            predict.AFTtvstan(stan_fit = x, t_seq = t_seq_temp,
                                                              oldXmat = Xmat, type ="marginal",
                                                              newXtv = 0, newXtv_time = 0,
                                                              thin = thin_amt,verbose=TRUE)
                                          }))
  saveRDS(baseline_survs_marg_death_list,file = paste0(ROSMAPtemp,"ROSMAP_deathtv_baseline_survs_marg.RDS"))
  
  for(Xtv_time_temp in c(5,20)){
    print(Xtv_time_temp)
    temp_list <- lapply(X = fit_list[["deathtv"]][paste0(base,"_",c("invariant","piece","spline"))],
                        function(x){
                          predict.AFTtvstan(stan_fit = x, t_seq = t_seq_temp,
                                            oldXmat = Xmat, type ="marginal",
                                            newXtv = 1, newXtv_time = Xtv_time_temp,
                                            thin = thin_amt,verbose=TRUE)
                        })
    names(temp_list) <- paste0(names(temp_list),"_time",Xtv_time_temp)
    exposed_survs_marg_death_list <- c(exposed_survs_marg_death_list,temp_list)
    saveRDS(exposed_survs_marg_death_list,file = paste0(ROSMAPtemp,"ROSMAP_deathtv_exposed_survs_marg.RDS"))
  }
}

#### Death: Acceleration Factors ####
p_seq_temp <- c(seq(from=0.01,to=0.99,by=0.02),0.995,0.999,0.9995)

AF_marg_death_list <- list()
for(base in c("wb","ln", "tbp",
              NULL)){
  thin_amt <- 1
  print(base)
  for(Xtv_time_temp in c(5,20)){
    print(Xtv_time_temp)
    temp_list <- lapply(X = fit_list[["deathtv"]][paste0(base,"_",c("invariant","piece","spline"))],
                        function(x){
                          AF.AFTtvstan_tvcov(stan_fit = x, p_seq = p_seq_temp, t_Xtv_seq = Xtv_time_temp, 
                                             oldXmat = Xmat, type = "marginal",thin=thin_amt,verbose=TRUE)
                        })
    names(temp_list) <- paste0(names(temp_list),"_time",Xtv_time_temp)
    AF_marg_death_list <- c(AF_marg_death_list,temp_list)
    saveRDS(AF_marg_death_list,file = paste0(ROSMAPtemp,"ROSMAP_deathtv_AF_marg.RDS"))  
  }
}

#### Death: Acceleration Surface ####
t_seq_temp <- seq(from=0, to=45, length=90)[-1]
p_seq_temp <- seq(from=0.01,to=0.99,by=0.01)

#invariant only
AF_marg_surface_long_list <- list()
for(base in c("wb","ln", "tbp",
              NULL)){
  thin_amt <- 10
  AF_marg_surface_long_list[[paste0(base,"_",c("invariant"))]] <- AF.AFTtvstan_tvcov(stan_fit = fit_list[["deathtv"]][[paste0(base,"_",c("invariant"))]],
                                  p_seq = p_seq_temp,type = "conditional",
                                  # oldXmat = Xmat,
                                  t_Xtv_seq = t_seq_temp,
                                  # p_Xtv_seq = p_seq_temp,
                                  thin=thin_amt,verbose = TRUE,long_out = TRUE)
  saveRDS(AF_marg_surface_long_list,file = paste0(ROSMAPtemp,"ROSMAP_deathtv_AFsurface_marg.RDS"))
}





####PLOT RESULTS####
setEPS(horizontal = FALSE, onefile = FALSE, paper = "special")

#use this to add little tags to the corners of figures
#https://www.r-bloggers.com/2017/03/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

t_seq_surface <- seq(from=0, to=45, length=90)[-1]
t_seq_temp <- seq(from=0, to=40, length=80)[-1]
p_seq_temp <- c(seq(from=0.01,to=0.99,by=0.01),0.995,0.999,0.9995)
n_samp <- 30000

#### AD: Survivor Curves ####
for(base in c("wb","ln", "tbp",
              NULL)){
  baseline_survs_marg <- do.call(what = cbind, args = baseline_survs_marg_AD_list[paste0(base,"_",c("invariant","piece","spline"))])
  exposed_survs_marg <- do.call(what = cbind, args = exposed_survs_marg_AD_list[paste0(base,"_",c("invariant","piece","spline"))])
  #marginal curves
  # cairo_pdf(file = paste0(figurepath,"ROSMAP_apoe_marg_surv_",base,".pdf"), width = 5,height=6)
  postscript(file = paste0(figurepath,"ROSMAP_apoe_marg_surv_",base,".eps"), width = 5,height=6)
  plot(S0_cox_AD_m,col=c("grey80",four_color_qual[1]),lwd=2, ylim=c(0,1),
       xlab="Time to AD/Dementia without Death, Years from Age 65",
       ylab="Survivor Function",conf.int = FALSE)
  matplot(x = t_seq_temp, y = exposed_survs_marg[,c(1,4,7)], type = "l",
          col=four_color_qual[2:4],lty=1,lwd = 2, add=TRUE)
  matplot(x = t_seq_temp,y = baseline_survs_marg[,c(1,4,7)],type = "l",
          col=c("grey60","grey40","grey20"),lty=1,lwd=2,add=TRUE)
  legend("topright", title="No APOE-e4",fill=c("grey80","grey60","grey40","grey20"),cex=0.75,
         legend=c("PH, Constant", "AFT, Constant","AFT, Piecewise Linear", "AFT, Spline"))
  legend("bottomleft", title="APOE-e4",fill=four_color_qual,cex=0.75,
         legend=c("PH, Constant", "AFT, Constant","AFT, Piecewise Linear", "AFT, Spline"))
  fig_label(text = "(a)",region = "figure",pos = "topleft")
  dev.off()
}

#### AD: Acceleration Factors ####
for(base in c("wb","ln", "tbp",
              NULL)){
  AF_marg <- do.call(what = cbind,args = AF_marg_AD_list[paste0(base,"_",c("invariant","piece","spline"))])
  baseline_survs_marg <- do.call(what = cbind, args = baseline_survs_marg_AD_list[paste0(base,"_",c("invariant","piece","spline"))])
  # cairo_pdf(file = paste0(figurepath,"ROSMAP_apoe_marg_af_",base,".pdf"), width = 5,height=6)
  postscript(file = paste0(figurepath,"ROSMAP_apoe_marg_af_",base,".eps"), width = 5,height=6)
  plot(x=p_seq_temp, y=rep(1,length(p_seq_temp)), type="l",
       col="white",xlim=c(1,0),ylim = c(0.2,1.2),
       xlab = "Quantile (p)", ylab="Acceleration Factor") #placeholder plot to set up frame
  rect(xleft=min(baseline_survs_marg),xright=0,ybottom=0,ytop=10,col = "grey96",border = NA)
  lines(x=p_seq_temp, rep(1,length(p_seq_temp)),type="l", lty=3,lwd=2,col="grey20")
  matplot(x = p_seq_temp,y = AF_marg,type = "l",lty=c(1,2,2),lwd=c(2,1,1),
          col=rep(four_color_qual[2:4],each=3),add=TRUE)
  legend("bottomright",
         legend=c("Constant","Piecewise Linear","Spline"),
         fill=four_color_qual[2:4],cex=0.8, bg="white")
  fig_label(text = "(b)",region = "figure",pos = "topleft")
  dev.off()
}

#### Death: Survivor Curves ####
for(base in c("wb","ln", "tbp",
              NULL)){
  #plot with all three flexible effects in the same plot
  baseline_survs_marg_death <-
    do.call(what = cbind,args = baseline_survs_marg_death_list[paste0(base,"_",c("invariant","piece","spline"))])
  # cairo_pdf(file = paste0(figurepath,"ROSMAP_deathtv_marg_surv_",base,".pdf"), width = 4,height=5.8)
  postscript(file = paste0(figurepath,"ROSMAP_deathtv_marg_surv_",base,".eps"), width = 4,height=5.8)
  par(mfrow=c(2,1))
  par(mar = c(4,4.1,0.5,0.1)) #bottom, left, top, right
  # par(mar = c(5.1, 4.1, 4.1, 2.1)) #default
  for(Xtv_time_temp in c(5,20)){
    exposed_survs_marg_death <-
      do.call(what = cbind,args = exposed_survs_marg_death_list[paste0(base,"_",c("invariant","piece","spline"),"_time",Xtv_time_temp)])
    matplot(x = t_seq_temp,y = baseline_survs_marg_death[,c(1,4,7)],type = "l",
            col=c("grey60","grey40","grey20"),lty=1,lwd=2,ylim=c(0,1),
            main="",ylab="Survivor Function",xlab=if(Xtv_time_temp == 20) "Time to Death, Years from Age 65" else "")
    matplot(x = t_seq_temp,y = exposed_survs_marg_death[,c(1,4,7)],type = "l",
            col=four_color_qual[2:4],lty=1,lwd=2, ylim=c(0,1), add=TRUE)
    legend("bottomleft", title="AD/Dementia Onset",legend=c("Constant","Piecewise Linear", "Spline"),
           fill=four_color_qual[2:4],cex=0.5)
    legend("topright", title="No AD/Dementia Onset",legend=c("Constant","Piecewise Linear", "Spline"),
           fill=c("grey60","grey40","grey20"),cex=0.5)
    fig_label(text = if(Xtv_time_temp==5) "(a)" else "(c)", region = "figure", pos = "topleft")
  }
  par(mfrow=c(1,1))
  dev.off()
  
  #plot with just constant effect
  baseline_survs_marg_death <- baseline_survs_marg_death_list[[paste0(base,"_",c("invariant"))]]
  # cairo_pdf(file = paste0(figurepath,"ROSMAP_deathtv_marg_surv_",base,"_invariant.pdf"), width = 4,height=5.8)
  postscript(file = paste0(figurepath,"ROSMAP_deathtv_marg_surv_",base,"_invariant.eps"), width = 4,height=5.8)
  par(mfrow=c(2,1))
  par(mar = c(4,4.1,0.5,0.1)) #bottom, left, top, right
  # par(mar = c(5.1, 4.1, 4.1, 2.1)) #default
  for(Xtv_time_temp in c(5,20)){
    exposed_survs_marg_death <- exposed_survs_marg_death_list[[paste0(base,"_",c("invariant"),"_time",Xtv_time_temp)]]
    plot(x = t_seq_temp,y = baseline_survs_marg_death[,1],type = "l", ylim=c(0,1),
         col=c("grey60"),lty=1,lwd=2,main="",ylab="Survivor Function",
         xlab=if(Xtv_time_temp == 20) "Time to Death, Years from Age 65" else "")
    lines(x = t_seq_temp,y = exposed_survs_marg_death[,1],
          type = "l",col=four_color_qual[2],lty=1,lwd=2)
    legend("bottomleft",
           legend=c("no AD Onset","AD Onset"),
           # legend=c("no AD/Dementia Onset","AD/Dementia Onset"),
           fill=c("grey60",four_color_qual[2]),cex=0.75)
  }
  par(mfrow=c(1,1))
  dev.off()
}

#### Death: Acceleration Factors ####
for(base in c("wb","ln", "tbp",
              NULL)){
  baseline_survs_marg_death <-
    do.call(what = cbind,args = baseline_survs_marg_death_list[paste0(base,"_",c("invariant","piece","spline"))])
  #plot with all three flexible effects in the same plot
  # cairo_pdf(file = paste0(figurepath,"ROSMAP_deathtv_marg_af_",base,".pdf"), width = 4,height=5.8)
  postscript(file = paste0(figurepath,"ROSMAP_deathtv_marg_af_",base,".eps"), width = 4,height=5.8)
  par(mfrow=c(2,1))
  par(mar = c(4,4.1,0.5,0.1)) #bottom, left, top, right
  # par(mar = c(5.1, 4.1, 4.1, 2.1)) #default
  for(Xtv_time_temp in c(5,20)){
    AF_marg_death <-
      do.call(what = cbind,args = AF_marg_death_list[paste0(base,"_",c("invariant","piece","spline"),"_time",Xtv_time_temp)])
    plot(x=p_seq_temp, y=rep(1,length(p_seq_temp)), type="l",
         col="white",xlim=c(1,0),ylim = c(0.2,1.2),
         xlab = if(Xtv_time_temp == 20) "Survival Quantile (p)" else "", 
         ylab="Acceleration Factor") #placeholder plot to set up frame
    rect(xleft=min(baseline_survs_marg_death),xright=0,ybottom=0,ytop=10,col = "grey96",border = NA)
    lines(x=p_seq_temp, rep(1,length(p_seq_temp)),type="l", lty=3,lwd=2,col="grey20")
    matplot(x = p_seq_temp,y = AF_marg_death,type = "l",lty=c(1,2,2),lwd=c(2,1,1),
            col=rep(four_color_qual[2:4],each=3),add=TRUE)
    legend("topright", legend=c("Constant","Piecewise Linear","Spline"),
           fill=four_color_qual[2:4],cex=0.5, bg="white")
    fig_label(text = if(Xtv_time_temp==5) "(b)" else "(d)", region = "figure", pos = "topleft")
  }
  par(mfrow=c(1,1))
  dev.off()
  #plot with just constant effect
  # cairo_pdf(file = paste0(figurepath,"ROSMAP_deathtv_marg_af_",base,"_invariant.pdf"), width = 4,height=5.8)
  postscript(file = paste0(figurepath,"ROSMAP_deathtv_marg_af_",base,"_invariant.eps"), width = 4,height=5.8)
  par(mfrow=c(2,1))
  par(mar = c(4,4.1,0.5,0.1)) #bottom, left, top, right
  # par(mar = c(5.1, 4.1, 4.1, 2.1)) #default
  for(Xtv_time_temp in c(5,20)){
    AF_marg_death <- AF_marg_death_list[[paste0(base,"_",c("invariant"),"_time",Xtv_time_temp)]]
    plot(x=p_seq_temp, y=rep(1,length(p_seq_temp)), type="l",
         col="white",xlim=c(1,0),ylim = c(0.5,1.2),
         xlab = if(Xtv_time_temp == 20) "Survival Quantile (p)" else "", 
         ylab="Acceleration Factor") #placeholder plot to set up frame
    rect(xleft=min(baseline_survs_marg_death),xright=0,ybottom=0,ytop=10,col = "grey96",border = NA)
    lines(x=p_seq_temp, rep(1,length(p_seq_temp)),type="l", lty=3,lwd=2,col="grey60")
    matplot(x = p_seq_temp,y = AF_marg_death,type = "l",lty=c(1,2,2),lwd=c(2,1,1),
            col=rep(four_color_qual[2:4],each=3),add=TRUE)
    legend("topright",
           # legend=c("no AD/Dementia Onset","AD/Dementia Onset Onset"),
           legend=c("no AD Onset","AD Onset"),
           fill=c("grey60",four_color_qual[2]),cex=0.5)
  }
  par(mfrow=c(1,1))
  dev.off()
}

#### Death: Acceleration Surface ####
breaks_plot <- seq(0.2,1.25,by = 0.05)
#invariant only
for(base in c("wb",
              "ln",
              "tbp",
              NULL)){
  for(eff in c("invariant",
               "piece",
               "spline",
               NULL)){
    label_temp <- switch(eff, "invariant"="(a)","piece"="(b)","spline"="(c)")
    
    AF_marg_surface_long <- AF_marg_surface_long_list[[paste0(base,"_",eff)]]
    if(eff=="piece" & base=="tbp"){
      #there's a single very small numerical value above 1 that is disrupting the 
      #coloring for the tbp piecewise model
      AF_marg_surface_long$AF[AF_marg_surface_long$AF>1] <- 1
    }
    
    if(eff =="spline"){
      n_blue_bins <- max(AF_marg_surface_long$AF) %/% .05 - 19
      n_red_bins <- 20 - min(AF_marg_surface_long$AF) %/% .05
      contour_colors <-
        rev(c(if(n_blue_bins > 0) rev(colorRampPalette(RColorBrewer::brewer.pal(6, "Blues"))(n_red_bins)[1:n_blue_bins]),
              colorRampPalette(RColorBrewer::brewer.pal(8, "Reds"))(n_red_bins)))
      #fix the coloring of the figures
      #specifically, we want to:
      #replace all 1's before the time of AD onset with NAs, but not after
      #without throwing out the values very early on when there are never any 1's in the row
      AF_marg_surface_long <- AF_marg_surface_long %>% group_by(Sinvx) %>%
        mutate(flag = ifelse(AF>1.000001,1,NA)) %>%
        fill(flag, .direction="up") %>% add_tally(wt=flag) %>%
        mutate(AF = ifelse(flag==1 | (n==0 & Sinvx < 10),AF, NA))
      
    } else {
      #how many steps of 5 get from the minimum value to 1?
      n_red_bins <- 20 - min(AF_marg_surface_long$AF) %/% .05
      contour_colors <-
        rev(colorRampPalette(RColorBrewer::brewer.pal(8, "Reds"))(n_red_bins))
      #https://www.datanovia.com/en/blog/easy-way-to-expand-color-palettes-in-r/
      # n_contour_bins <- 11
      # contour_colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(8, "Reds"))(n_contour_bins))
      AF_marg_surface_long$AF[AF_marg_surface_long$AF==1] <- NA      
    }
    
    width_temp <- 4; height_temp <- 6
    cairo_pdf(file = paste0(figurepath,"ROSMAP_deathtv_marg_afsurface_",base,"_",eff,"_qt.pdf"),
              width = width_temp,height=height_temp)
    # postscript(file = paste0(figurepath,"ROSMAP_deathtv_marg_afsurface_",base,"_",eff,"_qt.eps"),
    #         width = width_temp,height=height_temp)
    print(ggplot(data = AF_marg_surface_long,
                 mapping = aes(x = p,y = Sinvx,z = AF)) + theme_classic() +
            geom_contour_filled(breaks = breaks_plot) +
            # geom_segment(aes(x=0,xend=1,y=5,yend=5),linetype=3,color="grey20") +
            # geom_segment(aes(x=0,xend=1,y=20,yend=20),linetype=3,color="grey20") +
            xlim(1,0) +
            scale_fill_manual(values=contour_colors) +
            ylab("Years since 65 at AD Onset") + labs(tag=label_temp) + 
            xlab("Survival Quantile (p)") + #ggtitle(label = base,subtitle = eff) +
            guides(fill=guide_legend(title="AF",reverse=TRUE)))
    dev.off()
  }
}


##### TBP Prior Illustration ####

#define simple survivor function
S0_wb_tbp <- function(x,inter,alpha,w){
  J <- length(w)
  S0_wb <- pweibull(q=x,scale = exp(inter),shape = alpha,lower.tail = FALSE)
  S0_tbp <- 0
  for(j in 1:J){
    S0_tbp <- S0_tbp + w[j]*pbeta(q = S0_wb,shape1 = j,shape2 = J-j+1,ncp = 0)
  }
  S0_tbp
}

#verify that this is the same as the one above!
S0_wb_tbp2 <- function(x,inter,alpha,w){
  J <- length(w)
  S0_wb <- pweibull(q=x,scale = exp(inter),shape = alpha,lower.tail = FALSE)
  S0_tbp <- splines2::bernsteinPoly(x = S0_wb,
                                    degree = J-1,
                                    intercept = TRUE,
                                    Boundary.knots = c(0,1),
                                    integral = TRUE) %*% w
  return(S0_tbp * J)
}


#show the curves of the basis
pdf(file = paste0(figurepath,"tbp_basis_J5.pdf"),width = 4.5,height=3)
J <- 5
curve(pbeta(x,1,J),from = 0, to = 1,xlab=expression(p),ylab = expression(G(p)))
for(j in 2:J){
  curve(pbeta(x,j,J-j+1),add = TRUE)
}
dev.off()

#show the curves of the basis
z_list <- list(1:J,J:1,c(1,2,3,2,1),1 - c(1,4,5,2,1))
w_list <- lapply(z_list,function(x){exp(x)/sum(exp(x))})
# pdf(file = paste0(figurepath,"tbp_basis_ex.pdf"),width = 4,height=4)
curve(S0_wb_tbp(x,inter=0,alpha=1,w=rep(1/J,J)),from = 0,to = 4,lwd=2,
      xlab=expression(t),ylab = expression(S[0](t)))
for(i in 1:4){
  w_temp <- w_list[[i]]
  # curve(S0_wb_tbp(x,inter=0,alpha=1,w=w_temp),
  #       from = 0,to = 4,add = TRUE,col=four_color_qual[i])
  curve(S0_wb_tbp2(x,inter=0,alpha=1,w=w_temp),
        from = 0,to = 4,add = TRUE,col=four_color_qual[i],lty=2,lwd=2)
}
curve(S0_wb_tbp(x,inter=0,alpha=1,w=rep(1/J,J)),from = 0,to = 4,lwd=2,add=TRUE)
legend(x="topright",legend = 1:4,fill = four_color_qual)
dev.off()
print(xtable::xtable(do.call(what = cbind,args = w_list)))


test_vec <- rexp(n=1e4,rate = 1)
microbenchmark::microbenchmark(S0_wb_tbp(x = test_vec,inter=0,alpha=1,w=w_temp),
                               S0_wb_tbp2(x = test_vec,inter=0,alpha=1,w=w_temp))




set.seed(1234)
curve(S0_wb_tbp(x,inter=0,alpha=1,w=rep(1/J,J)),from = 0,to = 4)
for(i in 1:20){
  z <- rnorm(J)
  w_temp <- exp(z)/sum(exp(z))
  curve(S0_wb_tbp(x,inter=0,alpha=1,w=w_temp),
        from = 0,to = 4,add = TRUE,col="grey70")
}
curve(S0_wb_tbp(x,inter=0,alpha=1,w=rep(1/J,J)),from = 0,to = 4,add=TRUE)


#### Example Curves for Introductory Section ####
source("AFT_functions.R")

#setting up for some plotting
t_seq <- seq(from=0,to=10,by=0.01)

t_star <- 1
beta_temp <- 0.5

beta_temp_multi1 <- c(0.5,0.5,-0.6)
knots_temp1 <- c(0,0,0.5,2.5)

beta_temp_multi2 <- c(0.5,0,0)
knots_temp2 <- c(0,1.5,10,15)

beta_temp_multi3 <- c(-0.4,0,0)
knots_temp3 <- c(0,0,0.5,3)
S0_plot <- function(x){pexp(q=x,rate=0.3,lower.tail = FALSE)}
S0_inv_plot <- function(x){qexp(p=x,rate=0.3,lower.tail = FALSE)}

leg <- function(x="topright",plot=TRUE,cex=0.8){
  legend(x = x,title = expression(beta(t)),
         fill=c("grey","black","red","blue"),
         cex = cex,
         legend = c("0","1", "log(t+1)","piecewise"))
}

#survivor curve
cairo_pdf(paste0(figurepath,"S_ex.pdf"),width=5,height=6)
plot(t_seq,S0_plot(V(x = t_seq, beta = beta_temp, type = "baseline")), 
     type = "l",col="grey20",lwd=2,#lty=3,
     ylab="Survivor Function",
     # ylab=expression("Survivor Function,"~S(t~"|"~X)),
     # ylab=expression(S(t~"|"~X)~"="~S[0](V(t~"|"~X))), 
     xlab = "Time (t)", 
     ylim=c(0,1),yaxp=c(0,1,4))
lines(t_seq,S0_plot(V(x = t_seq, beta = beta_temp, type = "constant")), 
      type = "l", col=four_color_qual[2],lwd=2)
lines(t_seq,S0_plot(V(x = t_seq, beta = beta_temp_multi1, type = "piecewise",
                 knots=knots_temp1)), type = "l", col=four_color_qual[3],lwd=2)
lines(t_seq,S0_plot(V(x = t_seq, beta = beta_temp, type = "logplusone")), 
      type = "l", col=four_color_qual[4],lwd=2)
# leg(cex=1.2)
legend(x = "topright",title = expression(V(t)),
       fill=c("grey20",four_color_qual[2:4]),
       cex = 1,
       legend = c(expression(t),
                  expression(t %*% exp(-0.5*X)),
                  expression(t %*% exp((-0.5+0.6*I(t>2.5))*X)),
                  expression((t^{1-0.5*X} - 1) /(1-0.5*X))))
dev.off()

# now, let's look at ratio of inverse survivor function to inverse baseline survivor function
# this is a sort of "acceleration factor over the percentiles" visual
#here, because we know S0 is exponential, then inverse is -log(y)
step <- 0.01
p_seq <- seq(from=step,to=1-step,by=step)
S0p_seq <- S0_inv_plot(p_seq)

cairo_pdf(paste0(figurepath,"AF_ex.pdf"),width=5,height=6)
plot(p_seq,Vinv(x = S0p_seq, beta = beta_temp, type = "baseline")/S0p_seq,
     type = "l", xlim=c(1,0), ylim = c(0,3), xaxp=c(0,1,4),
     ylab="Acceleration Factor", 
     xlab = "Survival Quantile (p)", col="grey20",
     lwd=2,lty=3)
lines(p_seq,Vinv(x = S0p_seq, beta = beta_temp, type = "constant")/S0p_seq, 
      type = "l", col=four_color_qual[2],lwd=2)
lines(p_seq,Vx_inv(t_obj = S0p_seq, beta_base = rep(0,length(p_seq)), 
                   beta_tv = beta_temp_multi1,knots = knots_temp1,
                   tv_type = "piecewise")/S0p_seq, 
      type = "l", col=four_color_qual[3], lwd=2)
lines(p_seq,Vinv(x = S0p_seq, beta = beta_temp, type = "logplusone")/S0p_seq, 
      type = "l", col=four_color_qual[4], lwd=2)
dev.off()
