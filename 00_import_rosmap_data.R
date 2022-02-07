#read in packages
sapply(c("readxl","tidyverse","tidylog"),
       function(x){if(!suppressWarnings(require(package=x,character.only = T)))
       {install.packages(pkgs=x);require(package=x,character.only = T)}})

##*******************************##
####import baseline ROSMAP data####
##*******************************##

rosmap_person <- read_xlsx(path = "dataset_1101_basic_11-10-2021.xlsx",sheet="Sheet0")
rosmap_visit <- read_xlsx(path = "dataset_1101_long_11-10-2021.xlsx",sheet="Sheet0")

#INSPECT AVAILABLE COVARIATES
colnames(rosmap_person)
colnames(rosmap_visit)
summary(rosmap_person)
summary(rosmap_visit)
# #look for any weirdo values
# apply(rosmap_person,MARGIN = 2,FUN = unique)
# #look for any weirdo values
# apply(rosmap_visit,MARGIN = 2,FUN = unique)

##VISUALIZE MISSINGNESS
# library(visdat)
# vis_miss(rosmap_person, sort_miss = FALSE,cluster = TRUE)
# vis_miss(rosmap_person, sort_miss = TRUE,cluster = TRUE)
# # there are too many covariates to visualize missingness in visit-level data
# vis_miss(rosmap_visit)

#3 visits are missing an age at visit...
apply(rosmap_visit,MARGIN = 2,FUN = function(x){sum(is.na(x))})

# #VISUALIZING OUTCOME INFORMATION##
# ggplot(data=rosmap_person, mapping = aes(x=age_bl, y=age_death)) +
#   geom_jitter() + theme_classic() +
#   labs(title = "Age at Death by Age at Enrollment")
# ggplot(data=rosmap_person, mapping = aes(x=age_bl, y=age_first_ad_dx)) +
#   geom_jitter() + theme_classic() +
#   labs(title = "Age at AD diagnosis by Age at Enrollment")
# plot(ecdf(rosmap_person$age_bl))

# #Yikes turns out several hundred people are enrolled with AD, and then are missing a time of AD dx
# rosmap_visit %>% filter(fu_year == 0 & dcfdx > 1) %>% View
# #make lookup of all people who are enrolled with AD
ad_baseline_ids <- rosmap_visit %>% filter(fu_year==0 & dcfdx>=4) %>%
  mutate(studyid = paste0(study,"_",projid)) %>% dplyr::select(studyid) %>% as.matrix %>% as.vector

##*************************************************##
####MAKE SIMPLE CLEANED DATASETS FOR FURTHER WORK####
##*************************************************##
#note that I resolve missingness in smoking status and marital status with missingness indicators!!

rosmap_person_cleaned <- rosmap_person %>%
  mutate(studyid = paste0(study,"_",projid),
         #replace missing values for certain baseline factors with an unknown indicator
         marital_now_bl = replace_na(marital_now_bl,6), #46% missing
         smoking = replace_na(smoking,3), #less than 1% missing
         marital_bl_fct= recode(marital_now_bl,`1`="never",`2`="married",`3`="widowed",
                          `4`="divorced",`5`="separated", `6`="missing"),
         evermarried = ifelse(marital_bl_fct %in% c("married","divorced","separated","widowed"),1,0),
         wasmarried = ifelse(marital_bl_fct %in% c("divorced","separated","widowed"),1,0),
         smoking_fct= recode(smoking,`0`="never",`1`="former",`2`="current", `3`="missing"),
         education_15plus = ifelse(educ>=15,1,0),
         eversmoking = ifelse(smoking_fct %in% c("current","former"),1,0),
         diag_death_fct = recode(cogdx,`1`="NCI",`2`="MCI",`3`="MCI+",`4`="AD",`5`="AD+",`6`="Other"),
         race_fct= recode(race7,`1`="white",`2`="black",`3`="american_indian_alaskan_native",
                          `4`="native_hawaiian_pacific_islander",`5`="asian",`6`="other",`7`="unknown"),
         hisp_bin=ifelse(spanish==2,0,1),
         ceradsc_fct = recode(ceradsc,`1`="definiteAD",`2`="probableAD",`3`="possibleAD",`4`="noAD"),
         ceradsc_ad_bin = ifelse(ceradsc<=2,1,0),
         niareagansc_fct = recode(niareagansc,`1`="highADlikelihood",`2`="intADlikelihood",`3`="lowADlikelihood",`4`="noAD"),
         niareagansc_ad_bin = ifelse(niareagansc<=2,1,0)) %>%
  separate(col=apoe_genotype,sep=1,into=c("apoe1","apoe2")) %>%
  mutate(apoe4_any = as.numeric(pmax(apoe1,apoe2)==4), #create meaningful APOE variables
         apoe4_ct = as.numeric(apoe1==4) + as.numeric(apoe2==4)) %>%
  dplyr::select(-study,-projid,-scaled_to,-apoe1,-apoe2) %>%
  # filter(studyid != "ROS_20723662") %>% #remove one subject with no visit-level data collected!
  filter(!(studyid %in% ad_baseline_ids)) #remove those 218 subjects with AD at baseline

rosmap_visit_cleaned <- rosmap_visit %>% 
  filter(!is.na(age_at_visit)) %>% #filter the 3 visit obs with missing age, and also like everything is missing
  mutate(studyid = paste0(study,"_",projid),
         diag_long_fct = recode(dcfdx,`1`="NCI",`2`="MCI",`3`="MCI+",`4`="AD",`5`="AD+",`6`="Other")) %>%
  dplyr::select(-study,-projid,-scaled_to) %>%
  filter(!(studyid %in% ad_baseline_ids)) #remove those with AD at baseline, totalling 1000 observations (3%)

##********************************##
####MAKING A LONG FORMAT DATASET####
##********************************##

#first, merge on a flag indicating which visit was the enrollment visit
#based on matching "age at enrollment" in the two files
#ONE discrepancy who is in person level data but not visit data, so keep track of that
  #ROS_20723662
rosmap_visit_enroll <- left_join(x = rosmap_visit_cleaned, 
                                  y = rosmap_person_cleaned %>% 
                                    dplyr::select(studyid, age_bl) %>% 
                                    mutate(enroll = 1) %>%
                                    rename(age_at_visit=age_bl), 
                                  by = c("studyid","age_at_visit"))

#finding the discrepant participant with person data but no visit data
# rosmap_visit_enroll %>% group_by(studyid) %>% 
#   summarize(n_enroll=sum(enroll,na.rm=TRUE)) %>% filter(n_enroll == 0)
  
  
#next, merge on a flag for the visit at which diagnosis occurred
  #for the single case where a person has a diagnosis age with no corresponding visit age, that row is just added.
  #hence, the "full" join
rosmap_visit_diag <- full_join(x = rosmap_visit_enroll, 
                                y = rosmap_person_cleaned %>% 
                                      dplyr::select(studyid, age_first_ad_dx) %>% 
                                      filter(!is.na(age_first_ad_dx)) %>%
                                      mutate(ad_diagnosis = 1) %>%
                                      rename(age_at_visit=age_first_ad_dx), 
                                by = c("studyid","age_at_visit"))

#finally, merge on a flag for which visit death occurs
#when the death age does not correspond with a visit age, add it in as a new row.
#again, hence the "full" join
rosmap_visit_death <- full_join(x=rosmap_visit_diag,
                             y=rosmap_person_cleaned %>% 
                               dplyr::select(studyid, age_death) %>% 
                               filter(!is.na(age_death)) %>%
                               mutate(death = 1) %>% 
                               rename(age_at_visit=age_death),
                             by =c("studyid","age_at_visit"))

#now, let's also make a more refined person-level dataset with additional baseline
#and diagnosis interval information
rosmap_baseline_temp <- rosmap_person_cleaned %>%
  #for patients with a diagnosis, add on the age at 'previous' visit
  left_join(rosmap_visit_diag %>% 
              mutate(age_prev_first_ad_dx = lag(age_at_visit)) %>%
              filter(ad_diagnosis==1) %>% 
              dplyr::select(studyid,age_prev_first_ad_dx) %>% distinct,
            by = c("studyid")) %>%
  #add on age at final visit/encounter
  left_join(rosmap_visit_death %>% 
              group_by(studyid) %>% arrange(studyid,age_at_visit) %>% 
              mutate(index=row_number()) %>% filter(index==n()) %>% ungroup %>% 
              rename(age_last = age_at_visit) %>% 
              dplyr::select(studyid,age_last) %>% distinct,
            by = c("studyid")) %>%
  #add on visit-level covariates measured at baseline
  left_join(rosmap_visit_enroll %>% filter(enroll==1) %>% 
              rename_with(~paste0(., "_bl"), !studyid), #add _bl suffix to longitudinals measured at entry
            by = c("studyid")) %>% 
  dplyr::select(-enroll_bl,-fu_year_bl,-age_at_visit_bl) %>% #ditch the unhelpful "baseline" variables
  mutate(mci_bl = ifelse(dcfdx_bl>1,1,0), #indicator for MCI at baseline
         death = ifelse(is.na(age_death),0,1),
         ad_diagnosis = ifelse(is.na(age_first_ad_dx),0,1),
         age_first_ad_dx_mid = (age_prev_first_ad_dx + age_first_ad_dx) / 2,
         #create categorical outcome variable for ad, death, both, neither
         outcome_cat = recode(death + ad_diagnosis, `2`="both", `1`="death_only", `0`="neither"),
         outcome_cat = as.factor(ifelse(outcome_cat=="death_only" & ad_diagnosis==1,
                                        "ad_diagnosis_only", outcome_cat))) %>%
  relocate(studyid, age_bl,age_prev_first_ad_dx,age_first_ad_dx_mid,age_first_ad_dx,age_last,ad_diagnosis,death)

#no one is diagnosed and dead at the exact time which is good
# rosmap_baseline_temp %>% filter(age_first_ad_dx >= age_last_visit & death==1) %>% View

# #also, age previous to diagnosis age is NA iff age of diagnosis is NA, which is good!
# stopifnot(is.na(rosmap_baseline_temp$age_prev_first_ad_dx) == is.na(rosmap_baseline_temp$age_first_ad_dx))


##*******************************************************##
####MAKE FINAL BASELINE DATASET WITH dplyr::selectION SETTINGS####
##*******************************************************##

exclusion_age_low <- 65 #all subjects must be at least this age 
exclusion_age_high <- 86 #all subjects must be enrolled by this age
# censoring_age <- 100 #all followup is censored at this age

rosmap_baseline <- rosmap_baseline_temp %>%
  mutate(time_last = age_last - exclusion_age_low,
         time_ad_dx = age_first_ad_dx - exclusion_age_low,
         time_prev_ad_dx = age_prev_first_ad_dx - exclusion_age_low,
         time_ad_dx_mid = age_first_ad_dx_mid - exclusion_age_low,
         time_ad_dx = ifelse(is.na(time_ad_dx),time_last,time_ad_dx),
         time_ad_dx_mid = ifelse(is.na(time_ad_dx_mid),time_last,time_ad_dx_mid),
         time_bl = age_bl - exclusion_age_low,
         time_lt = as.numeric(time_bl>0)) %>% 
  filter(time_bl >= 0) %>% #set `baseline year` as dplyr::selected above
  filter(age_bl < exclusion_age_high) %>% #set maximum enrollment age as dplyr::selected above
  filter(time_bl < time_last) %>% #filter to only those whose baseline visit is not also their last visit
  relocate(studyid, time_bl,time_prev_ad_dx,time_ad_dx_mid,time_ad_dx,time_last,ad_diagnosis,death)

#add dummy variables for categorical outcomes
rosmap_baseline <- bind_cols(rosmap_baseline,
                           model.matrix(~ 0 + marital_bl_fct, rosmap_baseline) %>% as.data.frame(),
                           model.matrix(~ 0 + smoking_fct, rosmap_baseline) %>% as.data.frame(),
                           model.matrix(~ 0 + race_fct, rosmap_baseline) %>% as.data.frame())

# #the visit intervals containing AD onset are almost all 1 year,
# #so interval censoring does not seem a major issue here.
# #nevertheless, I do midpoint imputation just to make it all official.
# hist(rosmap_baseline$age_first_ad_dx - rosmap_baseline$age_prev_first_ad_dx,breaks=100)
