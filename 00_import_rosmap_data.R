#read in packages
sapply(c("readxl","tidyverse","tidylog"),
       function(x){if(!suppressWarnings(require(package=x,character.only = T)))
       {install.packages(pkgs=x);require(package=x,character.only = T)}})

ROSMAPpath <- "/Users/reederh/Documents/Projects/ROSMAP/"
ROSMAPinput <- paste0(ROSMAPpath,"Input/")
ROSMAPtemp <- paste0(ROSMAPpath,"Temp/")
ROSMAPscripts <- "/Users/reederh/Dropbox/Harrison Reeder/ROSMAP/Scripts/"
# ROSMAPoutput <- paste0(path,"Output/")

##*******************************##
####import baseline ROSMAP data####
##*******************************##

rosmap_person <- read_xlsx(path = paste0(ROSMAPinput,"dataset_1101_basic_11-10-2021.xlsx"),sheet="Sheet0")
rosmap_visit <- read_xlsx(path = paste0(ROSMAPinput,"dataset_1101_long_11-10-2021.xlsx"),sheet="Sheet0")

#INSPECT AVAILABLE COVARIATES
colnames(rosmap_person)
colnames(rosmap_visit)
summary(rosmap_person)
summary(rosmap_visit)

#make lookup of all people who are enrolled with AD
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
  filter(!(studyid %in% ad_baseline_ids)) #remove those 218 subjects with AD at baseline

rosmap_visit_cleaned <- rosmap_visit %>% 
  filter(!is.na(age_at_visit)) %>% #filter the 3 visit obs with missing age, and also like everything is missing
  mutate(studyid = paste0(study,"_",projid),
         diag_long_fct = recode(dcfdx,`1`="NCI",`2`="MCI",`3`="MCI+",`4`="AD",`5`="AD+",`6`="Other")) %>%
  dplyr::select(-study,-projid,-scaled_to) %>%
  filter(!(studyid %in% ad_baseline_ids)) #remove those with AD at baseline, totalling 1000 observations (3%)

#fix discrepancy with baseline age of a single individual (slightly off)
rosmap_person_cleaned$age_bl[rosmap_person_cleaned$studyid=="ROS_20723662"] <-
  rosmap_visit_cleaned$age_at_visit[rosmap_visit_cleaned$studyid=="ROS_20723662" & rosmap_visit_cleaned$fu_year==1]

#fix discrepancy with subject who has diagnosis age after last visit age
rosmap_person_cleaned$age_first_ad_dx[rosmap_person_cleaned$studyid=="ROS_74105189"] <- NA

#fix discrepancy with age of subjects that seem off by 10
rosmap_visit_cleaned$age_at_visit[rosmap_visit_cleaned$studyid=="MAP_06152191" & 
                                    rosmap_visit_cleaned$fu_year==13] <- 
  rosmap_visit_cleaned$age_at_visit[rosmap_visit_cleaned$studyid=="MAP_06152191" & 
                                      rosmap_visit_cleaned$fu_year==13] + 10
rosmap_visit_cleaned$age_at_visit[rosmap_visit_cleaned$studyid=="MAP_94734696" & 
                                    rosmap_visit_cleaned$fu_year==14] <- 
  rosmap_visit_cleaned$age_at_visit[rosmap_visit_cleaned$studyid=="MAP_94734696" & 
                                      rosmap_visit_cleaned$fu_year==14] + 10

##********************************##
####MAKING A LONG FORMAT DATASET####
##********************************##

#first, merge on a flag indicating which visit was the enrollment visit
#based on matching "age at enrollment" in the two files
rosmap_visit_enroll <- left_join(x = rosmap_visit_cleaned, 
                                  y = rosmap_person_cleaned %>% 
                                    dplyr::select(studyid, age_bl) %>% 
                                    mutate(enroll = 1) %>%
                                    rename(age_at_visit=age_bl), 
                                  by = c("studyid","age_at_visit"))

#next, merge on a flag for the visit at which diagnosis occurred
  #there was previously a single case where a person has a diagnosis age with no corresponding visit age, and that row is just added.
  #that doesn't happen anymore, so we're all good!
rosmap_visit_diag <- left_join(x = rosmap_visit_enroll, 
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
         outcome_cr = ifelse(ad_diagnosis,1,ifelse(death,2,0)),
         age_first_ad_dx_mid = (age_prev_first_ad_dx + age_first_ad_dx) / 2,
         #create categorical outcome variable for ad, death, both, neither
         outcome_cat = recode(death + ad_diagnosis, `2`="both", `1`="death_only", `0`="neither"),
         outcome_cat = as.factor(ifelse(outcome_cat=="death_only" & ad_diagnosis==1,
                                        "ad_diagnosis_only", outcome_cat))) %>%
  relocate(studyid, age_bl,age_prev_first_ad_dx,age_first_ad_dx_mid,age_first_ad_dx,age_last,ad_diagnosis,death,outcome_cr)

rosmap_baseline_temp <- bind_cols(rosmap_baseline_temp,
                             model.matrix(~ 0 + marital_bl_fct, rosmap_baseline_temp) %>% as.data.frame(),
                             model.matrix(~ 0 + smoking_fct, rosmap_baseline_temp) %>% as.data.frame(),
                             model.matrix(~ 0 + race_fct, rosmap_baseline_temp) %>% as.data.frame())

#no one is diagnosed and dead at the exact time which is good
# rosmap_baseline_temp %>% filter(age_first_ad_dx >= age_last_visit & death==1) %>% View

# #also, age previous to diagnosis age is NA iff age of diagnosis is NA, which is good!
stopifnot(is.na(rosmap_baseline_temp$age_prev_first_ad_dx) == is.na(rosmap_baseline_temp$age_first_ad_dx))


##*******************************************************##
####MAKE FINAL BASELINE DATASET WITH SELECTION SETTINGS####
##*******************************************************##

exclusion_age_low <- 65 #all subjects must be at least this age 
exclusion_age_high <- 86 #all subjects must be enrolled by this age

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

# #the visit intervals containing AD onset are almost all 1 year,
# #so interval censoring does not seem a major issue here.
# #nevertheless, I do midpoint imputation just to make it all official.
# hist(rosmap_baseline$age_first_ad_dx - rosmap_baseline$age_prev_first_ad_dx,breaks=100)

##************************************************************##
####MAKE SIMPLE LONGITUDINAL DATASET WITH SELECTION SETTINGS####
##************************************************************##

rosmap_longitudinal_temp <- rosmap_visit_cleaned %>% 
  left_join(x = rosmap_baseline_temp, by="studyid") %>%
  mutate(
    ad_at_risk=ifelse(!is.na(age_first_ad_dx),
                                  as.numeric(age_at_visit < age_first_ad_dx), 1),
         death_at_risk=ifelse(!is.na(age_first_ad_dx),
                           as.numeric(age_at_visit <= age_first_ad_dx), 1),
         time_last = age_last - exclusion_age_low,
         time_ad_dx = age_first_ad_dx - exclusion_age_low,
         time_ad_dx = ifelse(is.na(time_ad_dx),time_last,time_ad_dx),
         time_sojourn = ifelse(ad_at_risk==0,age_at_visit-age_first_ad_dx,0),
         time_sojourn_last = pmax(0,time_last - time_ad_dx),
         time_bl = age_bl - exclusion_age_low,
         time_at_visit = age_at_visit - exclusion_age_low) %>% 
  filter(time_bl >= 0) %>% #set `baseline year` as dplyr::selected above
  filter(time_bl < time_last) %>% #filter to only those whose baseline visit is not also their last visit
  dplyr::select(studyid, time_bl, time_at_visit,time_ad_dx,
                time_sojourn,time_sojourn_last,
                time_last,ad_at_risk,
                #death_at_risk,
                ad_diagnosis,death,outcome_cr,outcome_cat,
                race_fctwhite,msex,education_15plus,cogn_global) %>%
  group_by(studyid) %>% fill(cogn_global, .direction="down") %>%
  mutate(na_cogn_global=sum(is.na(cogn_global)),
         na_education_15plus=sum(is.na(education_15plus))) %>%
  filter(na_cogn_global==0,na_education_15plus==0) %>% 
  dplyr::select(-na_cogn_global,-na_education_15plus) %>% ungroup %>% 
  mutate(studyid_num = as.numeric(as.factor(studyid)))
View(rosmap_longitudinal_temp)

rosmap_split_cr <- rosmap_longitudinal_temp %>% 
  filter(time_at_visit < time_ad_dx) %>% 
  group_by(studyid) %>%
  mutate(Tstart = time_at_visit,
         Tstop = lead(time_at_visit),
         Tstop = ifelse(is.na(Tstop),time_ad_dx,Tstop),
         status = ifelse(Tstop==time_ad_dx,outcome_cr,0)) %>%
  filter(Tstart < Tstop) %>% ungroup

rosmap_split_sojourn <- rosmap_longitudinal_temp %>% 
  filter(ad_at_risk==0 & time_sojourn_last>0) %>% 
  group_by(studyid) %>%
  mutate(Tstart = time_sojourn,
         Tstop = lead(time_sojourn),
         Tstop = ifelse(is.na(Tstop),time_sojourn_last,Tstop),
         status = ifelse(Tstop==time_sojourn_last,1,0)) %>%
  filter(Tstart < Tstop) %>% ungroup
