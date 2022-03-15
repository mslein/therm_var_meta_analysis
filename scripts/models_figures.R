#####################################################
## Manuscript coding, models, and figures ###########
#####################################################
##loading libraries 
pacman::p_load(metafor, tidyverse, viridis, visreg, forcats, 
               devtools, patchwork, MuMIn)

#load data
acclim <- read_csv("data/acclim_cleaned.csv")%>%
  mutate(environment=as.factor(environment), size=as.factor(size), 
         exp_age=as.factor(exp_age), tpc_zone = as.factor(tpc_zone))

acute <- read_csv("data/acute_cleaned.csv")%>%
  mutate(metric= as.factor(metric), environment=as.factor(environment), size=as.factor(size), 
         exp_age=as.factor(exp_age), org_level=as.factor(org_level))

######manuscript statistics on papers##########
#determining overlap between the WOS and SCOPUS search
wos <- read_csv("wos176.csv")
scopus <- read_csv("scopus_405.csv")
comp<-inner_join(scopus,wos, by = "title")

#studies that featured data in both analyses
merge(acclim, acute, by = "study_id") %>%
  count(study_id)
count(acute,environment,genus)
count(acclim,environment,genus)
count(acute,org_level,study_id)

##################################
############# MODELS #############
##################################

##ACCLIMATION
#model selection: adding in moderators to see if they improve AICc by more than 2
#best model by AIC 

#################################################################
##ACCLIMATION
#model selection: sequentially adding in moderators to see if they improve AICc by more than 2

acclim_int<-rma.mv(yi, vi, data=acclim, 
                   random = ~1 |  study_id/unique_species/response_id,
                   method="REML")
#temp
acclim_1<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 
#tpc zone + temp
acclim_2<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(tpc_zone), ref = "neutral"), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML")
#age + temp
acclim_3<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(exp_age), ref="1"), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML")
#size + temp
acclim_4<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(size), ref="1"), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML")
#environment + temp
acclim_5<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +environment, 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML")
#size + tpc zone + temp
acclim_6<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(tpc_zone), ref = "neutral") +relevel(as.factor(size), ref="1"), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 
#tpc zone + age + temp
acclim_7<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(tpc_zone), ref = "neutral") +relevel(as.factor(exp_age), ref="1"), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML")
#tpc zone + environment + temp
acclim_8<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(tpc_zone), ref = "neutral") +environment, 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML")
#age + size + temp
acclim_9<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(exp_age), ref="1") +relevel(as.factor(size), ref="1"), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 
#size + environment + temp
acclim_9<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +environment +relevel(as.factor(size), ref="1"), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 
#environment + age + temp
acclim_10<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +environment +relevel(as.factor(exp_age), ref="1"), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML")
#tpc zone + size + age + temp
acclim_11<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(tpc_zone), ref = "neutral") +relevel(as.factor(exp_age), ref="1")
                 +relevel(as.factor(size), ref="1"), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 
#tpc zone + size + environment + temp 
acclim_12<-rma.mv(yi, vi, data=acclim, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                  +relevel(as.factor(tpc_zone), ref = "neutral") +relevel(as.factor(size), ref="1") 
                  +environment, 
                  random = ~1 |  study_id/unique_species/response_id,
                  method="REML") 
#age + size + environment + temp 
acclim_13<-rma.mv(yi, vi, data=acclim, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                  +relevel(as.factor(exp_age), ref="1") +relevel(as.factor(size), ref="1")
                  +environment, 
                  random = ~1 |  study_id/unique_species/response_id,
                  method="REML") 
#tpc zone + age + environment + temp
acclim_14<-rma.mv(yi, vi, data=acclim, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                  +relevel(as.factor(tpc_zone), ref = "neutral") +relevel(as.factor(exp_age), ref="1")
                  +environment, 
                  random = ~1 |  study_id/unique_species/response_id,
                  method="REML") 
#tpc zone + age + size + environment + temp 
acclim_15<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(tpc_zone), ref = "neutral") +relevel(as.factor(exp_age), ref="1")
                 +relevel(as.factor(size), ref="1") +environment, 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 

#generating stats for table s3
#getting aic and df for each model
AIC<-AICc(acclim_int, acclim_1, acclim_2, acclim_3, acclim_4, acclim_5, acclim_6, acclim_7,acclim_8,
     acclim_9, acclim_10, acclim_11, acclim_12, acclim_13, acclim_14, acclim_15) %>%
  mutate(deltaAICc = AICc - min(AICc))
#getting loglikelihood from each model
log_like <- c((logLik(acclim_int)[1]), (logLik(acclim_1)[1]),
              (logLik(acclim_2)[1]),(logLik(acclim_3)[1]),
                      (logLik(acclim_4)[1]),(logLik(acclim_5)[1]), 
                      (logLik(acclim_6)[1]),(logLik(acclim_7)[1]), 
                   (logLik(acclim_8)[1]), (logLik(acclim_9)[1]), 
                   (logLik(acclim_10)[1]),(logLik(acclim_11)[1]),
                   (logLik(acclim_12)[1]),(logLik(acclim_13)[1]), 
                   (logLik(acclim_14)[1]),(logLik(acclim_15)[1])) %>%
  as.data.frame() %>%
  rename(loglike= ".")
#cleaned table with aic, df, and loglik all together
acclim_stats <- cbind(AIC, log_like)%>%
  mutate(k=332, 
         model = c(0:15)) %>%
  arrange(desc(deltaAICc)) %>%
  round(digits=3)
write.csv(acclim_stats, "acclim_stats.csv", row.names = F)
#best model by AIC 
acclim_best<-rma.mv(yi, vi, data=acclim, 
                    mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                    +relevel(as.factor(tpc_zone), ref = "neutral") +relevel(as.factor(exp_age), ref="1")
                    +relevel(as.factor(size), ref="1") +environment, 
                    random = ~1 |  study_id/unique_species/response_id,
                    method="REML") 


##ACUTE 
#model selection: sequentially adding in moderators to see if they improve AICc by more than 2
acute_int <- rma.mv(yi, vi, data=acute,
                    random = ~1 | study_id/species/response_id,
                    method="REML") 
#temp
acute_1 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant)),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
#size +temp
acute_2 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +relevel(as.factor(size), ref="1"),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
#metric +temp
acute_3 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +relevel(as.factor(metric), ref="energetics"),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
#age +temp
acute_4 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +relevel(as.factor(exp_age), ref="1"),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
#org level +temp
acute_5 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +as.factor(org_level),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
#environment +temp
acute_6 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +environment,
                  random = ~1 | study_id/species/response_id,
                  method="REML") 

#age + environment +temp
acute_7 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +environment +relevel(as.factor(exp_age), ref="1"),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
#org + age  +temp
acute_8 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +as.factor(org_level) +relevel(as.factor(exp_age), ref="1"),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
#size + age  +temp
acute_9 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +relevel(as.factor(size), ref="1") +relevel(as.factor(exp_age), ref="1"),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
#metric + age  +temp
acute_10 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +relevel(as.factor(metric), ref="energetics") +relevel(as.factor(exp_age), ref="1"),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
#metric + size  +temp
acute_11 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(metric), ref="energetics") +relevel(as.factor(size), ref="1"),
                   random = ~1 | study_id/species/response_id,
                   method="REML") 
#environment + size  +temp
acute_12 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +environment +relevel(as.factor(size), ref="1"),
                   random = ~1 | study_id/species/response_id,
                   method="REML") 
#org + size  +temp
acute_13 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +as.factor(org_level) +relevel(as.factor(size), ref="1"),
                   random = ~1 | study_id/species/response_id,
                   method="REML") 
#environment + org  +temp
acute_14 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +environment +as.factor(org_level),
                   random = ~1 | study_id/species/response_id,
                   method="REML") 
#metric + org  +temp
acute_15 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(metric), ref="energetics") +as.factor(org_level),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
#metric + environment  +temp
acute_16 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(metric), ref="energetics") +environment,
                   random = ~1 | study_id/species/response_id,
                   method="REML")
#size + age + environment  +temp
#convergence error; modified convergence
acute_17 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(size), ref="1") +relevel(as.factor(exp_age), ref="1") +environment,
                   random = ~1 | study_id/species/response_id,
                  control = list(optimizer="optim", optmethod="Nelder-Mead"))
#size + age + metric  +temp
acute_18 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(size), ref="1") +relevel(as.factor(exp_age), ref="1") +relevel(as.factor(metric), ref="energetics"),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
#size + age + org  +temp
acute_19 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(size), ref="1") +relevel(as.factor(exp_age), ref="1") +as.factor(org_level),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
#size + environment + org +temp
acute_20 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(size), ref="1") +environment +as.factor(org_level),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
#size + environment + metric  +temp
acute_21 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(size), ref="1") +environment +relevel(as.factor(metric), ref="energetics"),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
# age + environment + metric  +temp
acute_22 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(exp_age), ref="1") +environment +relevel(as.factor(metric), ref="energetics"),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
# age + environment + org +temp
acute_23 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(exp_age), ref="1") +as.factor(org_level) +environment,
                   random = ~1 | study_id/species/response_id,
                   method="REML")
# age + org + metric  +temp
acute_24 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(exp_age), ref="1") +as.factor(org_level) +relevel(as.factor(metric), ref="energetics"),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
# size + metric + org +temp
acute_25 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(size), ref="1") +as.factor(org_level) +relevel(as.factor(metric), ref="energetics"),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
# age + size + environment + metric  +temp
acute_26 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(size), ref="1") +relevel(as.factor(exp_age), ref="1") +environment
                   +relevel(as.factor(metric), ref="energetics"),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
# age + size + environment + org  +temp
acute_27 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(size), ref="1") +relevel(as.factor(exp_age), ref="1") +environment
                   +as.factor(org_level),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
# age + environment + metric + org  +temp
acute_28 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +environment +relevel(as.factor(exp_age), ref="1") +relevel(as.factor(metric), ref="energetics")
                   +as.factor(org_level),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
# age + size + metric + org  +temp
acute_29 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(size), ref="1") +relevel(as.factor(exp_age), ref="1") 
                   +relevel(as.factor(metric), ref="energetics")
                   +as.factor(org_level),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
# size + environment + metric + org  +temp
acute_30 <- rma.mv(yi, vi, data=acute, 
                   mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +relevel(as.factor(size), ref="1") +environment
                   +relevel(as.factor(metric), ref="energetics")
                   +as.factor(org_level),
                   random = ~1 | study_id/species/response_id,
                   method="REML")
# size + environment + metric + org  +temp + org_level
acute_31 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +environment +relevel(as.factor(size), ref="1") +
                    relevel(as.factor(metric), ref="energetics") +relevel(as.factor(exp_age), ref="1") + as.factor(org_level),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
#generating stats for table s4
#getting aic and df for each model
AIC_acute<-AICc(acute_int, acute_1, acute_2, acute_3, acute_4, acute_5, acute_6, acute_7, acute_8, 
                acute_9, acute_10, acute_11, acute_12, acute_13, acute_14, acute_15, acute_16, acute_17,
                acute_18, acute_19, acute_20, acute_21, acute_22, acute_23, acute_24, acute_25, acute_26, 
                acute_27, acute_28, acute_29, acute_30, acute_31) %>%
  mutate(deltaAICc = AICc - min(AICc))
#getting loglikelihood from each model
log_like_acute <- c((logLik(acute_int)[1]), (logLik(acute_1)[1]),
              (logLik(acute_2)[1]),(logLik(acute_3)[1]),
              (logLik(acute_4)[1]),(logLik(acute_5)[1]), 
              (logLik(acute_6)[1]),(logLik(acute_7)[1]), 
              (logLik(acute_8)[1]), (logLik(acute_9)[1]), 
              (logLik(acute_10)[1]),(logLik(acute_11)[1]),
              (logLik(acute_12)[1]),(logLik(acute_13)[1]), 
              (logLik(acute_14)[1]),(logLik(acclim_15)[1]), 
              (logLik(acute_16)[1]), (logLik(acute_17)[1]), 
              (logLik(acute_18)[1]), (logLik(acute_19)[1]), 
              (logLik(acute_20)[1]), (logLik(acute_21)[1]), 
              (logLik(acute_22)[1]), (logLik(acute_23)[1]), 
              (logLik(acute_24)[1]), (logLik(acute_25)[1]), 
              (logLik(acute_26)[1]), (logLik(acute_27)[1]), 
              (logLik(acute_28)[1]), (logLik(acute_29)[1]), 
              (logLik(acute_30)[1]), (logLik(acute_31)[1])) %>%
  as.data.frame() %>%
  rename(loglike= ".") 
#cleaned table with aic, df, and loglik all together
acute_stats <- cbind(AIC_acute, log_like_acute)%>%
  mutate(k=366, 
         model = c(0:31)) %>%
  arrange(desc(deltaAICc)) %>%
  round(digits=3)
write.csv(acute_stats, "acute_stats.csv",row.names = F)

#best acute model by AIC
acute_best <- rma.mv(yi, vi, data=acute, 
                     mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                     +environment +relevel(as.factor(size), ref="1") +
                       relevel(as.factor(metric), ref="energetics") +relevel(as.factor(exp_age), ref="1") + as.factor(org_level),
                     random = ~1 | study_id/species/response_id,
                     method="REML") 
#################################################################



#I2 calculations from best models 
#study level heterogeneity
acclim_best$sigma2[1]/sum(acclim_best$sigma2[1], acclim_best$sigma2[2], acclim_best$sigma2[3]) #0.1330037
#species level heterogeneity
acclim_best$sigma2[2]/sum(acclim_best$sigma2[1], acclim_best$sigma2[2], acclim_best$sigma2[3]) #0.000
#response heterogeneity 
acclim_best$sigma2[3]/sum(acclim_best$sigma2[1], acclim_best$sigma2[2], acclim_best$sigma2[3]) #0.8669962

#study level heterogeneity
acute_best$sigma2[1]/sum(acute_best$sigma2[1], acute_best$sigma2[2], acute_best$sigma2[3]) #0.1659221
#species level heterogeneity
acute_best$sigma2[2]/sum(acute_best$sigma2[1], acute_best$sigma2[2], acute_best$sigma2[3]) #0.01421167
#response heterogeneity 
acute_best$sigma2[3]/sum(acute_best$sigma2[1], acute_best$sigma2[2], acute_best$sigma2[3]) #0.8198662




##################################
###### MAIN TEXT FIGURES #########
##################################
#####coefficient figures 
##acclimation coefficients
#Figure 4
acclim_mod_estimate<- as.data.frame(coef(summary(acclim_best))) %>%
  tibble::rownames_to_column("covariates") %>%
  rename(SMD = "estimate") %>%
  mutate(type = case_when(covariates %in% c('relevel(as.factor(exp_age), ref = "1")0', 'relevel(as.factor(exp_age), ref = "1")2') ~ "Age", 
                          covariates %in% c('relevel(as.factor(size), ref = "1")2') ~ "Body size",
                          covariates %in% c('relevel(as.factor(tpc_zone), ref = "neutral")low', 
                                            'relevel(as.factor(tpc_zone), ref = "neutral")high') ~ "TPC Zone", 
                          covariates %in% c("I(flux_range - mean(flux_range))", "I(mean_temp_reared - mean(mean_temp_reared))",
                                            "I(flux_range - mean(flux_range)):I(mean_temp_reared - mean(mean_temp_reared))") ~ "Temperature", 
                          covariates %in% c("environmentterrestrial") ~ "Environment",
                          TRUE ~ "Intercept")) 

acclim_df <- data.frame(covariates=c("intrcpt", 
                                     "I(flux_range - mean(flux_range))", 
                                     "I(flux_range - mean(flux_range)):I(mean_temp_reared - mean(mean_temp_reared))", 
                                     "relevel(as.factor(tpc_zone), ref = \"neutral\")high", 
                                     "relevel(as.factor(tpc_zone), ref = \"neutral\")low"),
                        SMD=c(0.75, 0.75, 0.75, 3,3))

acclim_coef_fig <-acclim_mod_estimate %>%
  ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_errorbar(aes(x = factor(covariates, level = c("intrcpt", "I(flux_range - mean(flux_range))", 
                                                     "I(mean_temp_reared - mean(mean_temp_reared))", 
                                                     "I(flux_range - mean(flux_range)):I(mean_temp_reared - mean(mean_temp_reared))", 
                                                     "relevel(as.factor(tpc_zone), ref = \"neutral\")high", 
                                                     "relevel(as.factor(tpc_zone), ref = \"neutral\")low",
                                                     'relevel(as.factor(exp_age), ref = "1")0', 
                                                     'relevel(as.factor(exp_age), ref = "1")2', 
                                                     'relevel(as.factor(size), ref = "1")2', "environmentterrestrial")), 
                    y=SMD, ymin=ci.lb, ymax=ci.ub), width=.2)+
  geom_point(aes(x=covariates, y=SMD, fill = SMD>0), size = 6, show.legend = FALSE, pch=21)+
  scale_fill_manual(values= c("turquoise4", "paleturquoise"))+
  scale_x_discrete(breaks=c("intrcpt", "I(flux_range - mean(flux_range))", 
                            "I(mean_temp_reared - mean(mean_temp_reared))", 
                            "I(flux_range - mean(flux_range)):I(mean_temp_reared - mean(mean_temp_reared))", 
                            "relevel(as.factor(tpc_zone), ref = \"neutral\")high", 
                            "relevel(as.factor(tpc_zone), ref = \"neutral\")low",
                            'relevel(as.factor(exp_age), ref = "1")0', 
                            'relevel(as.factor(exp_age), ref = "1")2', 
                            'relevel(as.factor(size), ref = "1")2', "environmentterrestrial"),
                   labels=c("Intercept", "Fluctuation range", 
                            "Mean temperature", "Mean temperature:Fluctuation range", 
                            "TPC Zone High", "TPC Zone Low", 
                            "Experimental age Class 0", "Experimental age Class 2", 
                            "Size Class 2", "Environment Terrestrial"))+
  ylab("Coefficient estimate")+
  xlab("Covariates")+
  ggtitle("ACCLIMATION")+
  #ylim(-15,8)+
  coord_flip()+
  geom_text(data=acclim_df, aes(x=covariates, y=SMD), label = "*", color = "black", size = 6)+
  theme_bw()  

ggsave("figure4.png", acclim_coef_fig, dpi=700, width=8, height=6)#saving to high res png

##acute coefficients
#Figure 5
acute_mod_estimate<- as.data.frame(coef(summary(acute_best))) %>%
  tibble::rownames_to_column("covariates") %>%
  rename(SMD = "estimate") %>%
  mutate(type = case_when(covariates %in% c('relevel(as.factor(exp_age), ref = "1")0', 'relevel(as.factor(exp_age), ref = "1")2') ~ "Age", 
                          covariates %in% c('relevel(as.factor(size), ref = "1")0', '"relevel(as.factor(size), ref = "1")2', 
                                            'relevel(as.factor(size), ref = "1")3 ') ~ "Body size",
                          covariates %in% c('relevel(as.factor(metric), ref = "energetics")biochemistry ', 
                                            'relevel(as.factor(metric), ref = "energetics")development', 
                                            'relevel(as.factor(metric), ref = "energetics")fecundity ', 
                                            'relevel(as.factor(metric), ref = "energetics")size', 
                                            'relevel(as.factor(metric), ref = "energetics")survivorship') ~ "Metric", 
                          covariates %in% c("I(flux_range - mean(flux_range))", 
                                            "I(mean_temp_constant - mean(mean_temp_constant))",
                                            "I(flux_range - mean(flux_range)):I(mean_temp_constant - mean(mean_temp_constant))") ~ "Temperature", 
                          covariates %in% c("environmentterrestrial") ~ "Environment", 
                          covariates %in% c("as.factor(org_level)1") ~ "Organization level",
                          TRUE ~ "Intercept"))

acute_df <- data.frame(covariates=c("I(flux_range - mean(flux_range)):I(mean_temp_constant - mean(mean_temp_constant))", 
                                    'relevel(as.factor(metric), ref = "energetics")biochemistry'),
                       SMD=c(0.25, 3))


acute_coef_fig <- acute_mod_estimate %>%
  ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_errorbar(aes(x = factor(covariates, level = c("intrcpt", "I(flux_range - mean(flux_range))", 
                                                     "I(mean_temp_constant - mean(mean_temp_constant))", 
                                                     "I(flux_range - mean(flux_range)):I(mean_temp_constant - mean(mean_temp_constant))",
                                                     'relevel(as.factor(metric), ref = "energetics")biochemistry',
                                                     'relevel(as.factor(metric), ref = "energetics")development', 
                                                     'relevel(as.factor(metric), ref = "energetics")fecundity', 
                                                     'relevel(as.factor(metric), ref = "energetics")size', 
                                                     'relevel(as.factor(metric), ref = "energetics")survivorship',
                                                     'relevel(as.factor(exp_age), ref = "1")0', 
                                                     'relevel(as.factor(exp_age), ref = "1")2', 
                                                     'relevel(as.factor(size), ref = "1")0', 
                                                     'relevel(as.factor(size), ref = "1")2', 
                                                     'relevel(as.factor(size), ref = "1")3', "environmentterrestrial", 
                                                     "as.factor(org_level)1")), 
                    y=SMD, ymin=ci.lb, ymax=ci.ub), width=.2)+
  geom_point(aes(x=covariates, y=SMD, fill = SMD>0), size = 6, show.legend = FALSE, pch=21)+
  scale_fill_manual(values= c("violetred4", "rosybrown1"))+
  scale_x_discrete(breaks=c("intrcpt", "I(flux_range - mean(flux_range))", 
                            "I(mean_temp_constant - mean(mean_temp_constant))", 
                            "I(flux_range - mean(flux_range)):I(mean_temp_constant - mean(mean_temp_constant))",
                            'relevel(as.factor(metric), ref = "energetics")biochemistry',
                            'relevel(as.factor(metric), ref = "energetics")development', 
                            'relevel(as.factor(metric), ref = "energetics")fecundity', 
                            'relevel(as.factor(metric), ref = "energetics")size', 
                            'relevel(as.factor(metric), ref = "energetics")survivorship',
                            'relevel(as.factor(exp_age), ref = "1")0', 
                            'relevel(as.factor(exp_age), ref = "1")2', 
                            'relevel(as.factor(size), ref = "1")0', 
                            'relevel(as.factor(size), ref = "1")2', 
                            'relevel(as.factor(size), ref = "1")3', "environmentterrestrial", 
                            'as.factor(org_level)1'),
                   labels=c("Intercept", "Fluctuation range", 
                            "Mean temperature", "Mean temperature:Fluctuation range", "Metric Biochemistry",
                            "Metric Development", "Metric Fecundity", 
                            "Metric Size", "Metric Survivorship",
                            "Experimental age Class 0", "Experimental age Class 2", 
                            "Size Class 0","Size Class 2", "Size Class 3", "Environment Terrestrial", 
                            "Organization Level 1"))+
  xlab("Covariates")+
  ylab("Coefficient estimate")+
  ggtitle("ACUTE")+
  coord_flip()+
  geom_text(data=acute_df, aes(x=covariates, y=SMD), label = "*", color = "black", size = 6)+
  theme_bw()
ggsave("figure5.png", acute_coef_fig, dpi=700, width=8, height=6)
#Figure 6
normalized_acclim <- acclim %>% filter(between(yi, -5, 5))
normalized_acute <- acute  %>% filter(between(yi, -5, 5))
#calculating mean effect sizes by mean temperature and flux range for ggplots
mean_acclim <- normalized_acclim %>%
  group_by(mean_temp_reared, flux_range) %>%
  summarise(yi=mean(yi), 
            vi=mean(vi))
mean_acute<-normalized_acute %>%
  group_by(mean_temp_constant, flux_range) %>%
  summarise(yi=mean(yi), 
            vi=mean(vi))
#interaction of mean temp and flux range for acclimation 
interaction_acclim<-ggplot(normalized_acclim, aes(y=yi, x=mean_temp_reared, fill=flux_range))+
  geom_line(y=0, linetype = "dashed")+
  geom_jitter(alpha = 0.5, pch=21,size=1, width = 0.7, height = 0.7)+
  geom_errorbar(data=mean_acclim, aes(ymin=yi-vi, ymax=yi+vi), size=.7, width=.4, color = "black")+
  geom_point(data = mean_acclim, size=4, pch=21, colour="black")+
  theme_bw()+
  scale_fill_gradient(low = "turquoise4", high = "khaki1")+
  ggtitle("ACCLIMATION")+
  labs(x = "Mean temperature (°C)",
       y ="SMD",
       fill = "Temperature fluctuation range (°C)")+
  theme(legend.position="bottom")
#interaction of mean temp and flux range for acute
interaction_acute<-ggplot(normalized_acute, aes(y=yi, x=mean_temp_constant, fill=flux_range))+
  geom_jitter(alpha = 0.5, pch=21,size=1, width = 0.3, height = 0.7)+
  geom_line(y=0, linetype = "dashed")+
  geom_errorbar(data=mean_acute, aes(ymin=yi-vi, ymax=yi+vi), size=.7, width=.4, color = "black")+
  geom_jitter(data = mean_acute, size=4, pch=21, colour="black")+
  theme_bw()+
  scale_fill_gradient(low = "violetred4", high = "khaki1")+
  ggtitle("ACUTE")+
  labs(x = "Mean temperature (°C)",
       y ="SMD",
       fill = "Temperature fluctuation range (°C)")+
  theme(legend.position="bottom")
#subset of mean temp and flux range for acute
subset_acute<-normalized_acute %>%
  filter(mean_temp_constant %in% c(20,24,30)) %>%
  ggplot(aes(y=yi, fill=mean_temp_constant, x=flux_range))+
  geom_jitter(alpha = 0.5, pch=21,size=3, width = 0.3, height = 0.7)+
  geom_line(y=0, linetype = "dashed")+
  theme_bw()+
  scale_fill_gradient(low = "violetred4", high = "khaki1")+
  ggtitle("ACUTE")+
  labs(x = "Temperature fluctuation range (°C)",
       y ="SMD",
       fill = "Mean Temperature (°C)")+
  theme(legend.position="bottom")+
  geom_smooth(method=lm)
  facet_wrap(~mean_temp_constant)
#subset of mean temp and flux range for acclimation
subset_acclim <-normalized_acclim %>%
  filter(mean_temp_reared %in% c(15, 20, 24)) %>%
  ggplot(aes(y=yi, x=flux_range, fill=mean_temp_reared))+
  geom_jitter(alpha = 0.5, pch=21,size=3, width = 0.3, height = 0.7)+
  geom_line(y=0, linetype = "dashed")+
  theme_bw()+
  scale_fill_gradient(low = "turquoise4", high = "khaki1")+
  ggtitle("ACCLIMATION")+
  labs(x = "Temperature fluctuation range (°C)",
       y ="SMD",
       fill = "Mean Temperature (°C)")+
  theme(legend.position="bottom")+
  geom_smooth(method=lm)+
  facet_wrap(~mean_temp_reared)
#joining the acclimation plots
p_big <- subset_acclim + interaction_acclim 
#joining the acute plots 
p_big2 <- subset_acute + interaction_acute
#making the main text figure 
interaction_fig <- p_big/p_big2
ggsave("figure6.png", interaction_fig, dpi=700, width=10, height=9)

######SUPPLEMENTARY FIGURES#####
#funnel plots 
### draw funnel plots
#figure s1
dummy_acute <- rma.mv(yi, vi, data=normalized_acute, 
                      random = ~1 |  study_id/ response_id,
                      method="REML")
funnel(dummy_acute, main="Acute", col = "violetred4", alpha=0.5)
#figure 2
dummy_acclim <- rma.mv(yi, vi, data=normalized_acclim, 
                       random = ~1 |  study_id/ response_id,
                       method="REML")
funnel(dummy_acclim, main="Acclimation", col = "turquoise4")

#histogram of effect sizes to demonstrate normal distribution
#figure s3
ggplot()+
  geom_histogram(data= normalized_var, aes(x=yi), fill = "violetred4",col = "white")+
  geom_histogram(data= normalized_acclim, aes(x=yi), fill= "turquoise3",col = "white",alpha=0.5)+
  theme_bw()+
  labs(x="Effect size (SMD)", 
       y= "Count")

#full facet of fluctuation range and mean temp for acute
#figure s4
ggplot(normalized_acclim, aes(y=yi, x=flux_range, fill=mean_temp_reared))+
  geom_jitter(alpha = 0.5, pch=21,size=3, width = 0.3, height = 0.7)+
  geom_line(y=0, linetype = "dashed")+
  #geom_errorbar(data=mean_var, aes(ymin=yi-vi, ymax=yi+vi), size=.7, width=.4, color = "black")+
  #geom_jitter(data = mean_var, size=4, pch=21, colour="black")+
  theme_bw()+
  scale_fill_gradient(low = "turquoise4", high = "khaki1")+
  ggtitle("ACCLIMATION")+
  labs(x = "Temperature fluctuation range (°C)",
       y ="SMD",
       fill = "Mean Temperature (°C)")+
  theme(legend.position="bottom")+
  geom_smooth(method=lm)+
  facet_wrap(~mean_temp_reared)

#full facet of fluctuation range and mean temp for acute
#figure s5
ggplot(normalized_acute, aes(y=yi, fill=mean_temp_constant, x=flux_range))+
  geom_jitter(alpha = 0.5, pch=21,size=3, width = 0.3, height = 0.7)+
  geom_line(y=0, linetype = "dashed")+
  #geom_errorbar(data=mean_var, aes(ymin=yi-vi, ymax=yi+vi), size=.7, width=.4, color = "black")+
  #geom_jitter(data = mean_var, size=4, pch=21, colour="black")+
  theme_bw()+
  scale_fill_gradient(low = "violetred4", high = "khaki1")+
  ggtitle("ACUTE")+
  labs(x = "Temperature fluctuation range (°C)",
       y ="SMD",
       fill = "Mean Temperature (°C)")+
  theme(legend.position="bottom")+
  geom_smooth(method=lm)+
  facet_wrap(~mean_temp_constant)


