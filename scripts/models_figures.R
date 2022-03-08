#####################################################
## Manuscript coding, models, and figures ###########
#####################################################
##loading libraries 
pacman::p_load(metafor, tidyverse, viridis, visreg, forcats, 
               devtools, patchwork)

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
merge(dat_acclim_ES, dat_full_var_ES, by = "study_id") %>%
  count(study_id)
count(dat_full_var_ES, environment, genus)
count(dat_acclim_ES, environment, genus)
count(dat_full_var_ES,  org_level,study_id)

##################################
############# MODELS #############
##################################

##ACCLIMATION
#model selection: sequentially adding in moderators to see if they improve AICc by more than 2

acclim_int<-rma.mv(yi, vi, data=acclim, 
                   random = ~1 |  study_id/unique_species/response_id,
                   method="REML") 
acclim_1<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) , 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 
acclim_2<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) + I(mean_temp_reared-mean(mean_temp_reared)), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 
acclim_3<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 
acclim_4<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(tpc_zone), ref = "neutral"), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 
acclim_5<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(tpc_zone), ref = "neutral") +relevel(as.factor(exp_age), ref="1"), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 
acclim_6<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(tpc_zone), ref = "neutral") +relevel(as.factor(exp_age), ref="1")
                 +relevel(as.factor(size), ref="1"), 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 
acclim_7<-rma.mv(yi, vi, data=acclim, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_reared-mean(mean_temp_reared)) 
                 +relevel(as.factor(tpc_zone), ref = "neutral") +relevel(as.factor(exp_age), ref="1")
                 +relevel(as.factor(size), ref="1") +environment, 
                 random = ~1 |  study_id/unique_species/response_id,
                 method="REML") 

AICc(acclim_int, acclim_1, acclim_2, acclim_3, acclim_4, acclim_5, acclim_6, acclim_7)

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
acute_1 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
acute_2 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) + I(mean_temp_constant-mean(mean_temp_constant)),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
acute_3 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant)),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
acute_4 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +environment,
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
acute_5 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +environment +relevel(as.factor(size), ref="1"),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 
acute_6<- rma.mv(yi, vi, data=acute, 
                 mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                 +environment +relevel(as.factor(size), ref="1") +
                   relevel(as.factor(metric), ref="energetics"),
                 random = ~1 | study_id/species/response_id,
                 method="REML") 
acute_7 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +environment +relevel(as.factor(size), ref="1") +
                    relevel(as.factor(metric), ref="energetics") +relevel(as.factor(exp_age), ref="1"),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 

acute_8 <- rma.mv(yi, vi, data=acute, 
                  mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +environment +relevel(as.factor(size), ref="1") +
                    relevel(as.factor(metric), ref="energetics") +relevel(as.factor(exp_age), ref="1") + as.factor(org_level),
                  random = ~1 | study_id/species/response_id,
                  method="REML") 

AICc(acute_int, acute_1, acute_2, acute_3, acute_4, acute_5, acute_6, acute_7, acute_8)

#best model -- organization level doesn't improve AIC by more than 2
acute_best <- rma.mv(yi, vi, data=acute, 
                      mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                      +environment +relevel(as.factor(size), ref="1") +
                       relevel(as.factor(metric), ref="energetics") +relevel(as.factor(exp_age), ref="1"),
                      random = ~1 | study_id/species/response_id,
                      method="REML") 


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

acclim_mod_estimate %>%
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
                          TRUE ~ "Intercept"))

acute_df <- data.frame(covariates=c("I(flux_range - mean(flux_range)):I(mean_temp_constant - mean(mean_temp_constant))", 
                                    'relevel(as.factor(metric), ref = "energetics")biochemistry'),
                       SMD=c(0.25, 3))


acute_mod_estimate %>%
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
                                                     'relevel(as.factor(size), ref = "1")3', "environmentterrestrial")), 
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
                            'relevel(as.factor(size), ref = "1")3', "environmentterrestrial"),
                   labels=c("Intercept", "Fluctuation range", 
                            "Mean temperature", "Mean temperature:Fluctuation range", "Metric Biochemistry",
                            "Metric Development", "Metric Fecundity", 
                            "Metric Size", "Metric Survivorship",
                            "Experimental age Class 0", "Experimental age Class 2", 
                            "Size Class 0","Size Class 2", "Size Class 3", "Environment Terrestrial"))+
  xlab("Covariates")+
  ylab("Coefficient estimate")+
  ggtitle("ACUTE")+
  coord_flip()+
  geom_text(data=acute_df, aes(x=covariates, y=SMD), label = "*", color = "black", size = 6)+
  theme_bw()

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
p_big/p_big2



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


