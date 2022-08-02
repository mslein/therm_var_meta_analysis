#load data
pacman::p_load(tidyverse, metafor, rotl, taxadb, taxize, Hmisc, dmetar, viridis, phytools, patchwork, orchaRd)
acclim <- read_csv("acclimation_cleaned_24may.csv")%>%
  mutate(environment=as.factor(environment), size=as.factor(size), 
         exp_age=as.factor(exp_age), tpc_zone = as.factor(tpc_zone), 
         experiment_type = "acclimation", 
         mean_temp_constant = mean_temp_reared, 
         metric = tpc_zone) %>%
  select(response_id, species, study_id, org_level, 
         mean_temp_constant, exp_age, size, yi, vi,metric, experiment_type, 
         flux_range, environment, genus, resp_def, secondary_temp, duration, duration_units)

acute <- read_csv("acute_cleaned_24may.csv")%>%
  mutate(metric= as.factor(metric), environment=as.factor(environment), size=as.factor(size), 
         exp_age=as.factor(exp_age), org_level=as.factor(org_level), 
         experiment_type = "acute") %>%
  select(response_id, species, study_id, org_level, 
         mean_temp_constant, exp_age, size, yi, vi,metric, experiment_type, 
         flux_range, environment, genus, resp_def, secondary_temp, duration, duration_units)

#excluding the species that aren't in the open tree of life
#there are 668 rows now, there were 698 rows in the original df
# ~4% data loss by removing these species 
match_names <- read_csv("rotl_cor.csv") %>%
  mutate(original_name = tolower(original_name)) 
new_metrics <- read_csv("response_metric_sam_and_maggie.csv") %>%
  select(resp_def, metric_1, metric_2)
metrics_24may <- read_csv("metrics_24may.csv")
full <- rbind(acute, acclim) %>%
  unite(original_name, c(genus, species), sep = " ", remove = FALSE) %>%
  mutate(original_name = tolower(original_name)) %>%
  filter(!original_name %in% c("escherichia coli", 
                               "isotopenola loftyensis", 
                               "phanaeus vindex")) %>%
  left_join(match_names, by="original_name") %>%
  left_join(new_metrics, by="resp_def") %>%
  filter(yi <80) %>%
  mutate(duration_std_days = case_when(duration_units ==  c("weeks") ~ duration*7, 
                                       duration_units ==  c("hours") ~ duration/24, 
                                       TRUE ~ duration), 
         duration_std_days_neg = case_when(experiment_type == c("acclimation") ~ duration_std_days*(-1), 
         experiment_type == c("acute") ~ duration_std_days)) %>%
  filter(!size %in% c(0,3), 
         !org_level %in% c(1)) %>% #39 species instead of 42 
  left_join(metrics_24may, by="resp_def")


acclim_new <- full %>%
  filter(experiment_type == "acclimation")



acute_new <- full %>%
  filter(experiment_type == "acute")

#looks like there is quite a big phylogenetic effect? but this doesn't have the nested effects-- it's just a list 
full_phylo <- rma.mv(yi, vi, 
                     mods = ~ +I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                     +experiment_type +exp_age 
                     +size +environment +org_level
                     +relevel(as.factor(metric_2), ref="physiology"),
                     random = list(~1 | study_id, ~1 | rotl_name, ~1 | response_id),
                     R = list(rotl_name=cor),
                     method="REML", data=full)



#with interactions 
rma.mv(yi, vi, 
       mods = ~ +I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
       +I(secondary_temp - mean(secondary_temp))
       +experiment_type* stressful +experiment_type*duration_std_days_neg +experiment_type*exp_age 
       +experiment_type*size +experiment_type*environment,
       random = list(~1 | study_id, ~1 | rotl_name, ~1 | response_id),
       R = list(rotl_name=cor),
       method="REML", data=full)

#without interactions 
rma.mv(yi, vi, 
       mods = ~ +I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
       +I(secondary_temp - mean(secondary_temp))
       +experiment_type +stressful +duration_std_days_neg +exp_age 
       +size +environment,
       random = list(~1 | study_id, ~1 | rotl_name, ~1 | response_id),
       R = list(rotl_name=cor),
       method="REML", data=full)


#no overall effect of variability 
overall <- rma.mv(yi, vi,,
                  random = list(~1 | study_id, ~1 | rotl_name, ~1 | response_id),
                  R = list(rotl_name=cor),
                  method="REML", data=full)
#significant negative effect of variability 
m1 <- rma.mv(yi, vi, mods= ~experiment_type,
                  random = list(~1 | study_id, ~1 | rotl_name, ~1 | response_id),
                  R = list(rotl_name=cor),
                  method="REML", data=full)
summary(m1)

### correlation between different moderators
cor(as.numeric(full$size), as.numeric(full$exp_age)) #-0.04014608
cor(as.numeric(full$size), as.numeric(full$org_level)) #-0.1803301

full_phylo_df <- as.data.frame(coef(summary(full_phylo))) %>%
  tibble::rownames_to_column("covariates") %>%
  rename(SMD = "estimate") %>%
  mutate(type = case_when(covariates %in% c("exp_age1", "exp_age2") ~ "age", 
                          covariates %in% c("size1", "size2", "size3") ~ "size",
                          covariates %in% c("metricdevelopment", "metricenergetics", "metricfecundity", 
                                            "metricsize", "metricsurvivorship", "metrichigh", "metriclow") ~ "metric", 
                          covariates %in% c("I(flux_range - mean(flux_range))", "I(mean_temp_constant - mean(mean_temp_constant))",
                                            "I(flux_range - mean(flux_range)):I(mean_temp_constant - mean(mean_temp_constant))") ~ "Temperature", 
                          covariates %in% c("environmentterrestrial", "experiment_typeacute", "org_level1") ~ "experimental design",
                          covariates %in% c("vi") ~ "publication bias",
                          TRUE ~ "Intercept")) 

full_phylo_df %>%
  ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_errorbar(aes(x = covariates, y=SMD, ymin=ci.lb, ymax=ci.ub), width=.2)+
  geom_point(aes(x=covariates, y=SMD, fill= SMD > 0), size = 6, show.legend = FALSE, pch=21)+
  scale_fill_viridis_d()+
  ylab("SMD")+
  xlab("Covariates")+
  #ggtitle("ACCLIMATION")+
  #ylim(-15,8)+
  #facet_grid(~type)+
  coord_flip()+
  #geom_text(data=acclim_df, aes(x=covariates, y=SMD), label = "*", color = "black", size = 6)+
  theme_bw()

trait<-full[,7]
names(trait)<-rownames(full)
trait <- as.data.frame(as.numeric(full$exp_age), full$rotl_name)

x<-fastBM(tree)
phylosig(tree, as.matrix(mydat2)[,c("exp_age")])

full %>%
  ggplot(aes(x=yi))+
  geom_histogram()

full %>%
  ggplot(aes(x=yi))+
  geom_histogram()+
  facet_wrap(~metric_2)

#age
a<- full %>%
  filter(yi <80) %>%
  group_by(experiment_type, exp_age) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=exp_age, y=mean))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_linerange(aes(x = exp_age, y=mean, ymin=mean-se, ymax=mean+se))+
  geom_point(aes(x = exp_age, y=mean, fill= mean > 0), size = 6, show.legend = FALSE, pch=21)+
  scale_fill_viridis_d()+
  coord_flip()+
  facet_wrap(~experiment_type)
#size
b <- full %>%
  filter(yi <80) %>%
  group_by(experiment_type, size)%>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=size, y=mean))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_linerange(aes(x = size, y=mean, ymin=mean-se, ymax=mean+se))+
  geom_point(aes(x = size, y=mean, fill= mean > 0), size = 6, show.legend = FALSE, pch=21)+
  scale_fill_viridis_d()+
  coord_flip()+
  facet_wrap(~experiment_type)
#org level
c <- full %>%
  filter(yi <80) %>%
  group_by(org_level) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=org_level, y=mean))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_linerange(aes(x = org_level, y=mean, ymin=mean-se, ymax=mean+se))+
  geom_point(aes(x = org_level, y=mean, fill= mean > 0), size = 6, show.legend = FALSE, pch=21)+
  scale_fill_viridis_d()+
  coord_flip()
#metric
d <- full %>%
  filter(yi <80) %>%
  group_by(metric_2) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=metric_2, y=mean))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_linerange(aes(x = metric_2, y=mean, ymin=mean-se, ymax=mean+se))+
  geom_point(aes(x = metric_2, y=mean, fill= mean > 0), size = 6, show.legend = FALSE, pch=21)+
  scale_fill_viridis_d()+
  coord_flip()
#experiment type
e <- full %>%
  filter(yi <80) %>%
  group_by(experiment_type) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=experiment_type, y=mean))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_linerange(aes(x = experiment_type, y=mean, ymin=mean-se, ymax=mean+se))+
  geom_point(aes(x = experiment_type, y=mean, fill= mean > 0), size = 6, show.legend = FALSE, pch=21)+
  scale_fill_viridis_d()+
  coord_flip()+
  facet_wrap(~experiment_type)
#environment
f <-full %>%
  filter(yi <80) %>%
  group_by(experiment_type, environment) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=environment, y=mean))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_linerange(aes(x = environment, y=mean, ymin=mean-se, ymax=mean+se))+
  geom_point(aes(x = environment, y=mean, fill= mean > 0), size = 6, show.legend = FALSE, pch=21)+
  scale_fill_viridis_d()+
  coord_flip()+
  facet_wrap(~experiment_type)

g <- full %>%
  filter(yi <80) %>%
  group_by(experiment_type, stressful) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=stressful, y=mean))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_linerange(aes(x = stressful, y=mean, ymin=mean-se, ymax=mean+se))+
  geom_point(aes(x = stressful, y=mean, fill= mean > 0), size = 6, show.legend = FALSE, pch=21)+
  scale_fill_viridis_d()+
  coord_flip()+
  facet_wrap(~experiment_type)

e / a / b / f / g



metrics_24may <- count(full, resp_def, experiment_type) %>%
  as.data.frame()
write_csv(metrics_24may, "metrics_24may.csv")




#sd of the effect sizes
#cor plot between the different moderators 
# tpc could be energetics ? could be biochemistry?
# a tpc is fundamentally different 
# run metric separately from the full model?
# mean by study -- 
# run different random effects to play around with study and phylogenetic 
# life history traits (size, development, etc) vs state variables (size)
# rates?
#physiology vs. life history vs morphology

df_sam <- full %>%
  select(resp_def, experiment_type) %>%
  unique()
write_csv(df_sam, "response_metric.csv")
