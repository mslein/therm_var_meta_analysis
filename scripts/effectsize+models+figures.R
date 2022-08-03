#####################################################
## Manuscript coding, models, and figures ###########
#####################################################
##loading libraries 
pacman::p_load(tidyverse, metafor, patchwork, rotl, MuMIn, viridis, ape, wesanderson, glmulti, MuMIn)
#loading in the raw data file
dat_0<- read_csv("data/rawdata29jun22.csv")
##wrangling calculations for acclimation data
#make standardized SD columns based on N and SE vs SD
dat <-dat_0 %>% 
  mutate(SD_constant = case_when(variance_type %in% c(0) ~ constant_variance * sqrt(constant_samp), 
                                 variance_type %in% c(1) ~ constant_variance, 
                                 TRUE ~ (constant_variance/3.92)/(sqrt(constant_samp))), 
         SD_flux = case_when(variance_type %in% c(0) ~ flux_variance * sqrt(flux_samp), 
                             variance_type %in% c(1) ~ flux_variance, 
                             TRUE ~ (flux_variance/3.92)/(sqrt(flux_samp))))
dat_ES <-escalc(measure="SMD", m1i=flux_resp, m2i=constant_resp, 
                       sd1i=`SD_flux`, sd2i= `SD_constant`, n1i=flux_samp, n2i=constant_samp, 
                       data=dat) 
dat_ES2<-dat_ES %>% 
  filter (!is.na(yi)) %>%
  unite(original_name, c(genus, species), sep = " ", remove = FALSE) %>%
  mutate(original_name = tolower(original_name))
#selecting only organismal level responses, adult, and juvenile organisms, and species that did have phylogenies
dat_ES_final<- dat_ES2 %>%
  filter(org_level %in% c("0")) %>%
  filter(size %in% c("1", "2")) %>%
  filter(!original_name %in% c("phanaeus vindex", 
                              "cryptopygus sp.", 
                              "metacatharsius opacus")) %>%
  mutate(ecosystem = case_when(environment == "terrestrial" ~ "terrestrial", 
         TRUE ~ "aquatic/marine"))
#getting the phylogenies
taxa <- tnrs_match_names(as.character(unique(dat_ES_final$original_name)))
taxa_df <- taxa %>% 
  as.data.frame() %>% 
  select(search_string, unique_name) %>% 
  rename(original_name=search_string)
#creating the tree
tr <- tol_induced_subtree(ott_id(taxa)[is_in_tree(ott_id(taxa))]) #80 species
otl_tips <- strip_ott_ids(tr$tip.label, remove_underscores = FALSE) #remove ott_id's from tips
tr$tip.label <- otl_tips
tree <- compute.brlen(tr)
cor <- vcv(tree, cor = T) 
save(cor, file="phylo_cor.rda") #trying to save a stable version of this
#adding the correlations back into the dataframe for the model
dat_ES_final_2 <- left_join(dat_ES_final, taxa_df, by="original_name") %>%
  separate(unique_name, c('u_genus', 'u_species')) %>%
  unite(phylo,c("u_genus", "u_species"), sep="_") %>%
  filter(yi <650) %>%
  mutate(duration_standard = case_when(duration_units %in% c("weeks") ~ duration*7, 
                                       duration_units %in% c("hours") ~ duration/24,
                                       TRUE ~ duration))
####candidate models
#full model
full_mod<-rma.mv(yi, vi, 
            mods = ~ +I(flux_range-mean(flux_range)) + I(mean_temp_constant-mean(mean_temp_constant))
            +experiment_type
            +I(secondary_temp - mean(secondary_temp))
            +duration_standard +as.factor(exp_age) 
            +as.factor(size) +ecosystem +trait_directionality +vi,
            random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
            R = list(phylo=cor),
            method="REML", data=dat_ES_final_2)
#without experiment type
mod1 <- rma.mv(yi, vi, 
               mods = ~ +I(flux_range-mean(flux_range)) + I(mean_temp_constant-mean(mean_temp_constant))
               +I(secondary_temp - mean(secondary_temp))
               +duration_standard +as.factor(exp_age) 
               +as.factor(size) +ecosystem +trait_directionality+vi,
               random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
               R = list(phylo=cor),
               method="REML", data=dat_ES_final_2)
#without mean temp
mod2 <- rma.mv(yi, vi, 
               mods = ~ +I(flux_range-mean(flux_range)) 
               +experiment_type
               +I(secondary_temp - mean(secondary_temp))
               +duration_standard +as.factor(exp_age) 
               +as.factor(size) +ecosystem +trait_directionality+vi,
               random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
               R = list(phylo=cor),
               method="REML", data=dat_ES_final_2)
#without flux range
mod3 <- rma.mv(yi, vi, 
               mods = ~ + I(mean_temp_constant-mean(mean_temp_constant))
               +experiment_type
               +I(secondary_temp - mean(secondary_temp))
               +duration_standard +as.factor(exp_age) 
               +as.factor(size) +ecosystem +trait_directionality+vi,
               random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
               R = list(phylo=cor),
               method="REML", data=dat_ES_final_2)
#without both flux range + mean temp
mod4 <- rma.mv(yi, vi, 
               mods = ~+experiment_type
               +I(secondary_temp - mean(secondary_temp))
               +duration_standard +as.factor(exp_age) 
               +as.factor(size) +ecosystem +trait_directionality+vi,
               random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
               R = list(phylo=cor),
               method="REML", data=dat_ES_final_2)
#without age
mod9 <- rma.mv(yi, vi, 
               mods = ~ +I(flux_range-mean(flux_range)) + I(mean_temp_constant-mean(mean_temp_constant))
               +experiment_type
               +I(secondary_temp - mean(secondary_temp))
               +duration_standard
               +as.factor(size) +ecosystem +trait_directionality+vi,
               random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
               R = list(phylo=cor),
               method="REML", data=dat_ES_final_2)
#without size
mod10 <- rma.mv(yi, vi, 
                mods = ~ +I(flux_range-mean(flux_range)) + I(mean_temp_constant-mean(mean_temp_constant))
                +experiment_type
                +I(secondary_temp - mean(secondary_temp))
                +duration_standard +as.factor(exp_age) 
                +ecosystem +trait_directionality+vi,
                random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
                R = list(phylo=cor),
                method="REML", data=dat_ES_final_2)
#without ecosystem 
mod11 <- rma.mv(yi, vi, 
                mods = ~ +I(flux_range-mean(flux_range)) + I(mean_temp_constant-mean(mean_temp_constant))
                +experiment_type
                +I(secondary_temp - mean(secondary_temp))
                +duration_standard +as.factor(exp_age) 
                +as.factor(size) +trait_directionality+vi,
                random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
                R = list(phylo=cor),
                method="REML", data=dat_ES_final_2)
#without age and size
mod12 <-rma.mv(yi, vi, 
               mods = ~ +I(flux_range-mean(flux_range)) + I(mean_temp_constant-mean(mean_temp_constant))
               +experiment_type
               +I(secondary_temp - mean(secondary_temp))
               +duration_standard +trait_directionality +ecosystem+vi,
               random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
               R = list(phylo=cor),
               method="REML", data=dat_ES_final_2)
#without size, age, and ecosystem
mod13 <-rma.mv(yi, vi, 
               mods = ~ +I(flux_range-mean(flux_range)) + I(mean_temp_constant-mean(mean_temp_constant))
               +experiment_type
               +I(secondary_temp - mean(secondary_temp))
               +duration_standard +trait_directionality+vi,
               random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
               R = list(phylo=cor),
               method="REML", data=dat_ES_final_2)
#########################
#########Figures#########
#########################

#Figure 3
p0 <-dat_ES_final%>%
  filter(yi <80) %>%
  group_by(experiment_type) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=experiment_type, y=mean, ymax=mean+se, ymin=mean-se))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=7, position=position_dodge(width = .75), colour="#440154")+
  geom_linerange(position=position_dodge(width = .75), size=1,  colour="#440154")+
  coord_flip()+
  theme_bw()+
  xlab("Experiment type")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

p1 <-dat_ES_final%>%
  filter(yi <80) %>%
  group_by(exp_age) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=as.factor(exp_age), y=mean, ymax=mean+se, ymin=mean-se))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=7, position=position_dodge(width = .75), colour="#3b528b")+
  geom_linerange(position=position_dodge(width = .75), size=1, colour="#3b528b")+
  coord_flip()+
  theme_bw()+
  xlab("Experimental age")+
  ylab("")+
  scale_x_discrete(labels = c("larval",'juvenile','adult'))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

p2 <-dat_ES_final%>%
  filter(yi <80) %>%
  group_by(size) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=as.factor(size), y=mean, ymax=mean+se, ymin=mean-se))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=7, position=position_dodge(width = .75), colour= "#21918c")+
  geom_linerange(position=position_dodge(width = .75), size=1, colour= "#21918c")+
  coord_flip()+
  theme_bw()+
  xlab("Size")+
  ylab("")+
  scale_x_discrete(labels = c("small",'medium'))+
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20,face="bold"))


p3 <-dat_ES_final%>%
  filter(yi <80) %>%
  group_by(ecosystem) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=ecosystem, y=mean, ymax=mean+se, ymin=mean-se))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=7, position=position_dodge(width = .75), colour= "#5ec962")+
  geom_linerange(position=position_dodge(width = .75), size=1, colour= "#5ec962")+
  coord_flip()+
  scale_colour_manual(values = wes_palette("Cavalcanti1"))+
  theme_bw()+
  xlab("Ecosystem")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

p4 <-dat_ES_final%>%
  filter(yi <80) %>%
  group_by(trait_directionality) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=trait_directionality, y=mean, ymax=mean+se, ymin=mean-se))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=7, position=position_dodge(width = .75), colour="#fde725")+
  geom_linerange(position=position_dodge(width = .75), size=1, colour="#fde725")+
  coord_flip()+
  scale_colour_manual(values = wes_palette("Cavalcanti1"))+
  theme_bw()+
  xlab("Trait directionality")+
  ylab("Effect size")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

p0 / p1 / p2 / p3 / p4

#Figure 4
mean_temp<- dat_ES_final_2 %>%
  filter(between(yi, -10, 10)) %>%
  ggplot(aes(x=mean_temp_constant, y=yi, colour=mean_temp_constant))+
  geom_point(alpha=0.5, size=2)+
  scale_colour_viridis(option="D")+
  #geom_smooth(method="lm")+
  theme_bw()+
  xlab("Mean temperature (°C)")+
  ylab("Effect size")+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))+
  geom_abline(intercept = 0.7056, slope = 0.0139, color="black", size=1)


flux_range <- dat_ES_final_2 %>%
  filter(between(yi, -10, 10)) %>%
  ggplot(aes(x=flux_range, y=yi, colour=flux_range))+
  geom_point(alpha=0.5, size=2)+
  scale_colour_viridis(option="D")+
  #geom_smooth(method="lm")+
  theme_bw()+
  xlab("Fluctuation magnitude (°C)")+
  ylab("Effect size")+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))+
  geom_abline(intercept = 0.7056, slope = -0.0145, color="black", size=1)

mean_temp + flux_range

#Figure 5
duration<- dat_ES_final_2 %>%
  filter(between(yi, -10, 10)) %>%
  #filter(between(duration, 0, 200)) %>%
  ggplot(aes(x=duration_standard, y=yi, colour=duration_standard))+
  geom_point(alpha=0.5, size=2)+
  scale_colour_viridis(option="D")+
  theme_bw()+
  xlab("Time (days)")+
  ylab("Effect size")+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))+
  geom_abline(intercept = 0.7056, slope = -0.0164 , color="black", size=1)

#######################
########SI Figures#####
#######################


#model comparison
aic_table <- AICc(full_mod, mod1, mod2, mod3, mod4, mod9, mod10, mod11, mod12, mod13) %>%
  as.data.frame() %>%
  mutate(delta_aicc = AICc- 118404.9) %>%
  arrange(desc(delta_aicc))

#funnel plot
normalized <- dat_ES_final_2 %>%
  filter(between(yi, -10, 10))
dummy_full_mod <- rma.mv(yi, vi,
                         random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
                         R = list(phylo=cor),
                         method="REML", data=dat_ES_final_2)
funnel(dummy_full_mod)

##histogram
dat_ES_final_2 %>%
  filter(between(yi, -10, 10)) %>%
  ggplot(aes(x=yi))+
  geom_histogram(bins=100)+
  theme_bw()+
  xlab("SMD")

##I2 calculations from best models 
#phylo level heterogeneity
full_mod$sigma2[1]/sum(full_mod$sigma2[1], full_mod$sigma2[2], full_mod$sigma2[3]) #0.9309465
#study level heterogeneity
full_mod$sigma2[2]/sum(full_mod$sigma2[1], full_mod$sigma2[2], full_mod$sigma2[3]) #0.06874321
#response heterogeneity 
full_mod$sigma2[3]/sum(full_mod$sigma2[1], full_mod$sigma2[2], full_mod$sigma2[3]) #0.0003103212





