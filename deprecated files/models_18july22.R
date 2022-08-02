pacman::p_load(tidyverse, metafor, patchwork, rotl, MuMIn, viridis, ape, wesanderson)

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

# Note in the above, you loose ~11% of data because of Some 'yi' and/or 'vi' values equal to +-Inf. Recoded to NAs
# below these are removed for now

missingspec<-read_csv("missing_species.csv")

dat_ES2<-dat_ES %>% 
  filter (!is.na(yi)) %>%
  unite(original_name, c(genus, species), sep = " ", remove = FALSE) %>%
  mutate(original_name = tolower(original_name))

dat_ES_final<- dat_ES2 %>%
  filter(org_level %in% c("0")) %>%
  filter(size %in% c("1", "2")) %>%
  filter(!original_name %in% c("phanaeus vindex", 
                              "cryptopygus sp.", 
                              "metacatharsius opacus")) %>%
  mutate(ecosystem = case_when(environment == "terrestrial" ~ "terrestrial", 
         TRUE ~ "aquatic/marine"))

#83 species
length(unique(dat_ES_final$original_name))
  #missing 3 species
  #phanaeus vindex
  #cryptopygus sp.
  #metacatharsius opacus
#about 37 effect sizes lost in this (looks like these aren't in the tree of life?)
dat_ES_final %>%
  filter(original_name %in% c("phanaeus vindex", 
                              "cryptopygus sp.", 
                              "metacatharsius opacus")) %>%
  length()

taxa <- tnrs_match_names(as.character(unique(dat_ES_final$original_name)))
taxa_df <- taxa %>% 
  as.data.frame() %>% 
  select(search_string, unique_name) %>% 
  rename(original_name=search_string)

tr <- tol_induced_subtree(ott_id(taxa)[is_in_tree(ott_id(taxa))]) #80 species
otl_tips <- strip_ott_ids(tr$tip.label, remove_underscores = FALSE) #remove ott_id's from tips
tr$tip.label <- otl_tips
plot(tr, show.tip.label = TRUE)


tree <- compute.brlen(tr)
cor <- vcv(tree, cor = T) 
save(cor, file="phylo_cor.rda")
cor_stable <- 

dat_ES_final_2 <- left_join(dat_ES_final, taxa_df, by="original_name") %>%
  separate(unique_name, c('u_genus', 'u_species')) %>%
  unite(phylo,c("u_genus", "u_species"), sep="_") %>%
  filter(yi <650) %>%
  mutate(plus_yi = yi + 1)


mod<-rma.mv(yi, vi, 
            mods = ~ +I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
            +I(mean_temp_constant-mean(mean_temp_constant)) * experiment_type
            +I(flux_range-mean(flux_range)) * experiment_type
            +I(secondary_temp - mean(secondary_temp))*experiment_type 
            +experiment_type*duration +experiment_type*as.factor(exp_age) 
            +experiment_type*as.factor(size) +experiment_type*ecosystem +experiment_type*trait_directionality,
            random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
            R = list(phylo=cor),
            method="REML", data=dat_ES_final_2)

summary(mod)

mod_exp <- rma.mv(yi, vi, 
                  mods = ~experiment_type,
                  random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor),
                  method="REML", data=dat_ES_final_2)

mod_temp <- rma.mv(yi, vi, 
                   mods = ~ +I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                   +I(secondary_temp - mean(secondary_temp)) + experiment_type,
                   random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
                   R = list(phylo=cor),
                   method="REML", data=dat_ES_final_2)
mod_eco <- rma.mv(yi, vi, 
                  mods = ~ +I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                  +I(secondary_temp - mean(secondary_temp)) + experiment_type + ecosystem,
                  random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor),
                  method="REML", data=dat_ES_final_2)

AIC(mod, mod_exp, mod_temp, mod_eco)


1) acc/acute only, 2) 1 + temp stuff, 3) 2 + metrics, 4) 3+ ecosystems, 5) everything. Or something along those lines





dat_ES_scale <- dat_ES2 %>%
  filter(!experiment_type %in% c("acclimation")) %>%
  filter(!original_name %in% c("phanaeus vindex", 
                               "cryptopygus sp.", 
                               "metacatharsius opacus", 
                               "litopenaeus vannamei",
                               "enterococcus sp.", 
                               "escherichia coli")) %>%
  mutate(ecosystem = case_when(environment == "terrestrial" ~ "terrestrial", 
                               TRUE ~ "aquatic/marine"))

taxa2 <- tnrs_match_names(as.character(unique(dat_ES_scale$original_name)))
taxa_df2 <- taxa2 %>% 
  as.data.frame() %>% 
  select(search_string, unique_name) %>% 
  rename(original_name=search_string)


tr2 <- tol_induced_subtree(ott_id(taxa2)[is_in_tree(ott_id(taxa2))]) #80 species
otl_tips2 <- strip_ott_ids(tr2$tip.label, remove_underscores = FALSE) #remove ott_id's from tips
tr2$tip.label <- otl_tips2
plot(tr2, show.tip.label = TRUE)


tree2 <- compute.brlen(tr2)
cor2 <- vcv(tree2, cor = T) 


dat_ES_scale_2 <- left_join(dat_ES_scale, taxa_df2, by="original_name") %>%
  separate(unique_name, c('u_genus', 'u_species')) %>%
  unite(phylo,c("u_genus", "u_species"), sep="_") %>%
  filter(yi <650) %>%
  mutate(plus_yi = yi + 1) %>%
  filter(!original_name %in% c("litopenaeus vannamei"))

mod2<-rma.mv(yi, vi, 
            mods = ~I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
            +I(secondary_temp - mean(secondary_temp))
            +duration +as.factor(exp_age) 
            +as.factor(size) +ecosystem +as.factor(org_level),
            random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
            R = list(phylo=cor2),
            method="REML", data=dat_ES_scale_2)

summary(mod2)

anti_join(dat_ES_scale, taxa_df, by="original_name")


count(dat_ES_scale, experiment_type, org_level)

# there is only organismal level resps in acclim
#should i filter this out?
count(dat_ES2, experiment_type, org_level)

#only size 1 and 2 have overlap between acute + acclim
count(dat_ES2, experiment_type, size)

#binning effect sizes between -300 and 300
dat_ES_final %>%
  filter(between(yi, -300, 300)) %>%
  ggplot(aes(x=yi))+
  geom_histogram()
  
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
  ggplot(aes(x=experiment_type, y=mean, ymax=mean+se, ymin=mean-se, shape=experiment_type, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=7, position=position_dodge(width = .75))+
  geom_linerange(position=position_dodge(width = .75), size=1)+
  coord_flip()+
  scale_colour_manual(values = wes_palette("Cavalcanti1"))+
  theme_bw()+
  xlab("Experiment type")

p1 <-dat_ES_final%>%
  filter(yi <80) %>%
  group_by(experiment_type, exp_age) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=as.factor(exp_age), y=mean, ymax=mean+se, ymin=mean-se, shape=experiment_type, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=7, position=position_dodge(width = .75))+
  geom_linerange(position=position_dodge(width = .75), size=1)+
  coord_flip()+
  scale_colour_manual(values = wes_palette("Cavalcanti1"))+
  theme_bw()+
  xlab("Experimental age")+
  scale_x_discrete(labels = c("Larval",'Juvenile','Adult'))

p2 <-dat_ES_final%>%
  filter(yi <80) %>%
  group_by(experiment_type, size) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=as.factor(size), y=mean, ymax=mean+se, ymin=mean-se, shape=experiment_type, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=7, position=position_dodge(width = .75))+
  geom_linerange(position=position_dodge(width = .75), size=1)+
  coord_flip()+
  scale_colour_manual(values = wes_palette("Cavalcanti1"))+
  theme_bw()+
  xlab("Size")+
  scale_x_discrete(labels = c("Small",'Medium'))

p3 <-dat_ES_final%>%
  filter(yi <80) %>%
  group_by(experiment_type, stressful) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=stressful, y=mean, ymax=mean+se, ymin=mean-se, shape=experiment_type, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=7, position=position_dodge(width = .75))+
  geom_linerange(position=position_dodge(width = .75), size=1)+
  coord_flip()+
  scale_colour_manual(values = wes_palette("Cavalcanti1"))+
  theme_bw()+
  xlab("Stressful")

p4 <-dat_ES_final%>%
  filter(yi <80) %>%
  group_by(experiment_type, environment) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=environment, y=mean, ymax=mean+se, ymin=mean-se, shape=experiment_type, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=7, position=position_dodge(width = .75))+
  geom_linerange(position=position_dodge(width = .75), size=1)+
  coord_flip()+
  scale_colour_manual(values = wes_palette("Cavalcanti1"))+
  theme_bw()+
  xlab("Environment")

p5 <-dat_ES_final%>%
  filter(yi <80) %>%
  group_by(experiment_type, trait_directionality) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=trait_directionality, y=mean, ymax=mean+se, ymin=mean-se, shape=experiment_type, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=7, position=position_dodge(width = .75))+
  geom_linerange(position=position_dodge(width = .75), size=1)+
  coord_flip()+
  scale_colour_manual(values = wes_palette("Cavalcanti1"))+
  theme_bw()+
  xlab("Trait directionality")

p0 / p1 / p2 / p4 / p5


dat_ES_final_2 %>%
  filter(experiment_type == "acclimation") %>%
  filter(resp_def %in% c("upper lethal limit", 
                         "ULT50", "Tpref", "Topt", "Thermal preference",
                         "Scope of thermal tolerance (CTmax-Ctmin)", 
                         "Performance breadth (Tbr)", 
                         "LTmin", "LTmax", "LT50", 
                         "lethal thermal limit", "CTmin", 
                         "CTmax", "Ctmax", "CT max")) %>%
  mutate(rep_def_broad = case_when(resp_def %in% c("upper lethal limit", 
                                                   "ULT50", "LTmax", "LT50",
                                                   "CTmax", "Ctmax", "CT max") ~ "upper", 
                                   resp_def %in% c("Tpref", "Topt", "Thermal preference",
                                   "Scope of thermal tolerance (CTmax-Ctmin)", 
                                   "Performance breadth (Tbr)") ~ "neutral",
                                    TRUE ~ "lower")) %>%
  group_by(rep_def_broad) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=rep_def_broad, y=mean, ymax=mean+se, ymin=mean-se))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=4, position=position_dodge(width = .75))+
  geom_linerange(position=position_dodge(width = .75), size=1)+
  coord_flip()
  


#size
b <- dat_ES_final %>%
  filter(yi <80) %>%
  group_by(experiment_type, size)%>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=as.factor(size), y=mean))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_linerange(aes(x = as.factor(size), y=mean, ymin=mean-se, ymax=mean+se))+
  geom_point(aes(x = size, y=mean, fill= mean > 0), size = 6, show.legend = FALSE, pch=21)+
  scale_fill_viridis_d()+
  coord_flip()+
  facet_wrap(~experiment_type)
#stressful
c <- dat_ES_final %>%
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

#experiment type
e <- dat_ES_final %>%
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
f <-dat_ES_final %>%
  filter(yi <80) %>%
  group_by(experiment_type, ecosystem) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=ecosystem, y=mean))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_linerange(aes(x = ecosystem, y=mean, ymin=mean-se, ymax=mean+se))+
  geom_point(aes(x = ecosystem, y=mean, fill= mean > 0), size = 6, show.legend = FALSE, pch=21)+
  scale_fill_viridis_d()+
  coord_flip()+
  facet_wrap(~experiment_type)

e / a / b / f 

dat_ES_final_2 %>%
  filter(yi <80) %>%
  ggplot(aes(x=as.factor(exp_age), y=log(plus_yi), colour=experiment_type))+
  geom_jitter(width = 0.02, alpha=0.5)+
  coord_flip()+
  facet_wrap(~experiment_type)+
  scale_colour_viridis(option="magma", discrete=TRUE)




dat_ES_final_2 %>%
  ggplot(aes(x=mean_temp_constant, y=log(plus_yi), colour=log(plus_yi)))+
  geom_point()+
  scale_colour_viridis(discrete=FALSE, option= "D")+
  facet_wrap(~experiment_type)+
  geom_smooth(method="lm")

dat_ES_final_2 %>%
  ggplot(aes(x=flux_range, y=log(plus_yi), colour=log(plus_yi)))+
  geom_point(size=2)+
  scale_colour_viridis(discrete=FALSE, option= "A")


#removes 300 rows of data i think because of the log transform
dat_ES_final_2 %>%
  filter(between(yi, -50, 50)) %>%
  filter(flux_range %in% c(3,4,5,6,7,8,10,12,14,15,20)) %>%
  ggplot(aes(colour=flux_range, x=mean_temp_constant, y=yi, shape=experiment_type))+
  geom_point(alpha=0.5, size=2)+
  scale_colour_viridis(discrete=FALSE, option= "D")+
  facet_wrap(~flux_range)+
  geom_smooth(method="lm")+
  theme_bw()

dat_ES_final_2 %>%
  filter(between(yi, -50, 50)) %>%
  filter(flux_range %in% c(4,6,8,10,12,14)) %>%
  ggplot(aes(colour=flux_range, x=mean_temp_constant, y=yi), shape=experiment_type)+
  geom_point(alpha=0.5, size=2)+
  scale_colour_viridis(discrete=FALSE, option= "D")+
  facet_wrap(~flux_range)+
  geom_smooth(method="lm")+
  theme_bw()

dat_ES_final_2 %>%
  filter(between(yi, -50, 50)) %>%
  filter(flux_range %in% c(4,6,8,10,12,14)) %>%
  ggplot(aes(colour=experiment_type, x=mean_temp_constant, y=yi, shape=experiment_type))+
  geom_point(alpha=0.5, size=2)+
  scale_colour_manual(values = wes_palette("Cavalcanti1"))+
  facet_wrap(~flux_range)+
  geom_smooth(method="lm")+
  theme_bw()

dat_ES_final_2 %>%
  filter(between(yi, -50, 50)) %>%
  filter(flux_range %in% c(4,6,8,10,12,14)) %>%
  ggplot(aes(colour=experiment_type, x=mean_temp_constant, y=yi, shape=experiment_type))+
  geom_point(alpha=0.5, size=2)+
  scale_colour_manual(values = wes_palette("Cavalcanti1"))+
  facet_grid(rows = vars(experiment_type), cols = vars(flux_range))+
  geom_smooth(method="lm")+
  theme_bw()


#secondary temp figure 
dat_ES_final_2 %>%
  filter(between(yi, -50, 50)) %>%
  filter(flux_range %in% c(4,6,8,10,12,14)) %>%
  ggplot(aes(colour=experiment_type, x=secondary_temp, y=yi, shape=experiment_type))+
  geom_point(alpha=0.5, size=2)+
  scale_colour_manual(values = wes_palette("Cavalcanti1"))+
  facet_wrap(~experiment_type)+
  geom_smooth(method="lm")+
  theme_bw()

subset<- dat_ES_final_2 %>%
  filter(between(yi, -50, 50))

  




