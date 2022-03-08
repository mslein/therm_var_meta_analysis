##############################################################
## EFFECT SIZE CALCULATIONS AND DATA FORMATTING ##############
##############################################################
##loading necesssary libraries 
pacman::p_load(metafor, tidyverse)
#load data
dat_acclim_0<- read_csv("acclimation_dummytherm_30jun.csv")
dat_full_var_0 <- read_csv("metafor.csv")

##wrangling calculations for acclimation data
#make standardized SD columns based on N and SE vs SD
dat_acclim <-dat_acclim_0 %>% 
  mutate ( SD_constant= if_else ( variance_type== 0, constant_variance * sqrt(constant_samp), constant_variance))  %>%
  mutate ( SD_variable= if_else ( variance_type== 0, flux_variance * sqrt(flux_variance) , flux_variance)) 

#calculate ES for data
dat_acclim_ES <-escalc(measure="SMD", m1i=flux_resp, m2i=constant_resp, 
                       sd1i=`SD_variable`, sd2i= `SD_constant`, n1i=flux_samp, n2i=constant_samp, 
                       data=dat_acclim) 

# Note in the above, you loose ~15% of data because of missing sample sizes and SEs reported, below these are removed for now
dat_acclim_ES<-dat_acclim_ES %>% 
  filter (!is.na(yi)) %>%
  mutate(tpc_zone = case_when(resp_def %in% c("Performance breadth (Tbr)",
                                            "Scope of thermal tolerance (CTmax - CTmin)",
                                            "Thermal preference", 
                                            "Topt", "Tpref", "Speed", "LDH activity", "CCO activity", 
                                            "CS activity","log CO2 production",
                                            "Oxygen consumption", 
                                            "RMR") ~ "neutral", 
         resp_def %in% c("Time to heat knockdown", "Maximum burst swim speed", 
                         "maximal swimming speed",
                         "Ucrit", "Umax", "CT max", "CTmax", "LTmax", "ULT50") ~ "high", 
         resp_def %in% c( "CTmin", "LTmin") ~ "low"), 
         metric_type = case_when(resp_units %in% c("C") ~ "trait", 
                                 TRUE ~ "rate"))


write_csv(dat_acclim_ES, file = 'acclimation.csv')

##wrangling calculations for full variability data
#make standardized SD columns based on N and SE vs SD
dat_full_var <-dat_full_var_0 %>% 
  mutate ( SD_constant= if_else ( variance_type== 0, constant_variance * sqrt(constant_samp), constant_variance))  %>%
  mutate ( SD_variable= if_else ( variance_type== 0, flux_variance * sqrt(flux_variance) , flux_variance)) 

#calculate ES for data
dat_full_var_ES <-escalc(measure="SMD", m1i=flux_resp, m2i=constant_resp, 
                         sd1i=`SD_variable`, sd2i= `SD_constant`, n1i=flux_samp, n2i=constant_samp, 
                         data=dat_full_var, slab=paste(study_id, experiment_id, response_id, sep=", ")) 

# Note in the above, you loose ~15% of data because of missing sample sizes and SEs reported, below these are removed for now
dat_full_var_ES<-dat_full_var_ES %>% 
  filter (!is.na(yi)) %>%
  mutate(metric = case_when(resp_def %in% c("abdomen length", "adult weight", 
                                            "blade length/total leaf length", 
                                            "body length", "body mass", "body size", 
                                            "body weight", "carapace height", 
                                            "carapace width", 
                                            "final length", "final mass",
                                            "fore-limb length", 
                                            "group mean body weight", 
                                            "head length", "head width", 
                                            "hind-limb length", 
                                            "hypocotyl length", "mass", 
                                            "maximal length", 
                                            "snout-vent length", "tail length",
                                            "total length", 
                                            "wing centroid", "body (centroid) size", 
                                            "egg mass", "ovary mass, dry", 
                                            "testes mass, dry", 
                                            "pupal shell weight", "pupal weight") ~ "size", 
                            resp_def %in%  c("average cumulative number of eggs laid per female",
                                             "hatching success", 
                                             "offspring per mating",
                                             "total offspring", 
                                             "percent females") ~ "fecundity",
                            resp_def %in% c("Adult eclosion to death", 
                                            "dessication tolerance", 
                                            "egg to adult viability",
                                            "infestation rate", 
                                            "Mortaility", 
                                            "Parasitism success", 
                                            "startvation tolerance", 
                                            "success of parasitism", 
                                            "survival", "Survival", 
                                            "Larval mortaility", "longevity") ~ "survivorship", 
                            resp_def %in% c("aquatic speed", 
                                            "daily energy expenditure", 
                                            "distance covered", 
                                            "Energy consumption", 
                                            "feeding efficiency", 
                                            "feeding rate",
                                            "RMR", "sprint speed",
                                            "Sqrt Cellular energy allocation (CEA)", 
                                            "Sqrt of available energy", 
                                            "population growth rate","productivity", 
                                            "rate of change", 
                                            "specific growth rate") ~ "energetics", 
                            resp_def %in% c("development time", 
                                            "development to stages 35-37", 
                                            "days to first slough", 
                                            "developmental time", 
                                            "egg to adult", 
                                            "egg to pupa", "incubation period", 
                                            "median days to bolt", 
                                            "Third instar to adult", 
                                            "Third instar to pupa") ~ "development", 
                            resp_def %in% c("Catalase", "cortisol",
                                            "Hexokinase", "Hsp70", 
                                            "oxidative damage", 
                                            "Pyruvate kinase", "TAC", 
                                            "SOD", "glucose") ~ "biochemistry"))

write_csv(dat_full_var_ES, file = 'acute_cleaned.csv')


