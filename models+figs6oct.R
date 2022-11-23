#####################################################
## Manuscript coding, models, and figures ###########
#####################################################


##loading libraries 
pacman::p_load(tidyverse, metafor, patchwork, rotl, MuMIn, viridis, ape, wesanderson, glmulti, MuMIn, metaAidR, corrplot)
#loading in the raw data file
dat_0<- read_csv("data/rawdata29jun22.csv") %>%
  select(-add_covariate, -add_covariate_value, -extra_covariate, -extra_value, -larger_group, -flux_pattern) 
##wrangling calculations for acclimation data
#make standardized SD columns based on N and SE vs SD
dat <-dat_0 %>% 
  mutate(SD_constant = case_when(variance_type %in% c(0) ~ constant_variance * sqrt(constant_samp), # sd = se / sqrt(n)
                                 variance_type %in% c(1) ~ constant_variance, # sd= sd
                                 TRUE ~ sqrt(flux_samp)*(flux_variance)/3.92), # sd = sqrt(n) * ci / 3.92
         SD_flux = case_when(variance_type %in% c(0) ~ flux_variance * sqrt(flux_samp), 
                             variance_type %in% c(1) ~ flux_variance, 
                             TRUE ~ sqrt(flux_samp)*(flux_variance)/3.92), 
         SE_constant = case_when(variance_type %in% c(0) ~ constant_variance, 
                                 variance_type %in% c(1) ~ constant_variance / sqrt(constant_samp), 
                                 TRUE ~ constant_variance/3.92), 
         SE_flux = case_when(variance_type %in% c(0) ~ flux_variance , 
                             variance_type %in% c(1) ~ flux_variance / sqrt(flux_samp), 
                             TRUE ~ flux_variance/3.92))
#calculating the effect size--> standardized mean different (SMD)
dat_ES <-escalc(measure="SMD", m1i=flux_resp, m2i=constant_resp, 
                sd1i=`SD_flux`, sd2i= `SD_constant`, n1i=flux_samp, n2i=constant_samp, 
                data=dat) 
#tidying the data to remove an NA's generating from the effect size calculation
dat_ES2<-dat_ES %>% 
  drop_na() %>% #4% data loss
  unite(original_name, c(genus, species), sep = " ", remove = FALSE) %>%
  mutate(original_name = tolower(original_name))
#selecting only organismal level responses and species that did have phylogenies, 
#excluding the one endotherm that made into the data by accident
#transforming environment into ecosystem to combine aquatic + marine studies
#transforming size into 3 categories, so that the XS organisms combined with small while XL combined with large
dat_ES_final<- dat_ES2 %>% #kramarz seems to be quite the outlier in the qqplot
  filter(org_level %in% c("0")) %>%
  filter(!original_name %in% c("phanaeus vindex", 
                               "cryptopygus sp.", 
                               "metacatharsius opacus", 
                               "coturnix japonica")) %>%
  mutate(ecosystem = case_when(environment == "terrestrial" ~ "terrestrial", 
                               TRUE ~ "aquatic/marine")) %>%
  mutate(size =case_when(size %in% c("0","1") ~ "1", 
                         size %in% c("2") ~ "2", 
                         size %in% c("3", "4") ~ "3")) %>%
  mutate(experiment_type = as.factor(experiment_type), 
         exp_age=as.factor(exp_age), 
         size=as.factor(size), 
         trait_directionality=as.factor(trait_directionality), 
         ecosystem=as.factor(ecosystem))

#specifying the reference level for each dummy variable, with the level with the most observations
dat_ES_final$experiment_type <- relevel(dat_ES_final$experiment_type, "acute")
dat_ES_final$exp_age <- relevel(dat_ES_final$exp_age, "1")
dat_ES_final$size <- relevel(dat_ES_final$size, "1")
dat_ES_final$trait_directionality <- relevel(dat_ES_final$trait_directionality, "positive")
dat_ES_final$ecosystem <- relevel(dat_ES_final$ecosystem, "terrestrial")

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
  mutate(duration_standard = case_when(duration_units %in% c("weeks") ~ duration*7, 
                                       duration_units %in% c("hours") ~ duration/24,
                                       TRUE ~ duration)) %>%
  mutate(sqrt_ni_1 = sqrt((flux_samp+constant_samp)/(flux_samp * constant_samp)), 
         ni_1 =(flux_samp+constant_samp)/(flux_samp * constant_samp)) 


#setting up neccesary packages to make the shared control matrix following Rodgers et al 2021
#install.packages("devtools")
#install_github("daniel1noble/metaAidR")
library(devtools)
# Error calc-residual variance at the observation level as metafor does not add this by default
dat_ES_final_2$obs <- 1:dim(dat_ES_final_2)[1]
# Make shared control matrix for ROM dataset
dat_ES_final_2$sc_cluster <- interaction(dat_ES_final_2$study_id, 
                                         dat_ES_final_2$shared_control)
# Create the Shared control V matrix
V <- metaAidR::make_VCV_matrix(dat_ES_final_2, "vi", cluster= "sc_cluster", 
                               type = "vcv", rho = 0.5)


#test models
int_mod <- rma.mv(yi=yi, V=vi, 
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)


##no random effects 
no_random_mod<- rma.uni(yi=yi, vi, 
                        mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                        +I(flux_range-mean(flux_range))*experiment_type
                        +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                        +duration_standard
                        +exp_age*experiment_type
                        +size*experiment_type
                        +ecosystem
                        +trait_directionality*experiment_type
                        +sqrt_ni_1,
                        method="ML", data=dat_ES_final_2)


####candidate models
#full model
full_mod<-rma.mv(yi=yi, V=V, 
                 mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                 +experiment_type
                 +I(secondary_temp - mean(secondary_temp))
                 +duration_standard
                 +exp_age
                 +size
                 +ecosystem
                 +trait_directionality
                 +sqrt_ni_1,
                 random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                 R = list(phylo=cor), test="t",
                 method="ML", data=dat_ES_final_2)
#without ecosystem
full_moda<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +experiment_type
                  +I(secondary_temp - mean(secondary_temp))
                  +duration_standard
                  +exp_age
                  +size
                  +trait_directionality
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

#without duration
full_modb<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +experiment_type
                  +I(secondary_temp - mean(secondary_temp))
                  +exp_age
                  +size
                  +ecosystem
                  +trait_directionality
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)
#without secondary temp
full_modc<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +experiment_type
                  +duration_standard
                  +exp_age
                  +size
                  +ecosystem
                  +trait_directionality
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)
#without publication bias
full_modd<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +experiment_type
                  +I(secondary_temp - mean(secondary_temp))
                  +duration_standard
                  +exp_age
                  +size
                  +ecosystem
                  +trait_directionality,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

model.sel(full_mod, full_moda, full_modb, full_modc, full_modd)
#full_mod, full_moda, full_modc are all within 2AIC, so i pulled the model with all the terms and built the subsequent models based on that 


# hypothesis 1
#full model with interaction between experiment_type and the mods
full_mode<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))*experiment_type
                  +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                  +I(secondary_temp - mean(secondary_temp))
                  +duration_standard
                  +exp_age*experiment_type
                  +size*experiment_type
                  +ecosystem
                  +trait_directionality*experiment_type
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor),test="t",
                  method="ML", data=dat_ES_final_2)

#no experiment type 
full_modf<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +duration_standard
                  +I(secondary_temp - mean(secondary_temp))
                  +exp_age
                  +size
                  +ecosystem
                  +trait_directionality
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

#with experiment type but no interaction 
full_modg<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +experiment_type
                  +I(secondary_temp - mean(secondary_temp))
                  +duration_standard
                  +exp_age
                  +size
                  +ecosystem
                  +trait_directionality
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

#no interaction between experiment type and flux range
full_modh<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))
                  +I(secondary_temp - mean(secondary_temp))
                  +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                  +duration_standard
                  +exp_age*experiment_type
                  +size*experiment_type
                  +ecosystem
                  +trait_directionality*experiment_type
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

#no interaction between experiment type and mean temp
full_modi<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))*experiment_type
                  +I(mean_temp_constant-mean(mean_temp_constant))
                  +duration_standard
                  +I(secondary_temp - mean(secondary_temp))
                  +exp_age*experiment_type
                  +size*experiment_type
                  +ecosystem
                  +trait_directionality*experiment_type
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

#no interaction between experiment type and age
full_modj<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))*experiment_type
                  +I(secondary_temp - mean(secondary_temp))
                  +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                  +duration_standard
                  +exp_age
                  +size*experiment_type
                  +ecosystem
                  +trait_directionality*experiment_type
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)
#no interaction between experiment type and size
full_modk<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))*experiment_type
                  +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                  +duration_standard
                  +I(secondary_temp - mean(secondary_temp))
                  +exp_age*experiment_type
                  +size
                  +ecosystem
                  +trait_directionality*experiment_type
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

#no interaction between experiment type and trait directionality
full_modl<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))*experiment_type
                  +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                  +duration_standard
                  +I(secondary_temp - mean(secondary_temp))
                  +exp_age*experiment_type
                  +size*experiment_type
                  +ecosystem
                  +trait_directionality
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

model.sel(full_mode, full_modf, full_modg, full_modh, full_modi, full_modj, full_modk, full_modl)
model.sel(full_mod, full_moda, full_modb, full_modc, full_modd)
#full mod e prevails by more than 2 AIC 

#h2 and h3 with interaction between flux range and experiment (based on full_mode as base)
# hypothesis 2 
# no mean temp
full_modaa<-rma.mv(yi=yi, V=V, 
                   mods = ~ I(flux_range-mean(flux_range))*experiment_type
                   +duration_standard
                   +I(secondary_temp - mean(secondary_temp))
                   +exp_age*experiment_type
                   +size*experiment_type
                   +ecosystem
                   +trait_directionality*experiment_type
                   +sqrt_ni_1,
                   random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                   R = list(phylo=cor), test="t",
                   method="ML", data=dat_ES_final_2)
#no flux range
full_modbb <- rma.mv(yi=yi, V=V, 
                     mods = ~ I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +duration_standard
                     +I(secondary_temp - mean(secondary_temp))
                     +exp_age*experiment_type
                     +size*experiment_type
                     +ecosystem
                     +trait_directionality*experiment_type
                     +sqrt_ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no interaction between flux range and mean temp
full_modcc <- rma.mv(yi=yi, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*experiment_type
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +duration_standard
                     +exp_age*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +size*experiment_type
                     +ecosystem
                     +trait_directionality*experiment_type
                     +sqrt_ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no flux range or mean temp
full_moddd <- rma.mv(yi=yi, V=V, 
                     mods = ~ duration_standard
                     +exp_age*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +size*experiment_type
                     +ecosystem
                     +trait_directionality*experiment_type
                     +sqrt_ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)

model.sel(full_mode, full_modh,
          full_modaa, full_modbb, full_modcc, full_moddd)
#full_mode prevails
model.sel(full_mode, full_modf, full_modg, full_modh, full_modi, full_modj, full_modk, full_modl)
model.sel(full_mod, full_moda, full_modb, full_modc, full_modd)

#hypothesis 3
#no interaction with age (based on full_mode)
full_modee <- rma.mv(yi=yi, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +duration_standard
                     +exp_age
                     +size*experiment_type
                     +ecosystem
                     +trait_directionality*experiment_type
                     +sqrt_ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no age
full_modff <- rma.mv(yi=yi, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +duration_standard
                     +size*experiment_type
                     +ecosystem
                     +trait_directionality*experiment_type
                     +sqrt_ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no interaction with size
full_modgg <- rma.mv(yi=yi, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +duration_standard
                     +exp_age*experiment_type
                     +size
                     +ecosystem
                     +trait_directionality*experiment_type
                     +sqrt_ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no size
full_modhh <- rma.mv(yi=yi, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +duration_standard
                     +exp_age*experiment_type
                     +ecosystem
                     +trait_directionality*experiment_type
                     +sqrt_ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)

#no interaction with trait 
full_modii <- rma.mv(yi=yi, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +duration_standard
                     +exp_age*experiment_type
                     +size*experiment_type
                     +ecosystem
                     +trait_directionality
                     +sqrt_ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no age or size
full_modjj <- rma.mv(yi=yi, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +duration_standard
                     +ecosystem
                     +trait_directionality*experiment_type
                     +sqrt_ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no interaction with age or size 
full_modkk <- rma.mv(yi=yi, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +duration_standard
                     +exp_age
                     +size
                     +ecosystem
                     +trait_directionality*experiment_type
                     +sqrt_ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no interaction with trait 
full_modll<- rma.mv(yi=yi, V=V, 
                    mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                    +I(flux_range-mean(flux_range))*experiment_type
                    +I(secondary_temp - mean(secondary_temp))
                    +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                    +duration_standard
                    +exp_age*experiment_type
                    +size*experiment_type
                    +ecosystem
                    +trait_directionality
                    +sqrt_ni_1,
                    random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                    R = list(phylo=cor), test="t",
                    method="ML", data=dat_ES_final_2)
#no trait 
full_modmm <- rma.mv(yi=yi, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +duration_standard
                     +exp_age*experiment_type
                     +size*experiment_type
                     +ecosystem
                     +sqrt_ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)

#no age, trait, size 

full_modnn <- rma.mv(yi=yi, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +duration_standard
                     +ecosystem
                     +sqrt_ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)


model.sel(full_mode, full_modee, full_modff, full_modgg, full_modhh, full_modii, full_modjj, 
          full_modkk, full_modll, full_modmm, full_modnn)

model.sel(full_mod, full_moda, full_modb, full_modc, full_modd, 
          full_mode, full_modf, full_modg, full_modh, full_modi, full_modj, full_modk, full_modl, 
          full_modaa, full_modbb, full_modcc, full_moddd, 
          full_modee, full_modff, full_modgg, full_modhh, full_modii, full_modjj, 
          full_modkk, full_modll, full_modmm, full_modnn)


#full_mode prevails
model.sel(full_mode, full_modf, full_modg, full_modh, full_modi, full_modj, full_modk, full_modl)
model.sel(full_mod, full_moda, full_modb, full_modc, full_modd)
#full_mode prevails 

#model averaging the top models
eval(metafor:::.MuMIn)
final_model <- summary(model.avg(model.sel(full_mode, full_modh)))
#generating CI's for the top model
confint(final_model)


#################################################
###################Figure 3######################
#################################################

# New facet label names for age variable
age.labs <- c("larval", "juvenile", "adult")
names(age.labs) <- c("0", "1", "2")
# New facet label names for size variable
size.labs <- c("small", "medium", "large")
names(size.labs) <- c("1", "2", "3")

#subsetting the data for the histograms
plot_dat <- dat_ES_final_2 %>%
  filter(between(yi, -2.5, 2.5))


vline_age <- dat_ES_final_2 %>%
  group_by(experiment_type, exp_age) %>%
  summarise(mean=mean(yi), 
            median=median(yi), 
            n=n())%>%
  cbind(coef=c(0.082753, 0.082753+0.1375499, 0.082753+0.2242981, 0.082753+0.0984037, 0.082753+0.0984037+0.1375499+-0.3834808 , 0.082753+0.0984037+0.2242981+-0.2975666))

hist_age <- ggplot()+
  geom_histogram(data=plot_dat, aes(x=yi, colour=experiment_type), fill="grey", bins = 100, alpha=0.5)+
  scale_fill_manual(values=c("#440154", "yellow3"))+
  scale_colour_manual(values=c("#440154", "yellow3"))+
  facet_wrap(~exp_age,  labeller = labeller(exp_age = age.labs))+
  geom_vline(data = vline_age,
             aes(xintercept = coef, colour=experiment_type), size=1)+
  theme_linedraw()+
  geom_vline(data = vline_age,
             aes(xintercept = median, colour=experiment_type), linetype="longdash")+
  geom_label(data=vline_age, aes(x=0.5, label=n, y=50, fill=experiment_type), 
             position = position_dodge(width=4), hjust=1, alpha=0.5, size=6)+
  xlab("Effect size")+
  ylab("count")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))+
  ggtitle("Age class")


vline_size <- dat_ES_final_2 %>%
  group_by(experiment_type, size) %>%
  summarise(mean=mean(yi), 
            median=median(yi), 
            n=n()) %>%
  cbind(coef=c(0.082753, 0.082753+0.2307199, 0.082753+0.7607518, 
               0.082753+0.0984037, 0.082753+0.0984037+0.2307199+0.2132700, 
               0.082753+0.0984037+0.7607518+3.8866970))

hist_size <- ggplot()+
  geom_histogram(data=plot_dat, aes(x=yi, colour=experiment_type), fill="grey", bins = 100, alpha=0.5)+
  xlim(-2.5, 2.5)+
  scale_fill_manual(values=c("#440154", "yellow3"))+
  scale_colour_manual(values=c("#440154", "yellow3"))+
  facet_wrap(~size,  labeller = labeller(size = size.labs))+
  geom_vline(data = vline_size,
             aes(xintercept = median, colour=experiment_type), linetype="longdash")+
  geom_vline(data = vline_size,
             aes(xintercept = coef, colour=experiment_type), size=1)+
  theme_linedraw()+
  geom_label(data=vline_size, aes(x=0.5, label=n, y=60, fill=experiment_type), 
             position = position_dodge(width=4), hjust=1, alpha=0.5, size=6)+
  xlab("")+
  ylab("count")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))+
  ggtitle("Size class")


vline_trait <- dat_ES_final_2 %>%
  group_by(experiment_type, trait_directionality) %>%
  summarise(mean=mean(yi), 
            median=median(yi), 
            n=n())%>%
  cbind(coef=c(0.082753, 0.082753+-0.0223745, 0.082753+0.0984037, 0.082753+-0.0223745+ 0.0984037+0.3651103))

hist_trait <- ggplot()+
  geom_histogram(data=plot_dat, aes(x=yi, colour=experiment_type), fill="grey", bins = 100, alpha=0.5)+
  scale_fill_manual(values=c("#440154", "yellow3"))+
  scale_colour_manual(values=c("#440154", "yellow3"))+
  facet_wrap(~trait_directionality)+
  geom_vline(data = vline_trait,
             aes(xintercept = median, colour=experiment_type), linetype="longdash")+
  geom_vline(data = vline_trait,
             aes(xintercept = coef, colour=experiment_type), size=1)+
  theme_linedraw()+
  theme_linedraw()+
  geom_label(data=vline_trait, aes(x=0.5, label=n, y=60, fill=experiment_type), 
             position = position_dodge(width=4), hjust=1, alpha=0.5, size=6)+
  xlab("")+
  ylab("count")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))+
  ggtitle("Trait directionality")


vline_exper <- dat_ES_final_2 %>%
  group_by(experiment_type) %>%
  summarise(mean=mean(yi), 
            median=median(yi), 
            n=n())%>%
  cbind(coef=c(0.082753, 0.082753+ 0.0984037))

hist_exper <- ggplot()+
  geom_histogram(data=plot_dat, aes(x=yi, colour=experiment_type), fill="grey", bins = 100, alpha=0.5)+
  scale_fill_manual(values=c("#440154", "yellow3"))+
  scale_colour_manual(values=c("#440154", "yellow3"))+
  facet_wrap(~experiment_type)+
  geom_vline(data = vline_exper,
             aes(xintercept = median, colour=experiment_type), linetype="longdash")+
  geom_vline(data = vline_exper,
             aes(xintercept = coef, colour=experiment_type), size=1)+
  theme_linedraw()+
  geom_label(data=vline_exper, aes(x=1, label=n, y=75, fill=experiment_type), 
             position = position_dodge(width=1), hjust=1, alpha=0.5, size=6)+
  xlab("")+
  ylab("count")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))+
  ggtitle("Experiment type")

hist_exper/ hist_trait / hist_size / hist_age 
ggsave("figures/fig3.png", width=15, height=13, dpi=300)




#################################################
###################Figure 4######################
#################################################

mean_temp<- dat_ES_final_2 %>%
  filter(between(yi, -5, 5)) %>%
  ggplot(aes(x=mean_temp_constant, y=yi))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.3, size=2, colour="gray")+
  scale_colour_viridis(option="D")+
  #geom_smooth(method="lm")+
  theme_bw()+
  xlab("Mean temperature (°C)")+
  ylab("Effect size")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))+
  geom_abline(intercept = 0.0829, slope = 0.0104, color="#440154", size=1)+
  geom_abline(intercept = 0.0829+0.0827, slope = -0.0064+0.0104, color="yellow3", size=1)
  #geom_label(x=15, y=4, aes(label="n=1653"), size=8, colour="black")


flux_range <- dat_ES_final_2 %>%
  filter(between(yi, -5, 5)) %>%
  ggplot(aes(x=flux_range, y=yi))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.3, size=2, colour="gray")+
  scale_colour_viridis(option="D")+
  #geom_smooth(method="lm")+
  theme_bw()+
  xlab("Fluctuation magnitude (°C)")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))+
  geom_abline(intercept = 0.0829, slope = -0.0108, color="#440154", size=1)+
  geom_abline(intercept = 0.0829+0.0827, slope = -0.0108+-0.0085, color="yellow3", size=1)
  #geom_label(x=5, y=4, aes(label="n=1653"), size=8, colour="black")

duration<- dat_ES_final_2 %>%
  filter(between(yi, -5, 5)) %>%
  ggplot(aes(x=duration_standard, y=yi))+
  xlim(0,250)+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(colour="gray",alpha=0.3, size=2)+
  scale_colour_viridis(option="D", discrete=TRUE)+
  theme_bw()+
  xlab("Duration (days)")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))+
  geom_abline(intercept = 0.0829, slope = -0.0222 , color="#440154", size=1)+
  geom_abline(intercept = 0.0829+0.0827, slope = -0.0222 , color="yellow3", size=1)
  #geom_label(x=50, y=4, aes(label="n=1653"), size=8, colour="black")

secondary_temp <- dat_ES_final_2 %>%
  filter(between(yi, -5, 5)) %>%
  #filter(between(duration, 0, 200)) %>%
  ggplot(aes(x=secondary_temp, y=yi))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.3, size=2, colour="gray")+
  scale_colour_viridis(option="D")+
  theme_bw()+
  xlab("Secondary temperature (°C)")+
  ylab("Effect size")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))+
  geom_abline(intercept = 0.082753, slope = 0.028821 , color="#440154", size=1)+
  geom_abline(intercept = 0.082753+0.0984037, slope = 0.028821 , color="yellow3", size=1)
  #geom_label(x=5, y=4, aes(label="n=1653"), size=8, colour="black")

mean_temp + flux_range +  secondary_temp + duration 
ggsave("figures/fig4.png", width=10, height=8.5, dpi=350)





###############################################
##############SI Tables and Plots##############
###############################################
#funnel plot
normalized <- dat_ES_final_2 %>%
  filter(between(yi, -10, 10))
dummy_full_mod <- rma.mv(yi, V,
                         random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
                         R = list(phylo=cor),
                         method="ML", data=dat_ES_final_2)
funnel(dummy_full_mod)

##histogram
dat_ES_final_2 %>%
  #filter(between(yi, -10, 10)) %>%
  ggplot(aes(x=yi))+
  geom_histogram(bins=100)+
  theme_bw()+
  xlab("SMD")

##I2 calculations from best models 
#phylo level heterogeneity
full_mode$sigma2[1]/sum(full_mode$sigma2[1], full_mode$sigma2[2], full_mode$sigma2[3]) #0.9808726
full_modh$sigma2[1]/sum(full_modh$sigma2[1], full_modh$sigma2[2], full_modh$sigma2[3]) #0.9835042
#study level heterogeneity
full_mode$sigma2[2]/sum(full_mode$sigma2[1], full_mode$sigma2[2], full_mode$sigma2[3]) #0.01353682
full_modh$sigma2[2]/sum(full_modh$sigma2[1], full_modh$sigma2[2], full_modh$sigma2[3]) #0.01131141
#response heterogeneity 
full_mode$sigma2[3]/sum(full_mode$sigma2[1], full_mode$sigma2[2], full_mode$sigma2[3]) #0.005590618
full_modh$sigma2[3]/sum(full_modh$sigma2[1], full_modh$sigma2[2], full_modh$sigma2[3]) #0.005184336

#calculating fail safe number
fsn(yi, vi, data=dat_ES_final_2,type="Rosenberg") #261304

#correlation between age and size
dat_cor <- dat_ES_final_2 %>%
  mutate(exp_age=as.numeric(exp_age), 
         size=as.numeric(size))

res0 <- cor.test(dat_cor$exp_age, dat_cor$size, 
                 method = "pearson")

#model coefficient output
model_df <- as.data.frame(coef(final_model)) %>% tibble::rownames_to_column("covariates") %>%
  mutate_if(is.numeric, ~round(., 3))
write_csv(model_df, "model_df.csv")

#frequency table
Var.tbl <- xtabs(~experiment_type+size+exp_age+trait_directionality+ecosystem, dat_ES_final_2)
Var.dbf <- as.data.frame.table(Var.tbl)

write_csv(Var.dbf, "var_dbf.csv")




