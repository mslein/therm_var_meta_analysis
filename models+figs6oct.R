#####################################################
## Manuscript coding, models, and figures ###########
#####################################################


##loading libraries 
pacman::p_load(tidyverse, metafor, patchwork, rotl, MuMIn, viridis, ape, wesanderson, glmulti, MuMIn, metaAidR, corrplot)
#loading in the raw data file
dat_0<- read_csv("data/rawdata23dec22.csv") %>%
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
                               "coturnix japonica"))  %>%
  mutate(ecosystem = case_when(environment == "terrestrial" ~ "terrestrial", 
                               TRUE ~ "aquatic/marine")) %>%
  mutate(size =case_when(size %in% c("0","1") ~ "1", 
                         size %in% c("2") ~ "2", 
                         size %in% c("3", "4") ~ "3")) %>%
  mutate(experiment_type = as.factor(experiment_type), 
         exp_age=as.factor(exp_age), 
         size=as.factor(size), 
         trait_directional=as.factor(trait_directional), 
         ecosystem=as.factor(ecosystem))

unique(dat_ES_final$original_name)#85

#specifying the reference level for each dummy variable, with the level with the most observations
dat_ES_final$experiment_type <- relevel(dat_ES_final$experiment_type, "acute")
dat_ES_final$exp_age <- relevel(dat_ES_final$exp_age, "1")
dat_ES_final$size <- relevel(dat_ES_final$size, "1")
dat_ES_final$trait_directional <- relevel(dat_ES_final$trait_directional, "positive")
dat_ES_final$ecosystem <- relevel(dat_ES_final$ecosystem, "terrestrial")

#getting the phylogenies
taxa <- tnrs_match_names(as.character(unique(dat_ES_final$original_name)))
taxa_df <- taxa %>% 
  as.data.frame() %>% 
  select(search_string, unique_name) %>% 
  rename(original_name=search_string)
#creating the tree
tr <- tol_induced_subtree(ott_id(taxa)[is_in_tree(ott_id(taxa))]) #82 species
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
  mutate(ni_1 =(flux_samp+constant_samp)/(flux_samp * constant_samp)) %>%
  mutate(scale_effect = case_when(trait_directional %in% c("negative") ~ yi*-1, 
                                  TRUE ~ yi)) %>%
  filter(trait_directional != "neutral")


dat_ES_final_2 %>%
  ggplot(aes(x=scale_effect, fill=trait_directional))+
  geom_histogram()+
  facet_wrap(~trait_directional)

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
int_mod <- rma.mv(yi=scale_effect, V=V, 
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)




####candidate models
#full model
full_mod<-rma.mv(yi=scale_effect, V=V, 
                 mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                 +experiment_type
                 +I(secondary_temp - mean(secondary_temp))
                 +log10(duration_standard)
                 +exp_age
                 +size
                 +ecosystem
                 +trait_directional
                 +ni_1,
                 random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                 R = list(phylo=cor), test="t",
                 method="ML", data=dat_ES_final_2)
#without ecosystem
full_moda<-rma.mv(yi=scale_effect, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +experiment_type
                  +I(secondary_temp - mean(secondary_temp))
                  +log10(duration_standard)
                  +exp_age
                  +size
                  +trait_directional
                  +ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

#without duration
full_modb<-rma.mv(yi=scale_effect, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +experiment_type
                  +I(secondary_temp - mean(secondary_temp))
                  +exp_age
                  +size
                  +ecosystem
                  +trait_directional
                  +ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)
#without secondary temp
full_modc<-rma.mv(yi=scale_effect, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +experiment_type
                  +log10(duration_standard)
                  +exp_age
                  +size
                  +ecosystem
                  +trait_directional
                  +ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)
#without publication bias
full_modd<-rma.mv(yi=scale_effect, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +experiment_type
                  +I(secondary_temp - mean(secondary_temp))
                  +log10(duration_standard)
                  +exp_age
                  +size
                  +ecosystem
                  +trait_directional,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

model.sel(full_mod, full_moda, full_modb, full_modc, full_modd)
#full_moda prevails

# hypothesis 1
#full model with interaction between experiment_type and the mods
full_mode<-rma.mv(yi=scale_effect, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))*experiment_type
                  +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                  +I(secondary_temp - mean(secondary_temp))
                  +log10(duration_standard)
                  +exp_age*experiment_type
                  +size*experiment_type
                  +trait_directional*experiment_type
                  +ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor),test="t",
                  method="ML", data=dat_ES_final_2)

#no experiment type 
full_modf<-rma.mv(yi=scale_effect, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +log10(duration_standard)
                  +I(secondary_temp - mean(secondary_temp))
                  +exp_age
                  +size
                  +trait_directional
                  +ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

#with experiment type but no interaction 
full_modg<-rma.mv(yi=scale_effect, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +experiment_type
                  +I(secondary_temp - mean(secondary_temp))
                  +log10(duration_standard)
                  +exp_age
                  +size
                  +trait_directional
                  +ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

#no interaction between experiment type and flux range
full_modh<-rma.mv(yi=scale_effect, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))
                  +I(secondary_temp - mean(secondary_temp))
                  +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                  +log10(duration_standard)
                  +exp_age*experiment_type
                  +size*experiment_type
                  +trait_directional*experiment_type
                  +ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

#no interaction between experiment type and mean temp
full_modi<-rma.mv(yi=scale_effect, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))*experiment_type
                  +I(mean_temp_constant-mean(mean_temp_constant))
                  +log10(duration_standard)
                  +I(secondary_temp - mean(secondary_temp))
                  +exp_age*experiment_type
                  +size*experiment_type
                  +trait_directional*experiment_type
                  +ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

#no interaction between experiment type and age
full_modj<-rma.mv(yi=scale_effect, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))*experiment_type
                  +I(secondary_temp - mean(secondary_temp))
                  +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                  +log10(duration_standard)
                  +exp_age
                  +size*experiment_type
                  +trait_directional*experiment_type
                  +ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)
#no interaction between experiment type and size
full_modk<-rma.mv(yi=scale_effect, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))*experiment_type
                  +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                  +log10(duration_standard)
                  +I(secondary_temp - mean(secondary_temp))
                  +exp_age*experiment_type
                  +size
                  +trait_directional*experiment_type
                  +ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

#no interaction between experiment type and trait directionality
full_modl<-rma.mv(yi=scale_effect, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))*experiment_type
                  +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                  +log10(duration_standard)
                  +I(secondary_temp - mean(secondary_temp))
                  +exp_age*experiment_type
                  +size*experiment_type
                  +trait_directional
                  +ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

model.sel(full_mode, full_modf, full_modg, full_modh, full_modi, full_modj, full_modk, full_modl)
model.sel(full_mod, full_moda, full_modb, full_modc, full_modd)
#full mod e + full mod j are the same, proceed with full model e

#h2 and h3 with interaction between flux range and experiment (based on full_mode as base)
# hypothesis 2 
# no mean temp
full_modaa<-rma.mv(yi=scale_effect, V=V, 
                   mods = ~ I(flux_range-mean(flux_range))*experiment_type
                   +log10(duration_standard)
                   +I(secondary_temp - mean(secondary_temp))
                   +exp_age*experiment_type
                   +size*experiment_type
                   +trait_directional*experiment_type
                   +ni_1,
                   random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                   R = list(phylo=cor), test="t",
                   method="ML", data=dat_ES_final_2)
#no flux range
full_modbb <- rma.mv(yi=scale_effect, V=V, 
                     mods = ~ I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +log10(duration_standard)
                     +I(secondary_temp - mean(secondary_temp))
                     +exp_age*experiment_type
                     +size*experiment_type
                     +trait_directional*experiment_type
                     +ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no interaction between flux range and mean temp
full_modcc <- rma.mv(yi=scale_effect, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*experiment_type
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +log10(duration_standard)
                     +exp_age*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +size*experiment_type
                     +trait_directional*experiment_type
                     +ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no flux range or mean temp
full_moddd <- rma.mv(yi=scale_effect, V=V, 
                     mods = ~ log10(duration_standard)
                     +exp_age*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +size*experiment_type
                     +trait_directional*experiment_type
                     +ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)

model.sel(full_mode, full_modj,
          full_modaa, full_modbb, full_modcc, full_moddd)
#full mod e + full mod j are the same, proceed with full model e
model.sel(full_mode, full_modf, full_modg, full_modh, full_modi, full_modj, full_modk, full_modl)
model.sel(full_mod, full_moda, full_modb, full_modc, full_modd)

#hypothesis 3
#no interaction with age (based on full_mode)
full_modee <- rma.mv(yi=scale_effect, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +log10(duration_standard)
                     +exp_age
                     +size*experiment_type
                     +trait_directional*experiment_type
                     +ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no age
full_modff <- rma.mv(yi=scale_effect, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +log10(duration_standard)
                     +size*experiment_type
                     +trait_directional*experiment_type
                     +ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no interaction with size
full_modgg <- rma.mv(yi=scale_effect, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +log10(duration_standard)
                     +exp_age*experiment_type
                     +size
                     +trait_directional*experiment_type
                     +ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no size
full_modhh <- rma.mv(yi=scale_effect, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +log10(duration_standard)
                     +exp_age*experiment_type
                     +trait_directional*experiment_type
                     +ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)

#no interaction with trait 
full_modii <- rma.mv(yi=scale_effect, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +log10(duration_standard)
                     +exp_age*experiment_type
                     +size*experiment_type
                     +trait_directional
                     +ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no age or size
full_modjj <- rma.mv(yi=scale_effect, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +log10(duration_standard)
                     +trait_directional*experiment_type
                     +ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no interaction with age or size 
full_modkk <- rma.mv(yi=scale_effect, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +log10(duration_standard)
                     +exp_age
                     +size
                     +trait_directional*experiment_type
                     +ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)
#no interaction with trait 
full_modll<- rma.mv(yi=scale_effect, V=V, 
                    mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                    +I(flux_range-mean(flux_range))*experiment_type
                    +I(secondary_temp - mean(secondary_temp))
                    +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                    +log10(duration_standard)
                    +exp_age*experiment_type
                    +size*experiment_type
                    +trait_directional
                    +ni_1,
                    random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                    R = list(phylo=cor), test="t",
                    method="ML", data=dat_ES_final_2)
#no trait 
full_modmm <- rma.mv(yi=scale_effect, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +log10(duration_standard)
                     +exp_age*experiment_type
                     +size*experiment_type
                     +ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)

#no age, trait, size 
full_modnn <- rma.mv(yi=scale_effect, V=V, 
                     mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                     +I(flux_range-mean(flux_range))*experiment_type
                     +I(secondary_temp - mean(secondary_temp))
                     +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                     +log10(duration_standard)
                     +ni_1,
                     random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                     R = list(phylo=cor), test="t",
                     method="ML", data=dat_ES_final_2)

#full mod e + full mod j are the same, proceed with full model e
model.sel(full_mode, full_modj, full_modee, full_modff, full_modgg, full_modhh, full_modii, full_modjj, 
          full_modkk, full_modll, full_modmm, full_modnn)

model.sel(full_mod, full_moda, full_modb, full_modc, full_modd, 
          full_mode, full_modf, full_modg, full_modh, full_modi, full_modj, full_modk, full_modl, 
          full_modaa, full_modbb, full_modcc, full_moddd, 
          full_modee, full_modff, full_modgg, full_modhh, full_modii, full_modjj, 
          full_modkk, full_modll, full_modmm, full_modnn)


#model averaging the top models
eval(metafor:::.MuMIn)
final_model <- summary(model.avg(model.sel(full_mode, full_modj)))
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
  filter(between(scale_effect, -2.5, 2.5))

pp1_labels <- data.frame(exp_age = c(0,1,2), labels=c("C", "D", "E"))

vline_age <- dat_ES_final_2 %>%
  group_by(experiment_type, exp_age) %>%
  summarise(mean=mean(scale_effect), 
            median=median(scale_effect), 
            n=n())%>%
  cbind(coef=c(-0.186, 
               -0.186+-0.154, 
               -0.186+-0.004, 
               -0.186+-0.073 , 
               -0.186+-0.073+-0.154+0.107, 
               -0.186+-0.073+-0.004+0.131))

pp1 <- ggplot()+
  geom_point(data=plot_dat, aes(x=experiment_type, y=scale_effect, colour=experiment_type), alpha=0.1)+
  #geom_violin(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, fill=experiment_type), alpha=0.3)+
  stat_summary(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, colour=experiment_type), 
               fun.data="mean_cl_boot", 
               geom = "crossbar", 
               width=0.75, 
               alpha=0.5)+
  scale_fill_manual(values=c("#440154", "yellow3"))+
  scale_colour_manual(values=c("#440154", "yellow3"))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(data = vline_age,
             aes(y = coef, x=experiment_type, fill=experiment_type),colour="black",pch=21, size=5)+
  theme_linedraw()+
  ylab("Effect size")+
  xlab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        #axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_label(data=vline_age, aes(x=experiment_type, label=n, y=-2.6, fill=experiment_type),
alpha=0.5, size=6)+
  facet_wrap(~exp_age,  labeller = labeller(exp_age = age.labs))+
  ggtitle("Age class")+
  geom_text(data=pp1_labels, aes(x=0.6, y=2.5,  label=labels), colour="black", size=10)


pp2_labels <- data.frame(experiment_type=c( "acute", "acclimation"), labels=c("B", "A"))


vline_exper <- dat_ES_final_2 %>%
  group_by(experiment_type) %>%
  summarise(mean=mean(scale_effect), 
            median=median(scale_effect), 
            n=n())%>%
  cbind(coef=c(-0.186, -0.186+ -0.073))

pp2 <- ggplot()+
  geom_point(data=plot_dat, aes(x=experiment_type, y=scale_effect, colour=experiment_type), alpha=0.1)+
  #geom_violin(data=plot_dat, aes(x=experiment_type, y=scale_effect, fill=experiment_type), alpha=0.3)+
  stat_summary(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, colour=experiment_type), 
               fun.data="mean_cl_boot", 
               geom = "crossbar", 
               width=0.75, 
               alpha=0.5)+
  scale_fill_manual(values=c("#440154", "yellow3"))+
  scale_colour_manual(values=c("#440154", "yellow3"))+
  geom_hline(yintercept = 0, linetype='dashed')+
  #facet_wrap(~experiment_type)+
  geom_point(data = vline_exper,
             aes(y = coef, x=experiment_type, fill=experiment_type),colour="black",pch=21, size=5)+
  theme_linedraw()+
  ylab("Effect size")+
  xlab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        #axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_label(data=vline_exper, aes(x=experiment_type, label=n, y=-2.7, fill=experiment_type),
             alpha=0.5, size=6)+
  ggtitle("Experiment type")+
  facet_grid(~experiment_type, scales = "free_x", space = "free_x")+
  geom_text(data=pp2_labels, aes(x=0.5, y=2.5,  label=labels), colour="black", size=10)

pp3_labels <- data.frame(size=c( "1", "2", "3"), labels=c("F", "G", "H"))

vline_size <- dat_ES_final_2 %>%
  group_by(experiment_type, size) %>%
  summarise(mean=mean(scale_effect), 
            median=median(scale_effect), 
            n=n()) %>%
  cbind(coef=c(-0.186 , 
               -0.186 +-0.311 , 
               -0.186 +-0.062, 
               -0.186 +-0.073 , 
               -0.186 +-0.073 +-0.311 +0.255 , 
               -0.186 +-0.073 +-0.062 +-0.226))

pp3 <- ggplot()+
  geom_point(data=plot_dat, aes(x=experiment_type, y=scale_effect, colour=experiment_type), alpha=0.1)+
  #geom_violin(data=plot_dat, aes(x=experiment_type, y=scale_effect, fill=experiment_type), alpha=0.3)+
  stat_summary(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, colour=experiment_type), 
               fun.data="mean_cl_boot", 
               geom = "crossbar", 
               width=0.75, 
               alpha=0.5)+
  scale_fill_manual(values=c("#440154", "yellow3"))+
  scale_colour_manual(values=c("#440154", "yellow3"))+
  geom_hline(yintercept = 0, linetype='dashed')+
  facet_wrap(~size,  labeller = labeller(size = size.labs))+
  geom_point(data = vline_size,
             aes(y = coef, x=experiment_type, fill=experiment_type),colour="black",pch=21, size=5)+
  theme_linedraw()+
  ylab("Effect size")+
  xlab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        #axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_label(data=vline_size, aes(x=experiment_type, label=n, y=-2.6, fill=experiment_type),
             alpha=0.5, size=6)+
  ggtitle("Size class")+
  geom_text(data=pp3_labels, aes(x=0.6, y=2.5,  label=labels), colour="black", size=10)


pp4_labels <- data.frame(trait_directional=c("positive", "negative"), labels=c("I", "J"))

vline_trait <- dat_ES_final_2 %>%
  group_by(experiment_type, trait_directional) %>%
  summarise(mean=mean(scale_effect), 
            median=median(scale_effect), 
            n=n())%>%
  cbind(coef=c(-0.186, 
               -0.186+-0.034, 
               -0.186+-0.073, 
               -0.186+-0.034+-0.073+-0.188))

pp4 <- ggplot()+
  geom_point(data=plot_dat, aes(x=experiment_type, y=scale_effect, colour=experiment_type), alpha=0.1)+
  #geom_violin(data=plot_dat, aes(x=experiment_type, y=scale_effect, fill=experiment_type), alpha=0.3)+
  stat_summary(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, colour=experiment_type), 
               fun.data="mean_cl_boot", 
               geom = "crossbar", 
               width=0.75, 
               alpha=0.5)+
  scale_fill_manual(values=c("#440154", "yellow3"))+
  scale_colour_manual(values=c("#440154", "yellow3"))+
  geom_hline(yintercept = 0, linetype='dashed')+
  facet_wrap(~trait_directional)+
  geom_point(data = vline_trait,
             aes(y = coef, x=experiment_type, fill=experiment_type),colour="black",pch=21, size=5)+
  theme_linedraw()+
  ylab("Effect size")+
  xlab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        #axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_label(data=vline_trait, aes(x=experiment_type, label=n, y=-2.6, fill=experiment_type),
             alpha=0.5, size=6)+
ggtitle("Trait directionality")+
geom_text(data=pp4_labels, aes(x=0.5, y=2.5,  label=labels), colour="black", size=10)

pp2 + pp1 + pp3+ pp4
ggsave("figures/fig3_mean.png", width=15, height=10, dpi=350)

#################################################
###################Figure 4######################
#################################################

mean_labels <- data.frame(experiment_type = c("acclimation", "acute"), labels=c("A", "B"))

mean_temp<- dat_ES_final_2 %>%
  filter(between(scale_effect, -2.5, 2.5)) %>%
  ggplot(aes(x=I(mean_temp_constant-mean(mean_temp_constant)), y=scale_effect, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.3, size=2)+
  scale_colour_manual(values =c("#440154", "yellow3"))+
  #geom_smooth(method="lm")+
  theme_bw()+
  xlab("Standardized mean temperature (°C)")+
  ylab("Effect size")+
  theme_linedraw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))+
  geom_abline(aes(intercept=i, slope=s, colour=experiment_type), 
              data=data.frame(experiment_type=c("acute","acclimation"), i=c(-0.186, -0.186 + -0.073), s=c(-0.061,-0.061+ 0.0572782))) +
  facet_wrap(~experiment_type)+
  geom_text(data=mean_labels, aes(x=-17, y=2.4,  label=labels), colour="black", size=10)
  


flux_labels <- data.frame(experiment_type = c("acclimation", "acute"), labels=c("C", "D"))


flux_range <- dat_ES_final_2 %>%
  filter(between(scale_effect, -2.5, 2.5)) %>%
  ggplot(aes(x=I(flux_range-mean(flux_range)), y=scale_effect, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.3, size=2)+
  scale_colour_manual(values =c("#440154", "yellow3"))+
  theme_linedraw()+
  #geom_smooth(method="lm")+
  xlab("Standardized fluctuation magnitude (°C)")+
  ylab("Effect size")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))+
  geom_abline(aes(intercept=i, slope=s, colour=experiment_type), 
              data=data.frame(experiment_type=c("acute","acclimation"), i=c(-0.186, -0.186 + -0.073), s=c(-0.010,-0.0108231+-0.0046560))) +
  facet_wrap(~experiment_type)+
  geom_text(data=flux_labels, aes(x=-8, y=2.4,  label=labels), colour="black", size=10)



duration_labels <- data.frame(experiment_type = c("acclimation", "acute"), labels=c("G", "H"))

duration<- dat_ES_final_2 %>%
  filter(between(scale_effect, -2.5, 2.5)) %>%
  ggplot(aes(x=log10(duration_standard), y=scale_effect, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.3, size=2)+
  scale_colour_manual(values =c("#440154", "yellow3"))+
  theme_linedraw()+
  #geom_smooth(method="lm")+
  xlab("Duration (log(days))")+
  ylab("Effect size")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))+
  geom_abline(aes(intercept=i, slope=s, colour=experiment_type), 
              data=data.frame(experiment_type=c("acute","acclimation"), i=c(-0.186, -0.186 + -0.073), s=c(-0.052,-0.052))) +
  geom_text(data=duration_labels, aes(x=-0.275, y=2.4,  label=labels), colour="black", size=10)+
  facet_wrap(~experiment_type)
  


secondary_labels <- data.frame(experiment_type = c("acclimation", "acute"), labels=c("E", "F"))

secondary_temp <- dat_ES_final_2 %>%
  filter(between(scale_effect, -2.5, 2.5)) %>%
  ggplot(aes(x=I(secondary_temp-mean(secondary_temp)), y=scale_effect, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.3, size=2)+
  scale_colour_manual(values =c("#440154", "yellow3"))+
  theme_linedraw()+
  #geom_smooth(method="lm")+
  xlab("Standardized secondary temperature (°C)")+
  ylab("Effect size")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))+
  geom_abline(aes(intercept=i, slope=s, colour=experiment_type), 
              data=data.frame(experiment_type=c("acute","acclimation"), i=c(-0.186, -0.186 + -0.073), s=c(0.0002018,0.0002018))) +
  geom_text(data=secondary_labels, aes(x=-33, y=2.4,  label=labels), colour="black", size=10)+
  #data=data.frame(experiment_type =c("acute", "acclimation"), y=3.35, x=10, label=c("*")), size=20)+
  facet_wrap(~experiment_type)
  

mean_temp + flux_range +  secondary_temp + duration 
ggsave("figures/fig4.png", width=14, height=8.5, dpi=350)




###############################################
##############SI Tables and Plots##############
###############################################
#funnel plot
dummy_full_mod <- rma.mv(scale_effect, V,
                         random = list(~1 | phylo, ~1 | study_id, ~1 | response_id),
                         R = list(phylo=cor),
                         method="ML", data=normalized)
funnel(dummy_full_mod)

##histogram
dat_ES_final_2 %>%
  #filter(between(yi, -10, 10)) %>%
  ggplot(aes(x=scale_effect))+
  geom_histogram(bins=100)+
  theme_bw()+
  xlab("SMD")

#flux range X mean temp interaction plot
dat_ES_final_2 %>%
  filter(between(scale_effect, -5, 5)) %>%
  filter(flux_range %in% c(4,6,8,10,12,15))%>%
  mutate(flux_range2 = round(I(flux_range-mean(flux_range)), 1), 
         mean_temp2= I(mean_temp_constant-mean(mean_temp_constant))) %>%
  ggplot(aes(y=scale_effect, x=mean_temp2, group=flux_range2, colour=flux_range))+
  geom_jitter()+
  scale_colour_viridis(name="Fluctuation magnitude (°C)")+
  geom_smooth(method="lm")+
  facet_grid(~flux_range2)+
  theme_bw()+
  xlab("Standardized mean temperature (°C)")+
  ylab("Effect size")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        #legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))

###figure s4
pp1S_labels <- data.frame(exp_age = c(0,1,2), 
                         labels=c("C", "D", "E"),
                         y=c(12.5, 33.5, 17.5),
                         x=c(0.77, 0.75, 0.75))

vline_age <- dat_ES_final_2 %>%
  group_by(experiment_type, exp_age) %>%
  summarise(mean=mean(scale_effect), 
            median=median(scale_effect), 
            n=n())%>%
  cbind(coef=c(-0.186, 
               -0.186+-0.154, 
               -0.186+-0.004, 
               -0.186+-0.073 , 
               -0.186+-0.073+-0.154+0.107, 
               -0.186+-0.073+-0.004+0.131))

pp1S <- ggplot()+
  geom_point(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, colour=experiment_type), alpha=0.1)+
  #geom_violin(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, fill=experiment_type), alpha=0.3)+
  stat_summary(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, colour=experiment_type), 
               fun.data="mean_cl_boot", 
               geom = "crossbar", 
               width=0.75, 
               alpha=0.5)+
  scale_fill_manual(values=c("#440154", "yellow3"))+
  scale_colour_manual(values=c("#440154", "yellow3"))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(data = vline_age,
             aes(y = coef, x=experiment_type, fill=experiment_type),colour="black",pch=21, size=5)+
  theme_linedraw()+
  ylab("Effect size")+
  xlab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        #axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_label(data=vline_age, aes(x=experiment_type, label=n, y=-10, fill=experiment_type),
             alpha=0.5, size=6)+
  facet_wrap(~exp_age,  labeller = labeller(exp_age = age.labs), scales="free")+
  ggtitle("Age class")+
  geom_text(data=pp1S_labels, aes(x=x, y=y,  label=labels), colour="black", size=10)


pp2S_labels <- data.frame(experiment_type=c( "acute", "acclimation"), 
                         labels=c("B", "A"),
                         y=c(33.5, 33.5),
                         x=c(0.5, 0.5))


vline_exper <- dat_ES_final_2 %>%
  group_by(experiment_type) %>%
  summarise(mean=mean(scale_effect), 
            median=median(scale_effect), 
            n=n())%>%
  cbind(coef=c(-0.186, -0.186+ -0.073))

pp2S <- ggplot()+
  geom_point(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, colour=experiment_type), alpha=0.1)+
  #geom_violin(data=plot_dat, aes(x=experiment_type, y=scale_effect, fill=experiment_type), alpha=0.3)+
  stat_summary(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, colour=experiment_type), 
               fun.data="mean_cl_boot", 
               geom = "crossbar", 
               width=0.75, 
               alpha=0.5)+
  scale_fill_manual(values=c("#440154", "yellow3"))+
  scale_colour_manual(values=c("#440154", "yellow3"))+
  geom_hline(yintercept = 0, linetype='dashed')+
  #facet_wrap(~experiment_type)+
  geom_point(data = vline_exper,
             aes(y = coef, x=experiment_type, fill=experiment_type),colour="black",pch=21, size=5)+
  theme_linedraw()+
  ylab("Effect size")+
  xlab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        #axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_label(data=vline_exper, aes(x=experiment_type, label=n, y=-5, fill=experiment_type),
             alpha=0.5, size=6)+
  ggtitle("Experiment type")+
  facet_grid(~experiment_type, scales = "free", space = "free")+
  geom_text(data=pp2S_labels, aes(x=x, y=y,  label=labels), colour="black", size=10)

pp3S_labels <- data.frame(size=c( "1", "2", "3"),
                         labels=c("F", "G", "H"),
                         y=c(17, 35, 2.7),
                         x=c(0.75, 0.75, 0.75))

vline_size <- dat_ES_final_2 %>%
  group_by(experiment_type, size) %>%
  summarise(mean=mean(scale_effect), 
            median=median(scale_effect), 
            n=n()) %>%
  cbind(coef=c(-0.186 , 
               -0.186 +-0.311 , 
               -0.186 +-0.062, 
               -0.186 +-0.073 , 
               -0.186 +-0.073 +-0.311 +0.255 , 
               -0.186 +-0.073 +-0.062 +-0.226))

pp3S <- ggplot()+
  geom_point(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, colour=experiment_type), alpha=0.1)+
  #geom_violin(data=plot_dat, aes(x=experiment_type, y=scale_effect, fill=experiment_type), alpha=0.3)+
  stat_summary(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, colour=experiment_type), 
               fun.data="mean_cl_boot", 
               geom = "crossbar", 
               width=0.75, 
               alpha=0.5)+
  scale_fill_manual(values=c("#440154", "yellow3"))+
  scale_colour_manual(values=c("#440154", "yellow3"))+
  geom_hline(yintercept = 0, linetype='dashed')+
  facet_wrap(~size,  labeller = labeller(size = size.labs), scales="free")+
  geom_point(data = vline_size,
             aes(y = coef, x=experiment_type, fill=experiment_type),colour="black",pch=21, size=5)+
  theme_linedraw()+
  ylab("Effect size")+
  xlab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        #axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_label(data=vline_size, aes(x=experiment_type, label=n, y=-3.5, fill=experiment_type),
             alpha=0.5, size=6)+
  ggtitle("Size class")+
  geom_text(data=pp3S_labels, aes(x=x, y=y,  label=labels), colour="black", size=10)


pp4S_labels <- data.frame(trait_directional=c("positive", "negative"), 
                         labels=c("I", "J"), 
                         x=c(0.5, 0.5),
                         y=c(16.5, 34))

vline_trait <- dat_ES_final_2 %>%
  group_by(experiment_type, trait_directional) %>%
  summarise(mean=mean(scale_effect), 
            median=median(scale_effect), 
            n=n())%>%
  cbind(coef=c(-0.186, 
               -0.186+-0.034, 
               -0.186+-0.073, 
               -0.186+-0.034+-0.073+-0.188))

pp4S <- ggplot()+
  geom_point(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, colour=experiment_type), alpha=0.1)+
  #geom_violin(data=plot_dat, aes(x=experiment_type, y=scale_effect, fill=experiment_type), alpha=0.3)+
  stat_summary(data=dat_ES_final_2, aes(x=experiment_type, y=scale_effect, colour=experiment_type), 
               fun.data="mean_cl_boot", 
               geom = "crossbar", 
               width=0.75, 
               alpha=0.5)+
  scale_fill_manual(values=c("#440154", "yellow3"))+
  scale_colour_manual(values=c("#440154", "yellow3"))+
  geom_hline(yintercept = 0, linetype='dashed')+
  facet_wrap(~trait_directional, scales="free")+
  geom_point(data = vline_trait,
             aes(y = coef, x=experiment_type, fill=experiment_type),colour="black",pch=21, size=5)+
  theme_linedraw()+
  ylab("Effect size")+
  xlab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        #axis.text.x = element_text( size = 20),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  geom_label(data=vline_trait, aes(x=experiment_type, label=n, y=-10, fill=experiment_type),
             alpha=0.5, size=6)+
  ggtitle("Trait directionality")+
  geom_text(data=pp4S_labels, aes(x=x, y=y,  label=labels), colour="black", size=10)

pp2S + pp1S + pp3S + pp4S


########fig s5

mean_labels_S <- data.frame(experiment_type = c("acclimation", "acute"), labels=c("A", "B"))

mean_temp<- dat_ES_final_2 %>%
  #filter(between(scale_effect, -2.5, 2.5)) %>%
  ggplot(aes(x=I(mean_temp_constant-mean(mean_temp_constant)), y=scale_effect, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.3, size=2)+
  scale_colour_manual(values =c("#440154", "yellow3"))+
  #geom_smooth(method="lm")+
  theme_bw()+
  xlab("Standardized mean temperature (°C)")+
  ylab("Effect size")+
  theme_linedraw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))+
  geom_abline(aes(intercept=i, slope=s, colour=experiment_type), 
              data=data.frame(experiment_type=c("acute","acclimation"), i=c(-0.186, -0.186 + -0.073), s=c(-0.061,-0.061+ 0.0572782))) +
  facet_wrap(~experiment_type)+
  geom_text(data=mean_labels_S, aes(x=-17, y=30,  label=labels), colour="black", size=10)



flux_labels_S <- data.frame(experiment_type = c("acclimation", "acute"), labels=c("C", "D"))


flux_range <- dat_ES_final_2 %>%
  #filter(between(scale_effect, -2.5, 2.5)) %>%
  ggplot(aes(x=I(flux_range-mean(flux_range)), y=scale_effect, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.3, size=2)+
  scale_colour_manual(values =c("#440154", "yellow3"))+
  theme_linedraw()+
  #geom_smooth(method="lm")+
  xlab("Standardized fluctuation magnitude (°C)")+
  ylab("Effect size")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))+
  geom_abline(aes(intercept=i, slope=s, colour=experiment_type), 
              data=data.frame(experiment_type=c("acute","acclimation"), i=c(-0.186, -0.186 + -0.073), s=c(-0.010,-0.0108231+-0.0046560))) +
  facet_wrap(~experiment_type)+
  geom_text(data=flux_labels_S, aes(x=-8, y=30,  label=labels), colour="black", size=10)



duration_labels_S <- data.frame(experiment_type = c("acclimation", "acute"), labels=c("G", "H"))

duration<- dat_ES_final_2 %>%
  #filter(between(scale_effect, -2.5, 2.5)) %>%
  ggplot(aes(x=log10(duration_standard), y=scale_effect, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.3, size=2)+
  scale_colour_manual(values =c("#440154", "yellow3"))+
  theme_linedraw()+
  #geom_smooth(method="lm")+
  xlab("Duration (log(days))")+
  ylab("Effect size")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))+
  geom_abline(aes(intercept=i, slope=s, colour=experiment_type), 
              data=data.frame(experiment_type=c("acute","acclimation"), i=c(-0.186, -0.186 + -0.073), s=c(-0.052,-0.052))) +
  geom_text(data=duration_labels_S, aes(x=-0.275, y=30,  label=labels), colour="black", size=10)+
  facet_wrap(~experiment_type)



secondary_labels_S <- data.frame(experiment_type = c("acclimation", "acute"), labels=c("E", "F"))

secondary_temp <- dat_ES_final_2 %>%
  #filter(between(scale_effect, -2.5, 2.5)) %>%
  ggplot(aes(x=I(secondary_temp-mean(secondary_temp)), y=scale_effect, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.3, size=2)+
  scale_colour_manual(values =c("#440154", "yellow3"))+
  theme_linedraw()+
  #geom_smooth(method="lm")+
  xlab("Standardized secondary temperature (°C)")+
  ylab("Effect size")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text( size = 20 ),
        legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20))+
  geom_abline(aes(intercept=i, slope=s, colour=experiment_type), 
              data=data.frame(experiment_type=c("acute","acclimation"), i=c(-0.186, -0.186 + -0.073), s=c(0.0002018,0.0002018))) +
  geom_text(data=secondary_labels_S, aes(x=-33, y=30,  label=labels), colour="black", size=10)+
  #data=data.frame(experiment_type =c("acute", "acclimation"), y=3.35, x=10, label=c("*")), size=20)+
  facet_wrap(~experiment_type)


mean_temp + flux_range +  secondary_temp + duration 




#######################################################
##I2 calculations from best models 
#phylo level heterogeneity
full_mode$sigma2[1]/sum(full_mode$sigma2[1], full_mode$sigma2[2], full_mode$sigma2[3]) #0.9637687
full_modj$sigma2[1]/sum(full_modj$sigma2[1], full_modj$sigma2[2], full_modj$sigma2[3]) #0.9580754
#study level heterogeneity
full_mode$sigma2[2]/sum(full_mode$sigma2[1], full_mode$sigma2[2], full_mode$sigma2[3]) #0.01272126
full_modj$sigma2[2]/sum(full_modj$sigma2[1], full_modj$sigma2[2], full_modj$sigma2[3]) #0.01668716
#response heterogeneity 
full_mode$sigma2[3]/sum(full_mode$sigma2[1], full_mode$sigma2[2], full_mode$sigma2[3]) #0.02351002
full_modj$sigma2[3]/sum(full_modj$sigma2[1], full_modj$sigma2[2], full_modj$sigma2[3]) #0.02523745

#calculating fail safe number
fsn(yi, vi, data=dat_ES_final_2,type="Rosenberg") #260582

#correlation between age and size
dat_cor <- dat_ES_final_2 %>%
  mutate(exp_age=as.numeric(exp_age), 
         size=as.numeric(size))

res0 <- cor.test(dat_cor$exp_age, dat_cor$size, 
                 method = "pearson") #-0.3950939

#model coefficient output
model_df <- as.data.frame(coef(final_model)) %>% tibble::rownames_to_column("covariates") %>%
  mutate_if(is.numeric, ~round(., 3))
write_csv(model_df, "model_df.csv")

#frequency table
Var.tbl <- xtabs(~experiment_type+size+exp_age+trait_directional, dat_ES_final_2)
Var.dbf <- as.data.frame.table(Var.tbl)

write_csv(Var.dbf, "var_dbf.csv")

#category table 
tables1 <- count(dat_ES_final_2, trait_directional, resp_def) %>%
  as.data.frame()
write_csv(tables1, "tables1.csv")

#prisma diagram 
prisma <- read_csv("metadata/litsearch_subgroups14sept21 copy.csv")
count(prisma, notes)





stand_dat<- dat_ES_final_2 %>%
  mutate(flux_range2 = round(I(flux_range-mean(flux_range)), 1), 
         mean_temp2= round(I(mean_temp_constant-mean(mean_temp_constant)), 1))

upper <- stand_dat %>%
  filter(mean_temp2 > 0)

lower <-  stand_dat %>%
  filter(mean_temp2 < 0)


upper %>%
  #group_by(experiment_type) %>%
  summarise(mean_yi = mean(scale_effect))

lower %>%
  #group_by(experiment_type) %>%
  summarise(mean_yi = mean(scale_effect))








