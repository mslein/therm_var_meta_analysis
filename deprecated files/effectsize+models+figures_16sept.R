#####################################################
## Manuscript coding, models, and figures ###########
#####################################################
##clean workspace
# Clean workspace
rm(list = ls()) 


##loading libraries 
pacman::p_load(tidyverse, metafor, patchwork, rotl, MuMIn, viridis, ape, wesanderson, glmulti, MuMIn, metaAidR, corrplot)
#loading in the raw data file
dat_0<- read_csv("data/rawdata29jun22.csv") %>%
  select(-add_covariate, -add_covariate_value, -extra_covariate, -extra_value, -larger_group, -flux_pattern) 
##wrangling calculations for acclimation data
#make standardized SD columns based on N and SE vs SD
#need to double check the 3.92 as being appropriate or how i would get larger numbers from the t-dist?
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


dat_ES <-escalc(measure="SMD", m1i=flux_resp, m2i=constant_resp, 
                       sd1i=`SD_flux`, sd2i= `SD_constant`, n1i=flux_samp, n2i=constant_samp, 
                       data=dat) 



dat_ES2<-dat_ES %>% 
  drop_na() %>% #4% data loss
  unite(original_name, c(genus, species), sep = " ", remove = FALSE) %>%
  mutate(original_name = tolower(original_name))
#selecting only organismal level responses, adult, and juvenile organisms, and species that did have phylogenies
dat_ES_final<- dat_ES2 %>% #kramarz seems to be quite the outlier in the qqplot
  filter(org_level %in% c("0")) %>%
  filter(size %in% c("1", "2")) %>%
  filter(!original_name %in% c("phanaeus vindex", 
                              "cryptopygus sp.", 
                              "metacatharsius opacus")) %>%
  mutate(ecosystem = case_when(environment == "terrestrial" ~ "terrestrial", 
         TRUE ~ "aquatic/marine")) %>%
  mutate(experiment_type = factor(experiment_type, ordered = FALSE), 
         exp_age=factor(exp_age, ordered = FALSE), 
         size=factor(size, ordered = FALSE), 
         trait_directionality=factor(trait_directionality, ordered = FALSE), 
         ecosystem=factor(ecosystem, ordered = FALSE))


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

  

install.packages("devtools")
library(devtools)
install_github("daniel1noble/metaAidR")

# Error calc-residual variance at the observation level as metafor does not add this by default
dat_ES_final_2$obs <- 1:dim(dat_ES_final_2)[1]

# Make shared control matrix for ROM dataset
dat_ES_final_2$sc_cluster <- interaction(dat_ES_final_2$study_id, 
                                         dat_ES_final_2$shared_control)

# Create the Shared control V matrix
V <- metaAidR::make_VCV_matrix(dat_ES_final_2, "vi", cluster= "sc_cluster", 
                               type = "vcv", rho = 0.5)





####candidate models
int_mod <- rma.mv(yi=yi, V=V, 
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)

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



#full_modc prevails -- no secondary temperature
# hypothesis 1

#full model with interaction between experiment_type and the mods
full_mode<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                  +I(flux_range-mean(flux_range))*experiment_type
                  +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                  +duration_standard
                  +exp_age*experiment_type
                  +size*experiment_type
                  +ecosystem
                  +trait_directionality*experiment_type
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor),
                  method="ML", data=dat_ES_final_2)

#no experiment type 
full_modf<-rma.mv(yi=yi, V=V, 
                 mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                 +duration_standard
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
                  +exp_age*experiment_type
                  +size*experiment_type
                  +ecosystem
                  +trait_directionality
                  +sqrt_ni_1,
                  random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                  R = list(phylo=cor), test="t",
                  method="ML", data=dat_ES_final_2)



model.sel(full_mode, full_modf, full_modg, full_modh, full_modi, full_modj, full_modk, full_modl)


#full mod h prevails by more than 2 AIC when I ran two different versions with and without the interaction in hypothesis 2 


#h2 and h3 without interaction flux_range (based on full_modh)
# hypothesis 2 
# no mean temp
full_modm<- rma.mv(yi=yi, V=V, 
                   mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                   +I(flux_range-mean(flux_range))
                   +duration_standard
                   +exp_age*experiment_type
                   +size*experiment_type
                   +ecosystem
                   +trait_directionality*experiment_type
                   +sqrt_ni_1,
                   random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                   R = list(phylo=cor), test="t",
                   method="ML", data=dat_ES_final_2)
#no flux_range
full_modn <- rma.mv(yi=yi, V=V, 
                    mods = ~ I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                    +duration_standard
                    +exp_age*experiment_type
                    +size*experiment_type
                    +ecosystem
                    +trait_directionality*experiment_type
                    +sqrt_ni_1,
                    random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                    R = list(phylo=cor), test="t",
                    method="ML", data=dat_ES_final_2)


#no interaction between mean_temp and flux range
full_modo <-rma.mv(yi=yi, V=V, 
                   mods = ~ I(flux_range-mean(flux_range))
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
#no flux range, mean temp, or interaction
full_modp <- rma.mv(yi=yi, V=V, 
                    mods = ~ duration_standard
                    +exp_age*experiment_type
                    +size*experiment_type
                    +ecosystem
                    +trait_directionality*experiment_type
                    +sqrt_ni_1,
                    random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                    R = list(phylo=cor), test="t",
                    method="ML", data=dat_ES_final_2)

model.sel(full_modh, full_modm, full_modn, full_modo, full_modp)



#hypothesis 3 - full_modh as base from h2
#no interaction with age
full_modq <- rma.mv(yi=yi, V=V, 
                    mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                    +I(flux_range-mean(flux_range))
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
full_modr <- rma.mv(yi=yi, V=V, 
                    mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                    +I(flux_range-mean(flux_range))
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
full_mods <- rma.mv(yi=yi, V=V, 
                    mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                    +I(flux_range-mean(flux_range))
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
full_modt <- rma.mv(yi=yi, V=V, 
                    mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                    +I(flux_range-mean(flux_range))
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
full_modu <- rma.mv(yi=yi, V=V, 
                   mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                   +I(flux_range-mean(flux_range))
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
#no age or size
full_modv <- rma.mv(yi=yi, V=V, 
                    mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                    +I(flux_range-mean(flux_range))
                    +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                    +duration_standard
                    +ecosystem
                    +trait_directionality*experiment_type
                    +sqrt_ni_1,
                    random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                    R = list(phylo=cor), test="t",
                    method="ML", data=dat_ES_final_2)
#no interaction with age or size 
full_modw <- rma.mv(yi=yi, V=V, 
                    mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                    +I(flux_range-mean(flux_range))
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
full_modx <- rma.mv(yi=yi, V=V, 
                    mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                    +I(flux_range-mean(flux_range))
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
full_mody <- rma.mv(yi=yi, V=V, 
                    mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                    +I(flux_range-mean(flux_range))
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

full_modz <- rma.mv(yi=yi, V=V, 
                    mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                    +I(flux_range-mean(flux_range))
                    +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                    +duration_standard
                    +ecosystem
                    +sqrt_ni_1,
                    random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                    R = list(phylo=cor), test="t",
                    method="ML", data=dat_ES_final_2)

model.sel(full_modh, full_modq, full_modr, full_mods, full_modt, full_modu, full_modv, full_modw, full_modx, full_mody, full_modz)


#h2 and h3 with interaction between flux range and experiment (based on full_mode as base)
# hypothesis 2 
# no mean temp
full_modaa<-rma.mv(yi=yi, V=V, 
                  mods = ~ I(flux_range-mean(flux_range))*experiment_type
                  +duration_standard
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
       +size*experiment_type
       +ecosystem
       +trait_directionality*experiment_type
       +sqrt_ni_1,
       random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
       R = list(phylo=cor), test="t",
       method="ML", data=dat_ES_final_2)

model.sel(full_mode,
          full_modaa, full_modbb, full_modcc, full_moddd)



#hypothesis 3
#no interaction with age (based on full_mode)
full_modee <- rma.mv(yi=yi, V=V, 
                    mods = ~ I(flux_range-mean(flux_range))*I(mean_temp_constant-mean(mean_temp_constant))
                    +I(flux_range-mean(flux_range))*experiment_type
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
                    +I(mean_temp_constant-mean(mean_temp_constant))*experiment_type
                    +duration_standard
                    +ecosystem
                    +sqrt_ni_1,
                    random = list( ~1 | phylo, ~1 | study_id, ~1 | response_id),
                    R = list(phylo=cor), test="t",
                    method="ML", data=dat_ES_final_2)


model.sel(full_mode, full_modee, full_modff, full_modgg, full_modhh, full_modii, full_modjj, 
          full_modkk, full_modll, full_modmm, full_modnn)


#full_mode is best in this version


model.sel(full_mode, full_modh)


eval(metafor:::.MuMIn)
model_avg<-model.avg(model.sel(full_mode, full_modh))

summary(model_avg)


#getting confidence intervals 
confint(model_avg, full = T)


#########################
#########Figures#########
#########################

mod_coef<- as.data.frame(coef(summary(model_avg))) %>%
  tibble::rownames_to_column("covariates") %>%
  rename(SMD = "estimate") %>%
  mutate(type = case_when(covariates %in% c("as.factor(exp_age)1", "as.factor(exp_age)2") ~ "Age", 
                          covariates %in% c("as.factor(size)2") ~ "Body size",
                          covariates %in% c("trait_directionalitypositive") ~ "Trait directionality", 
                          covariates %in% c("experiment_typeacute") ~ "acute", 
                          covariates %in% c("I(flux_range - mean(flux_range))", 
                                            "I(mean_temp_constant - mean(mean_temp_constant))", 
                                            "I(secondary_temp - mean(secondary_temp))", 
                                            "I(secondary_temp - mean(secondary_temp))") ~ "Temperature", 
                          covariates %in% c("ecosystemterrestrial") ~ "Ecosystem", 
                          covariates %in% c("vi") ~ "Variance",
                          covariates %in% c("duration_standard") ~ "Duration", 
                          covariates %in% c("intrcpt") ~ "Intercept")) %>%
  mutate(estimate_add = case_when(type %in% c("Intercept", "Duration", "Temperature") ~ SMD, 
                                  TRUE ~ SMD + 1.7464),
         ciub_add = case_when(type %in% c("Intercept", "Duration", "Temperature") ~ ci.ub, 
                                  TRUE ~ ci.ub + 1.7464), 
         cilb_add = case_when(type %in% c("Intercept", "Duration", "Temperature") ~ ci.lb, 
                                  TRUE ~ ci.lb + 1.7464))

exp_mod <- mod_coef %>%
  filter(type %in% c("acute"))

dat_ES_final_2_exp <- dat_ES_final_2 %>%
  filter(experiment_type == "acute")

a <- ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  #geom_violin(data=dat_ES_final_2_exp, aes(x=experiment_type, y=yi), width = 0.15, fill="#440154")+
  geom_jitter(data=dat_ES_final_2_exp, aes(x=experiment_type, y=yi), width = 0.15, alpha=0.1, colour="#440154")+
  geom_pointrange(data=exp_mod, aes(x=type, y=estimate_add, ymin=cilb_add, ymax=ciub_add), colour="red")+
  coord_flip()+
  theme_bw()+
  xlab("Experiment")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))
  #geom_label(aes(label=n),nudge_y=1.5, size=8)

age_mod <- mod_coef %>%
  filter(type %in% c("Age")) %>%
  mutate(type2 = case_when(covariates %in% c("as.factor(exp_age)2") ~ "2",
                           covariates %in% c("as.factor(exp_age)1") ~ "1"))

dat_ES_final_2_age <- dat_ES_final_2 %>%
  mutate(yi_plus=yi+1) %>%
  filter(exp_age %in% c("1", "2"))

b <- ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  #geom_violin(data=dat_ES_final_2_age, aes(x=as.factor(exp_age), y=yi), width = 0.15, fill="#414487")+
  geom_jitter(data=dat_ES_final_2_age, aes(x=as.factor(exp_age), y=yi), width = 0.15, alpha=0.1, colour="#414487")+
  geom_pointrange(data=age_mod, aes(x=type2, y=estimate_add, ymin=cilb_add, ymax=ciub_add), colour="red")+
  coord_flip()+
  theme_bw()+
  xlab("Age")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

size_mod <- mod_coef %>%
  filter(covariates == "as.factor(size)2") %>%
  mutate(type2 = case_when(covariates %in% c("as.factor(size)2") ~ "2"))

dat_ES_final_2_size <- dat_ES_final_2 %>%
  mutate(yi_plus=yi+1) %>%
  filter(size %in% c("2"))


c <- ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  #geom_violin(data=dat_ES_final_2_size, aes(x=as.factor(size), y=yi), width = 0.15, fill="#2a788e")+
  geom_jitter(data=dat_ES_final_2_size, aes(x=as.factor(size), y=yi), width = 0.15, alpha=0.1, colour="#2a788e")+
  geom_pointrange(data=size_mod, aes(x=type2, y=estimate_add, ymin=cilb_add, ymax=ciub_add), colour="red")+
  coord_flip()+
  theme_bw()+
  xlab("Size")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))


ecosystem_mod <- mod_coef %>%
  filter(covariates == "ecosystemterrestrial") %>%
  mutate(type2 = case_when(covariates %in% c("ecosystemterrestrial") ~ "terrestrial"))

dat_ES_final_2_eco <- dat_ES_final_2 %>%
  mutate(yi_plus=yi+1) %>%
  filter(ecosystem %in% c("terrestrial"))

d <- ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  #geom_violin(data=dat_ES_final_2_eco, aes(x=ecosystem, y=yi), width = 0.15, fill="#22a884")+
  geom_jitter(data=dat_ES_final_2_eco, aes(x=ecosystem, y=yi), width = 0.15, alpha=0.1, colour="#22a884")+
  geom_pointrange(data=ecosystem_mod, aes(x=type2, y=estimate_add, ymin=cilb_add, ymax=ciub_add), colour="red")+
  coord_flip()+
  theme_bw()+
  xlab("Ecosystem")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

trait_mod <- mod_coef %>%
  filter(covariates == "trait_directionalitypositive") %>%
  mutate(type2 = case_when(covariates %in% c("trait_directionalitypositive") ~ "positive"))

dat_ES_final_2_trait <- dat_ES_final_2 %>%
  mutate(yi_plus=yi+1) %>%
  filter(trait_directionality %in% c("positive"))

e <- ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  #geom_violin(data=dat_ES_final_2_trait, aes(x=trait_directionality, y=yi), width = 0.15, fill="#7ad151")+
  geom_jitter(data=dat_ES_final_2_trait, aes(x=trait_directionality, y=yi), width = 0.15, alpha=0.1, colour="#7ad151")+
  geom_pointrange(data=trait_mod, aes(x=type2, y=estimate_add, ymin=cilb_add, ymax=ciub_add), colour="red")+
  coord_flip()+
  theme_bw()+
  xlab("Trait")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

overall_mod <- mod_coef %>%
  filter(covariates == "intrcpt") %>%
  mutate(type2 = case_when(covariates %in% c("intrcpt") ~ "intercept"))

dat_ES_final_2_intr <- dat_ES_final_2 %>%
  mutate(yi_plus=yi+1, 
         intercept="intercept") 

f <- ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  #geom_violin(data=dat_ES_final_2_intr, aes(x=intercept, y=yi), width = 0.15, fill="#fde725")+
  geom_jitter(data=dat_ES_final_2_intr, aes(x=intercept, y=yi), width = 0.15, alpha=0.1, colour="#fde725")+
  geom_pointrange(data=overall_mod, aes(x=type2, y=estimate_add, ymin=cilb_add, ymax=ciub_add), colour="red")+
  coord_flip()+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

f/a/b/c/d/e



#interaction model
redo_coef<- as.data.frame(coef(summary(redo_mod))) %>%
  tibble::rownames_to_column("covariates") %>%
  rename(SMD = "estimate") %>%
  mutate(type = case_when(covariates %in% c("intrcpt", "I(flux_range - mean(flux_range))", 
                                                    "I(mean_temp_constant - mean(mean_temp_constant))",
                                                    "I(secondary_temp - mean(secondary_temp))", 
                                                    "duration_standard", 
                                                    "vi", "I(flux_range - mean(flux_range)):experiment_typeacute", 
                                                    "experiment_typeacute:I(mean_temp_constant - mean(mean_temp_constant))", 
                                                    "experiment_typeacute:I(secondary_temp - mean(secondary_temp))", 
                                                    "experiment_typeacute:duration_standard", 
                                                    "experiment_typeacute:vi") ~ "continuous", 
                                  TRUE ~ "categorical"),
         smd_add = case_when(type == "categorical" ~ SMD+1.3743852210, 
                              TRUE ~ SMD), 
         cilb_add = case_when(type == "categorical" ~ ci.lb +1.3743852210, 
         TRUE ~ ci.lb), 
         ciub_add = case_when(type == "categorical" ~ ci.ub +1.3743852210, 
                              TRUE ~ ci.ub))

redo_coef_exp <- redo_coef %>%
  filter(covariates == "experiment_typeacute") %>%
  mutate(type2 = case_when(covariates %in% c("experiment_typeacute") ~ "acute"))

redo_acute <- dat_ES_final_2 %>%
  filter(experiment_type=="acute")



ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  #geom_violin(data=redo_acute, aes(x=experiment_type, y=yi), width = 0.15, fill="#440154")+
  geom_jitter(data=redo_acute, aes(x=experiment_type, y=yi), width = 0.15, alpha=0.1, colour="#440154")+
  geom_pointrange(data=redo_coef_exp, aes(x=type2, y=smd_add, ymin=cilb_add, ymax=ciub_add), colour="red")+
  geom_pointrange(data=redo_coef_exp, aes(x=type2, y=SMD, ymin=ci.lb, ymax=ci.ub), colour="green")+
  coord_flip()+
  theme_bw()+
  xlab("Experiment")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))



ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_violin(data=dat_ES_final_2, aes(x=as.factor(exp_age), y=yi, fill=experiment_type), width = 0.15)+
  #geom_jitter(data=dat_ES_final_2, aes(x=as.factor(exp_age), y=yi, shape=experiment_type), width = 0.15, alpha=0.1, colour="#414487")+
  geom_pointrange(data=age_mod, aes(x=type2, y=estimate_add, ymin=cilb_add, ymax=ciub_add), colour="red")+
  scale_fill_viridis_d()+
  coord_flip()+
  theme_bw()+
  xlab("Age")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))


ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_violin(data=dat_ES_final_2, aes(x=as.factor(size), y=yi, fill=experiment_type), width = 0.15)+
  #geom_jitter(data=dat_ES_final_2, aes(x=as.factor(size), y=yi, fill=experiment_type), width = 0.15, alpha=0.1)+
  geom_pointrange(data=size_mod, aes(x=type2, y=estimate_add, ymin=cilb_add, ymax=ciub_add), colour="red")+
  scale_fill_viridis_d()+
  coord_flip()+
  theme_bw()+
  xlab("Size")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_violin(data=dat_ES_final_2, aes(x=ecosystem, y=yi, fill=experiment_type), width = 0.15)+
  #geom_jitter(data=dat_ES_final_2, aes(x=ecosystem, y=yi, shape=experiment_type), width = 0.15, alpha=0.1, colour="#22a884")+
  geom_pointrange(data=ecosystem_mod, aes(x=type2, y=estimate_add, ymin=cilb_add, ymax=ciub_add), colour="red")+
  scale_fill_viridis_d()+
  coord_flip()+
  theme_bw()+
  xlab("Ecosystem")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_violin(data=dat_ES_final_2, aes(x=trait_directionality, y=yi, fill=experiment_type), width = 0.15)+
  #geom_jitter(data=dat_ES_final_2, aes(x=trait_directionality, y=yi, shape=experiment_type), width = 0.15, alpha=0.1, colour="#7ad151")+
  geom_pointrange(data=trait_mod, aes(x=type2, y=estimate_add, ymin=cilb_add, ymax=ciub_add), colour="red")+
  scale_fill_viridis_d()+
  coord_flip()+
  theme_bw()+
  xlab("Trait")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

ggplot()+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_violin(data=dat_ES_final_2_intr, aes(x=intercept, y=yi, fill=experiment_type), width = 0.15)+
  #geom_jitter(data=dat_ES_final_2_intr, aes(x=intercept, y=yi), width = 0.15, alpha=0.1, colour="#fde725")+
  geom_pointrange(data=overall_mod, aes(x=type2, y=estimate_add, ymin=cilb_add, ymax=ciub_add), colour="red")+
  coord_flip()+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))



p0 <-dat_ES_final_2%>%
  group_by(experiment_type) %>%
  summarise(mean=mean(yi), 
            sd=sd(yi), 
            n=n()) %>%
  mutate(margin = qt(0.975,df=n-1)*sd/sqrt(n), 
         lower = mean - margin, 
         upper = mean + margin, 
         se = sd/sqrt(n)) %>%
  ggplot(aes(x=experiment_type, y=mean, ymax=mean+se, ymin=mean-se, colour=experiment_type))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(size=7, position=position_dodge(width = .75))+
  geom_linerange(position=position_dodge(width = .75), size=1)+
  coord_flip()+
  theme_bw()+
  xlab("Experiment type")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))+
  geom_label(aes(label=n),nudge_y=1.5, size=8)

p1 <-dat_ES_final_2%>%
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
  geom_point(size=7, position=position_dodge(width = .75))+
  geom_linerange(position=position_dodge(width = .75), size=1)+
  coord_flip()+
  theme_bw()+
  #xlab("Age")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))+
  geom_label(aes(label=n),nudge_y=1.5, size=8)

p2 <-dat_ES_final_2%>%
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
  geom_point(size=7, position=position_dodge(width = .75))+
  geom_linerange(position=position_dodge(width = .75), size=1)+
  coord_flip()+
  theme_bw()+
  #xlab("Age")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))+
  geom_label(aes(label=n),nudge_y=1.5, size=8)


p3 <-dat_ES_final_2%>%
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
  geom_point(size=7, position=position_dodge(width = .75))+
  geom_linerange(position=position_dodge(width = .75), size=1)+
  coord_flip()+
  theme_bw()+
  #xlab("Age")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))+
  geom_label(aes(label=n),nudge_y=1.5, size=8)

p4 <-dat_ES_final_2%>%
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
  geom_point(size=7, position=position_dodge(width = .75))+
  geom_linerange(position=position_dodge(width = .75), size=1)+
  coord_flip()+
  theme_bw()+
  #xlab("Age")+
  ylab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))+
  geom_label(aes(label=n),nudge_y=1.5, size=8)

p0 / p1 / p2 / p3 / p4







library(patchwork)


#Figure 4
#1681 effect sizes vs the 1703 in the model
mean_temp<- dat_ES_final_2 %>%
  filter(between(yi, -5, 5)) %>%
  ggplot(aes(x=mean_temp_constant, y=yi, colour=mean_temp_constant))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.5, size=2)+
  scale_colour_viridis(option="D")+
  #geom_smooth(method="lm")+
  theme_bw()+
  xlab("Mean temperature (°C)")+
  ylab("Effect size")+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))+
  geom_abline(intercept = 0.4814732, slope = 0.0086691, color="black", size=1)+
  geom_abline(intercept = 0.4814732, slope = -0.0065049, color="red", size=1)
  #geom_label(x=28, y=9, aes(label="n=1620"), size=8, colour="black")


flux_range <- dat_ES_final_2 %>%
  filter(between(yi, -5, 5)) %>%
  ggplot(aes(x=flux_range, y=yi, colour=flux_range))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.5, size=2)+
  scale_colour_viridis(option="D")+
  #geom_smooth(method="lm")+
  theme_bw()+
  xlab("Fluctuation magnitude (°C)")+
  ylab("Effect size")+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))+
  geom_abline(intercept = 0.4814732, slope = -0.0109933, color="black", size=1)+
  geom_abline(intercept = 0.4814732, slope = -0.0032488, color="red", size=1)
  #geom_label(x=28, y=9, aes(label="n=1620"), size=8, colour="black")+
  #facet_wrap(~experiment_type)
  #geom_label(x=16, y=9, aes(label="n=1620"), size=8, colour="black")

mean_and_flux <- dat_ES_final_2 %>%
  filter(between(yi, -5, 5)) %>%
  ggplot(aes(x=mean_temp_constant, y=yi, colour=flux_range))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_jitter(alpha=0.5, size=2)+
  scale_colour_viridis(option="D")+
  #facet_wrap(~flux_range)+
  #geom_abline(intercept = 0.4814732, slope = 0.0089653, color="black", size=1)+
  #geom_smooth(method="lm")+
  theme_bw()+
  xlab("Mean temperature (°C)")+
  ylab("SMD")+
  #facet_wrap(~flux_range)
  #theme(legend.position="none")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))

#Figure 5
duration<- dat_ES_final_2 %>%
  filter(between(yi, -10, 10)) %>%
  #filter(between(duration, 0, 200)) %>%
  ggplot(aes(x=duration_standard, y=yi, colour=duration_standard))+
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_point(alpha=0.5, size=2)+
  scale_colour_viridis(option="D")+
  theme_bw()+
  xlab("Time (days)")+
  ylab("Effect size")+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))+
  geom_abline(intercept = 0.4814732, slope = -0.0282509 , color="black", size=1)
  #geom_label(x=150, y=9, aes(label="n=1620"), size=8, colour="black")

mean_temp + flux_range 

mean_and_flux 

duration

t <- dat_ES_final_2 %>%
  filter(between(yi, -10, 10))
length(t$yi)

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
                         method="ML", data=dat_ES_final_2)
funnel(dummy_full_mod)

##histogram
dat_ES_final_2 %>%
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

#calculating fail safe number
fsn(yi, vi, data=dat_ES_final_2,type="Rosenberg") #238691

#correlation between age and size
#-0.42
res0 <- cor.test(as.numeric(dat_ES_final_2$exp_age), as.numeric(dat_ES_final_2$size), 
                method = "pearson")
