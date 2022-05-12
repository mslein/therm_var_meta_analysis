#load data
pacman::p_load(tidyverse, metafor, rotl, taxadb, taxize, ape)
acclim <- read_csv("data/acclim_cleaned.csv")%>%
  mutate(environment=as.factor(environment), size=as.factor(size), 
         exp_age=as.factor(exp_age), tpc_zone = as.factor(tpc_zone), 
         experiment_type = "acclimation", 
         mean_temp_constant = mean_temp_reared, 
         metric = tpc_zone) %>%
  select(response_id, species, study_id, org_level, 
         mean_temp_constant, exp_age, size, yi, vi,metric, experiment_type, 
         flux_range, environment, genus)

acute <- read_csv("data/acute_cleaned.csv")%>%
  mutate(metric= as.factor(metric), environment=as.factor(environment), size=as.factor(size), 
         exp_age=as.factor(exp_age), org_level=as.factor(org_level), 
         experiment_type = "acute") %>%
  select(response_id, species, study_id, org_level, 
         mean_temp_constant, exp_age, size, yi, vi,metric, experiment_type, 
         flux_range, environment, genus)

full <- rbind(acute, acclim) %>%
  mutate() %>% 
  unite(genus_species, c(genus, species), sep = " ", remove = FALSE)

#first attempt at running a model with all the data together
full_best <- metafor::rma.mv(yi, vi, data=full, 
                     mods = ~ I(flux_range-mean(flux_range)) * I(mean_temp_constant-mean(mean_temp_constant))
                     +environment +size  +exp_age +org_level 
                     + experiment_type +metric,
                     random = ~1 | study_id/species/response_id,
                     method="REML") 

#42 genuses
length(unique(full$genus))

#trying to get family for each genus from rotl
family_taxa <- tax_name(query = unique(full$genus), get = "family", db = "ncbi") 
#still missing 5 species here at the family level
#Culex -- Culicidae (wiki)
#Drosophila -- Drosophilidae (wiki)
#Calyptocephalella -- Calyptocephalellidae (wiki)
#Isotopenola -- Isotomidae (original article)
#Mucrosomia -- Isotomidae (original article)
family_taxa_df <- family_taxa %>%
  as.data.frame() %>%
  drop_na()
#adding the missing families + querying rotl
family_taxa_missing_df <- data.frame(db=c("wiki", "wiki", "wiki", "hoskins2020", "hoskins2020"), 
                                     query= c("Culex", "Drosophila", "Calyptocephalella", "Isotopenola", 
                                              "Mucrosomia"), 
                                     family=c("Culicidae", "Drosophilidae", "Calyptocephalellidae", 
                                              "Isotomidae", "Isotomidae"))
#combining for the full set of families
family_total <- rbind(family_taxa_df, family_taxa_missing_df) #32 unique families
#rotl computation to get phylogenetic tree
taxa_family <- tnrs_match_names(unique(family_total$family))
#it would appear there are ott_ids for each of the families
#wondering if the tree is being pruned because there is no way to compute phylogenies between all of these?
tr_family <- tol_induced_subtree(ott_id(taxa_family)[is_in_tree(ott_id(taxa_family))])
plot(tr_family, show.tip.label = TRUE)
tree_family <- compute.brlen(tr_family)
cor_family <- vcv(tree_family, cor = T)


