pacman::p_load(tidyverse)
prisma <- read_csv("metadata/litsearch_subgroups14sept21 copy.csv") %>%
  mutate(notes_tidy = case_when(notes %in% c("background", "review", 
                                                    "review/synthesis", "book") ~ "reviews/background", 
                                       notes %in% c("modeling", "modelling") ~ "modeling", 
                                       notes %in% c("no constant", 
                                                    "no constant/flux", 
                                                    "no constant/fux", 
                                                    "no flux pattern", 
                                                    "no flux treatment", 
                                                    "no flux trt") ~ "no flux/constant",
                                       notes %in% c("no error", 
                                                    "no error reported") ~ "no error", 
                                       notes %in% c("non-biologically relevant", 
                                                    "not extractable") ~ "not relevant",
                                       notes %in% c("non-lab") ~ "uncontrolled var", 
                                      TRUE  ~ "greater diff"))
tally <- count(prisma, notes_tidy) 
sum(tally$n)

