args <- commandArgs(trailingOnly=TRUE)

library(tidyverse)

all_sample_ped_file <- args[1]
index_ped_file <- args[2]
index_sample_file <- args[3]

#"Z:/resources/OGLsample/genome.2022-08.ped"

all_sample <- read_tsv(all_sample_ped_file, col_names = FALSE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert() 

index <- filter(all_sample, X6 == "2") %>% 
  separate(X2, c("temp_family", "temp_individual"), sep = "_", remove = FALSE) %>% 
  mutate(temp_individual = gsub("\\D", "", temp_individual)) %>% 
  mutate(temp_individual = as.integer(temp_individual)) %>% 
  group_by(X1) %>% 
  slice(which.min(temp_individual)) %>%
  ungroup() %>% 
  select(-starts_with("temp"))

index_sample <- select(index, X2)

write_tsv(index, file.path('.',  index_ped_file), col_names = FALSE, na=".")
write_tsv(index_sample, file.path('.',  index_sample_file), col_names = FALSE, na=".")

