args <- commandArgs(trailingOnly=TRUE)

library(tidyverse)

# args <- c("Z:/genome/RodYoung/prioritization/chileanMAC.ped",
#           "Z:/genome/RodYoung/clinSV/RY1.clinSV.RARE_PASS_GENE.annotated.tsv",
#           "Z:/development/genome/clinSV/D1596_1.RARE_PASS_GENE.eG.tsv",
#           "Z:/development/genome/clinSV/D1596_1.RARE_PASS_GENE.eG.filtered.xlsx")

all_sample_ped_file <- args[1]
index_ped_file <- args[2]
index_sample_file <- args[3]

#"Z:/resources/OGLsample/genome.2022-08.ped"

all_sample <- read_tsv(all_sample_ped_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert() 

index <- filter(all_sample, DiseaseStatus == "2") %>% 
  separate(SubjectID, c("temp_family", "temp_individual"), sep = "_|x", remove = FALSE) %>% 
  mutate(temp_individual = gsub("\\D", "", temp_individual)) %>% 
  mutate(temp_individual = as.integer(temp_individual)) %>%
  mutate(temp_individual = ifelse(is.na(temp_individual), "1", temp_individual)) %>% 
  group_by(`##FamilyID`) %>% 
  slice(which.min(temp_individual)) %>%
  ungroup() %>% 
  select(-starts_with("temp"))

index_sample <- select(index, SubjectID)

write_tsv(index, file.path('.',  index_ped_file), col_names = FALSE, na=".")
write_tsv(index_sample, file.path('.',  index_sample_file), col_names = FALSE, na=".")

