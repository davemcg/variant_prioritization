args <- commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(readxl)
#args <- c("Z:/resources/jaxCNV/genome.jaxCNV.tsv", "Z:/exome/BlueprintGenetics/scramble_anno/all.exome.scramble.v1.xlsx", "Z:/resources/SCRAMBLEvariantClassification.GRCh38.xlsx", "Z:/exome/BlueprintGenetics/scramble_anno/db.test.xlsx")
jaxcnv_file <- args[1] 
output_cohort_file <- args[2]

#"Z:/resources/manta/exome.manta.all.7column.tsv"
all_patient <- read_tsv(jaxcnv_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>% 
  type_convert() %>% 
  mutate(SV_start = round(SV_start + 1, -3), SV_end = round(SV_end, -3)) %>% 
  unite("variant", SV_chrom:SV_type, sep = "-", remove = FALSE)
  
#161960 observations

#finding duplicate
# manta_uniq <- manta %>% 
#   distinct(AnnotSV_ID, Samples_ID, .keep_all = TRUE) %>% 
#   unite("Four_ID", SV_chrom:Samples_ID, sep = "-", remove = FALSE) %>% 
#   group_by(Four_ID) %>% 
#   filter(n()>1)


family_count <- all_patient %>% select(Samples_ID) %>% 
  distinct() %>% 
  nrow()
#394 exome

variant_af_output <- all_patient %>% select(variant, Samples_ID) %>%  
  distinct() %>% 
  group_by(variant) %>% 
  summarise(CohortFreq = n()/family_count, AC = n(), AN = family_count, OGLsamples = paste(Samples_ID, collapse = ", ")) %>% 
  mutate(OGLsamples = ifelse(AC < 6, OGLsamples, "gt 5")) %>%
  unite("NaltP/NtotalP", AC, AN, sep = "/", remove = T) %>% 
  ungroup()

write_tsv(variant_af_output, file.path('.',  output_cohort_file), col_names = TRUE, na=".")
  