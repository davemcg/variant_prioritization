args <- commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(readxl)
#args <- c("Z:/exome/BlueprintGenetics/scramble_anno/all.exome.scramble.txt", "Z:/exome/BlueprintGenetics/scramble_anno/all.exome.scramble.v1.xlsx", "Z:/resources/SCRAMBLEvariantClassification.GRCh38.xlsx", "Z:/exome/BlueprintGenetics/scramble_anno/db.test.xlsx")
scramble_del_file <- args[1] 
output_cohort_file <- args[2]

#"Z:/resources/scramble_del/exome.scramble_del.all.7column.tsv"
scramble_del <- read_tsv(scramble_del_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>% 
  type_convert() %>% 
  distinct(AnnotSV_ID, Samples_ID, .keep_all = TRUE) %>% 
  mutate(SV_start = round(SV_start, -3), SV_end = round(SV_end, -3)) %>%
  unite("variant", SV_chrom:SV_type, sep = "-", remove = FALSE)
#161960 observations

#finding duplicate
# scramble_del_uniq <- scramble_del %>% 
#   distinct(AnnotSV_ID, Samples_ID, .keep_all = TRUE) %>% 
#   unite("Four_ID", SV_chrom:Samples_ID, sep = "-", remove = FALSE) %>% 
#   group_by(Four_ID) %>% 
#   filter(n()>1)


family_count <- scramble_del %>% select(Samples_ID) %>% 
  distinct() %>% 
  nrow()
#394 exome

variant_af <- scramble_del %>% 
  select(variant, Samples_ID) %>%  
  group_by(variant) %>% 
  summarise(CohortFreq = n()/family_count, AC = n(), AN = family_count) %>% 
  unite("NaltP/NtotalP", AC, AN, sep = "/", remove = TRUE) %>% 
  ungroup() 
#40329 observations, 37149 < 1% (found in 3 or less); Use 2% as cut-off. Useful strategy?


write_tsv(variant_af, file.path('.',  output_cohort_file), col_names = TRUE, na=".")
  