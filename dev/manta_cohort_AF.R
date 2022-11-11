args <- commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(readxl)
#args <- c("Z:/exome/BlueprintGenetics/scramble_anno/all.exome.scramble.txt", "Z:/exome/BlueprintGenetics/scramble_anno/all.exome.scramble.v1.xlsx", "Z:/resources/SCRAMBLEvariantClassification.GRCh38.xlsx", "Z:/exome/BlueprintGenetics/scramble_anno/db.test.xlsx")
manta_file <- args[1] 
manta_genome_file <- args[2]
output_cohort_file <- args[3]

#"Z:/resources/manta/exome.manta.all.7column.tsv"
manta <- read_tsv(manta_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>% 
  type_convert() %>% 
  distinct(AnnotSV_ID, Samples_ID, .keep_all = TRUE) %>% 
  mutate(SV_start = round(SV_start, -3), SV_end = round(SV_end, -3)) %>% 
  unite("variant", SV_chrom:SV_type, sep = "-", remove = FALSE)
  
#161960 observations

#finding duplicate
# manta_uniq <- manta %>% 
#   distinct(AnnotSV_ID, Samples_ID, .keep_all = TRUE) %>% 
#   unite("Four_ID", SV_chrom:Samples_ID, sep = "-", remove = FALSE) %>% 
#   group_by(Four_ID) %>% 
#   filter(n()>1)


family_count <- manta %>% select(Samples_ID) %>% 
  distinct() %>% 
  nrow()
#394 exome

variant_af <- manta %>% 
  select(variant, Samples_ID) %>%  
  distinct() %>% 
  group_by(variant) %>% 
  summarise(CohortFreq = n()/family_count, AC = n(), AN = family_count) %>% 
  unite("NaltP/NtotalP", AC, AN, sep = "/", remove = TRUE) %>% 
  ungroup() 
#40329 observations, 37149 < 1% (found in 3 or less); Use 2% as cut-off. Useful strategy?

#"Z:/resources/manta/genome.manta.all.7column.tsv"
manta_genome <- read_tsv(manta_genome_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>% 
  type_convert() %>% 
  distinct(AnnotSV_ID, Samples_ID, .keep_all = TRUE) %>% 
  mutate(SV_start = round(SV_start, -3), SV_end = round(SV_end, -3)) %>% 
  unite("variant", SV_chrom:SV_type, sep = "-", remove = FALSE)
#1.4 millions observations

#finding duplicate
# manta_uniq <- manta %>% 
#   distinct(AnnotSV_ID, Samples_ID, .keep_all = TRUE) %>% 
#   unite("Four_ID", SV_chrom:Samples_ID, sep = "-", remove = FALSE) %>% 
#   group_by(Four_ID) %>% 
#   filter(n()>1)


family_count_genome <- manta_genome %>% select(Samples_ID) %>% 
  distinct() %>% 
  nrow()
#135 exome

variant_af_genome <- manta_genome %>% 
  select(variant, Samples_ID) %>%
  distinct() %>% 
  group_by(variant) %>% 
  summarise(CohortFreq = n()/family_count_genome, AC = n(), AN = family_count_genome) %>% 
  unite("NaltP/NtotalP", AC, AN, sep = "/", remove = TRUE) %>% 
  ungroup()

#358435 observations, 312010 in 1 or 2 samples; 3/135=0.022 Use 0.025 as cut-off. Useful strategy?
#MAK exome scramble AF = 1.4%
variant_af_output <- rbind(variant_af, variant_af_genome) %>%
  distinct(variant, .keep_all = TRUE)
#keep the exome AF if common in both   

write_tsv(variant_af_output, file.path('.',  output_cohort_file), col_names = TRUE, na=".")
  