args <- commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(readxl)

# args <- c("C:/Users/guanb/OneDrive - National Institutes of Health/Projects/eG-AA/AA_Batch_Validator.xlsx", 
#           "C:/Users/guanb/OneDrive - National Institutes of Health/Projects/eG-AA/aa.vcf.tsv",
#           "C:/Users/guanb/OneDrive - National Institutes of Health/Projects/EYS/input/4.2.1.2 Participant Variant Data.txt",
#           "C:/Users/guanb/OneDrive - National Institutes of Health/Projects/EYS/input/4.2.1.1 Participant Data.txt")

variant_validator_file <- args[1]
vcf_file <- args[2]

variant_validator <- read_xlsx(variant_validator_file, sheet = "Batch Validator", na = c("NA", "", "None", "NONE", ".")) %>%
  select(Input,GRCh38_CHR:GRCh38_ALT) %>% 
  mutate(GRCh38_CHR = sub("^", "chr", GRCh38_CHR)) %>% 
  select( GRCh38_CHR, GRCh38_POS, Input, GRCh38_REF, GRCh38_ALT) %>% 
  rename(`#CHROM` = GRCh38_CHR, "POS" = GRCh38_POS, "ID" = Input, "REF" = GRCh38_REF, "ALT" = GRCh38_ALT) %>% 
  distinct() %>%
  filter(!is.na(`#CHROM`)) %>% 
  mutate(`#CHROM` = factor(`#CHROM`, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                                "chr20","chr21","chr22","chrX","chrY"))) %>% 
  arrange(`#CHROM`, POS)

write_tsv(variant_validator, vcf_file, na=".")