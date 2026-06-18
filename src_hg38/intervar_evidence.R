#This script take the OGL_variant excel file and produce a InterVar Evidence file. 

args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
#args <- c("Y:/resources/OGLvariantsClassification.xlsx",
#         "Y:/resources/evidence.input", "sorted.tsv", "filtered.tsv", "G05068", "filtered.xlsx")

library(tidyverse)
library(readxl)

evidence_vcf <-  read_xlsx(args[1], sheet = "OGLvariants", na =c("", ".", "NA", " ")) %>% 
  filter(!is.na(hg38_variant_id)) %>% 
  separate(hg38_variant_id, c("#CHROM", "POS", "REF", "ALT"), remove = TRUE, convert = TRUE) %>% 
  mutate(QUAL = NA, FILTER = NA) %>% # these are just names, Convert2annova does not care the contents or names
  select("#CHROM", "POS", "InterVar_evidence", "REF", "ALT", "QUAL", "FILTER", "annovarAnnotation")
write_tsv(evidence_vcf, file = args[2])
