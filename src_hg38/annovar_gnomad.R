args <- commandArgs(trailingOnly=TRUE)
#args <- c("temp/Prasov22-1__chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY.avinput.hg38_multianno.txt.intervar", "temp/Prasov22-1__chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY.avinput.hg38_multianno.modified.txt", "temp/spliceai.Prasov22-1__chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY.tsv.cut", "Z:/resources/HGMD/HGMDtranscript.txt", "temp/20200422.freebayes__1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y.spliceai_annovar_intervar")

library(tidyverse)
#library(vroom)

annovar_file <- args[1]
gnomad_file <- args[2]
output_file <- args[3]

annovar <- read_tsv(annovar_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>% type_convert() %>%
  mutate(Chr = as.factor(Chr), CHROM = as.factor(CHROM))

gnomad <- read_tsv(gnomad_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>% type_convert() %>%
  mutate(CHROM = as.factor(CHROM))

annovar_gnomad <- left_join(annovar, gnomad, by = c("CHROM", "POS", "REF", "ALT")) 
#  select(Chr:Alt, gno3_af_all:gno3_af_oth, everything()) 

write_tsv(annovar_gnomad, file.path('.',  output_file), na=".")
