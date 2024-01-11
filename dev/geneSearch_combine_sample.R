library(tidyverse)
#library(readxl)

args <- commandArgs(trailingOnly=TRUE)
#setwd("W:/abca4/clinvar.hgmd")
#args <- c("ABCA4.clinvar.hgmd.OGLanno.tsv", "ABCA4.clinvar.hgmd.OGLanno.select.xlsx", "crossmap.hg19.gene.hgmd.clinvar__chr1.tsv", "test.gene.hgmd.clinvar__chr1.ps.tsv")

Input_file <- args[1]
output_file <- args[2]
output_excel <- args[3]
#psOutput_file <- args[4]

OGLanno <- read_tsv(Input_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% mutate(exon = sub("^", " ", exon), intron = sub("^", " ", intron)) 

count <- OGLanno %>% mutate(temp_zygosity = ifelse(gt_types == "1", 1, 2)) %>% group_by(chr_variant_id) %>% summarise(AlleleCount = sum(temp_zygosity))

OGLanno_count <- left_join(OGLanno, count, by = "chr_variant_id")

filtered <- OGLanno_count %>% filter(priority_score > 3 , caller %in% c("dvFb", "fbDvg","fbDv"))

#write_tsv(OGLanno_count, output_file, na = ".")

openxlsx::write.xlsx(list("ps>3" = filtered, "all" = OGLanno_count), file = output_excel, firstRow = TRUE, firstCol = TRUE)

