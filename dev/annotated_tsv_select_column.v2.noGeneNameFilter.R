library(tidyverse)
library(readxl)

args <- commandArgs(trailingOnly=TRUE)
#setwd("W:/abca4/clinvar.hgmd")
#args <- c("ABCA4.clinvar.hgmd.OGLanno.tsv", "ABCA4.clinvar.hgmd.OGLanno.select.xlsx", "crossmap.hg19.gene.hgmd.clinvar__chr1.tsv", "test.gene.hgmd.clinvar__chr1.ps.tsv")

Input_file <- args[1]
output_excel_file <- args[3]
output_tsv_file <- args[2]
#psOutput_file <- args[4]



OGLanno <- read_tsv(Input_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% 
  separate_rows(CSQ, sep = "\\,") %>%
  separate(CSQ, c('allele','consequence','codons','amino_acids','gene','symbol','mane','mane_select','mane_plus_clinical',
                  'feature','exon','intron','hgvsc','hgvsp','max_af','max_af_pops','protein_position','biotype','canonical',
                  'domains','existing_variation','clin_sig','pick','pubmed','phenotypes','sift','polyphen','cadd_raw','cadd_phred',
                  'genesplicer','spliceregion','maxentscan_alt','maxentscan_diff','maxentscan_ref','existing_inframe_oorfs','existing_outofframe_oorfs','existing_uorfs','five_prime_utr_variant_annotation','five_prime_utr_variant_consequence',
                  'motif_name','motif_pos','high_inf_pos','motif_score_change','am_class','am_pathogenicity','refseq_match','bam_edit','source'),
           sep = "\\|", remove = TRUE, convert = TRUE) %>% 
  #gno2x_af_all,gno3_af_all,maxaf_annovar,gno2x_af_popmax,gno3_popmax,gno_gx_ratio,gno2x_an_all,gno3_an_all,gno2x_filter,gno3_filter AND max_af above from VEP.
  mutate(temp_genesplicer = case_when(grepl("High", genesplicer) ~ 3,
                                      grepl("Medium", genesplicer) ~ 1,
                                      TRUE ~ 0)) %>% 
  mutate(temp_maxentscan_diff =  case_when(abs(maxentscan_diff) > 6 ~ 6,
                                           abs(maxentscan_diff) > 3 ~ 3,
                                           TRUE ~ 0)) %>%
  replace_na(list(max_af=0)) %>%
  mutate(CSQ_score = ifelse(grepl("deleterious", sift), 0.5, 0) +
           ifelse(grepl("damaging", polyphen), 0.5, 0) +
           ifelse(is.na(cadd_phred), 0, ifelse(cadd_phred > 15, 0.5, 0) ) +
           temp_genesplicer + temp_maxentscan_diff +
           ifelse(five_prime_utr_variant_consequence == "" | max_af > 0.001, 0, 1) +
           ifelse(am_class == "likely_pathogenic", 0.5, 0) ) %>% 
  group_by(ID) %>% slice(which.max(CSQ_score)) %>% ungroup() %>%
  select(-temp_genesplicer, -temp_maxentscan_diff, -CSQ_score )

openxlsx::write.xlsx(list("OGLanno" = OGLanno), file = output_excel_file, firstRow = TRUE, firstCol = TRUE)
write_tsv(OGLanno, file.path('.', output_tsv_file), na="")

