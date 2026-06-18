library(tidyverse)
library(readxl)

args <- commandArgs(trailingOnly=TRUE)
#setwd("W:/abca4/clinvar.hgmd")
#args <- c("ABCA4.clinvar.hgmd.OGLanno.tsv", "ABCA4.clinvar.hgmd.OGLanno.select.xlsx", "crossmap.hg19.gene.hgmd.clinvar__chr1.tsv", "test.gene.hgmd.clinvar__chr1.ps.tsv")

Input_file <- args[1]
geneNames <- args[2]
output_file <- args[3]
#psOutput_file <- args[4]

print(geneNames)

OGLanno <- read_tsv(Input_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% 
  separate_rows(CSQ, sep = "\\,") %>%
  separate(CSQ, c('allele','consequence','codons','amino_acids','gene','symbol','mane','mane_select','mane_plus_clinical',
                  'feature','exon','intron','hgvsc','hgvsp','max_af','max_af_pops','protein_position','biotype','canonical',
                  'domains','existing_variation','clin_sig','pick','pubmed','phenotypes','sift','polyphen','cadd_raw','cadd_phred',
                  'genesplicer','spliceregion','maxentscan_alt','maxentscan_diff','maxentscan_ref','existing_inframe_oorfs','existing_outofframe_oorfs','existing_uorfs','five_prime_utr_variant_annotation','five_prime_utr_variant_consequence',
                  'motif_name','motif_pos','high_inf_pos','motif_score_change','am_class','am_pathogenicity','refseq_match','bam_edit','source'),
           sep = "\\|", remove = TRUE, convert = TRUE) %>% 
  filter(canonical == "YES", grepl(geneNames, symbol))
openxlsx::write.xlsx(list("OGLanno" = OGLanno), file = output_file, firstRow = TRUE, firstCol = TRUE)

