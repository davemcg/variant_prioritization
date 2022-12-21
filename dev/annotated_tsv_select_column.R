library(tidyverse)
library(readxl)

args <- commandArgs(trailingOnly=TRUE)
#setwd("W:/abca4/clinvar.hgmd")
#args <- c("ABCA4.clinvar.hgmd.OGLanno.tsv", "ABCA4.clinvar.hgmd.OGLanno.select.xlsx", "crossmap.hg19.gene.hgmd.clinvar__chr1.tsv", "test.gene.hgmd.clinvar__chr1.ps.tsv")

Input_file <- args[1]
output_file <- args[2]
crossmap_file <- args[3]
psOutput_file <- args[4]

OGLanno <- read_tsv(Input_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% 
  separate_rows(CSQ, sep = "\\,") %>%
  separate(CSQ, c("Allele","Consequence","Codons","Amino_acids","Gene","SYMBOL","MANE_SELECT","Feature","EXON","INTRON","HGVSc","HGVSp","MAX_AF","MAX_AF_POPS","Protein_position","BIOTYPE","CANONICAL","DOMAINS","Existing_variation","CLIN_SIG","PICK","PUBMED","Phenotypes","SIFT","PolyPhen","CADD_RAW","CADD_PHRED","GeneSplicer","SpliceRegion","MaxEntScan_alt","MaxEntScan_diff","MaxEntScan_ref","existing_InFrame_oORFs","existing_OutOfFrame_oORFs","existing_uORFs","five_prime_UTR_variant_annotation","five_prime_UTR_variant_consequence","Mastermind_counts","Mastermind_MMID3","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE"), sep = "\\|", remove = TRUE) %>% 
  filter(CANONICAL == "YES", SYMBOL == "ABCA4") %>% 
  mutate(HGVSc = sub("ENST\\d*\\.\\d", "NM_000350.3", HGVSc), 
         HGVSp = sub("ENSP\\d*\\.\\d:", "", HGVSp),
         Mastermind_MMID3 =ifelse(is.na(Mastermind_MMID3), Mastermind_MMID3, sub("^", "https://mastermind.genomenon.com/detail?mutation=", Mastermind_MMID3)) ) %>% 
  select(CHROM:Consequence,SYMBOL,EXON:MAX_AF_POPS,PUBMED, DOMAINS,Existing_variation,SIFT:MaxEntScan_ref,Mastermind_counts,Mastermind_MMID3,
         priority_score:spliceai_maxscore,grch37variant_id:InterVar_and_Evidence,gno2x_ac_all:hgmd_phen,ClinPred_score:hmc_score)
#mastermind better to use HGVSG & MMID3 as the ID. Next version.
openxlsx::write.xlsx(list("OGLanno" = OGLanno), file = output_file, firstRow = TRUE, firstCol = TRUE)

