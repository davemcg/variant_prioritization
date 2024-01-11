library(tidyverse)
#library(readxl)

args <- commandArgs(trailingOnly=TRUE)
#setwd("W:/abca4/clinvar.hgmd")
#args <- c("ABCA4.clinvar.hgmd.OGLanno.tsv", "ABCA4.clinvar.hgmd.OGLanno.select.xlsx", "crossmap.hg19.gene.hgmd.clinvar__chr1.tsv", "test.gene.hgmd.clinvar__chr1.ps.tsv")

Input_file <- args[1]
geneNames <- args[2]
sampleName <- args[3]
output_file <- args[4]
#psOutput_file <- args[4]

OGLanno <- read_tsv(Input_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% 
  rename_all(funs(str_replace(., sampleName, ""))) %>% rename_all(funs(str_replace(., "\\.$", ""))) %>% 
  select('gene', 'sample', 'chr_variant_id', 'grch37variant_id', 'gts', 'gt_types',
         'gt_depths', 'gt_alt_freqs', 'gt_quals', 'aaf', 'caller', 'panel_class',
         'priority_score', 'prscore_intervar', 'clinvar_hgmd_score', 'splice_score',
          'gno2e3g_af', 'gno2e3g_acan', 'pmaxaf', 'gno2x_af_all',
         'gno2x_filter', 'gno3_af_all', 'gno3_filter', 'max_af', 'max_af_pops',
         'gno2e3g_hom', 'ref_gene', 'func_refgenewithver', 'exonicfunc_refgenewithver',
         'refgenewithver', 'mane_select', 'hgvsc', 'hgvsp', 'exon', 'intron', 'aa_length',
         'omim_gene', 'omim_inheritance', 'omim_phen', 'pvs1', 'truncating_vep', 'hgmd_id',
         'hgmd_class', 'hgmd_phen', 'hgmd_overlap4aa', 'existing_variation', 'clnalleleid',
         'clnsig', 'clin_sig', 'clnrevstat', 'clndn', 'clndisdb') %>% 
 filter(gene == geneNames | ref_gene == geneNames)

write_tsv(OGLanno, output_file, na = ".")

#openxlsx::write.xlsx(list("OGLanno" = OGLanno), file = output_file, firstRow = TRUE, firstCol = TRUE)

