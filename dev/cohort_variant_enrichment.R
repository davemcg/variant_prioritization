args <- commandArgs(trailingOnly=TRUE)
# args <- c( "Z:/NextSeqAnalysis/221118_221121/prioritization/gemini_tsv_filtered/variant.enrich.input.txt", "Z://NextSeqAnalysis/221118_221121/prioritization/gemini_tsv_filtered/enriched.variant.xlsx", "Z:/genome/nisc_chuuk/bcmlocus/CORDX.tsv",
#            "Z:/resources/OGLpanelGeneDxORcandidate.xlsx", "rearranged.tsv", "filtered.tsv", "108976P", "filtered.xlsx", "0.5", "W:/ddl_nisc_custom_capture/042020/CoNVaDING/CNV_hiSens/108976P.b37.aligned.only.best.score.shortlist.txt")

input_file <- args[1]
output_file <- args[2]
# bcm_coverage_bed_file <- args[3]
# bcmlocusCN_file <- args[4]
# bcmlocusCN_wide_file <- args[5]

library(tidyverse)
library(readxl)
#library(vroom)

row_fisher <- function(row, alt = 'two.sided', cnf = 0.95) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(OR = round(f$estimate[[1]], 2),
           Conf_Int = paste(round(f$conf.int, 2), collapse = ";"),
           P_value = signif(f$p.value, 3)))
}
#use: p <- data.frame(t(apply(test_df, 1, row_fisher)))

input <- read_tsv(input_file, col_names = TRUE, na = c("NA", "", "None"), col_types = cols(.default = col_character())) %>% 
  select(ref_gene:hg38_pos, gt_types, aaf:cherry_sum_score,gno2g_ac_all:gno3_popmax,ac:an, chrom, gno3_nonpar) %>% 
  type_convert() %>% 
  unite(Sample_zygosity, sample, gt_types, sep = "-", remove = TRUE) %>% 
  rename(sample = Sample_zygosity) %>% 
  select(-af, ref_gene, sample, everything()) %>% 
  filter(priority_score > 2) %>% 
  unique() %>% 
  group_by(chr_variant_id) %>% 
  mutate(sample = paste(sample, collapse = ",")) %>% 
  unique() %>% 
  ungroup() %>% 
  mutate(gno3_an_all = case_when(is.na(gno3_an_all) & chrom %in% c("X", "chrX") & gno3_nonpar == "1" ~ 116830,
                                 is.na(gno3_an_all) & chrom %in% c("Y", "chrY") & gno3_nonpar == "1" ~ 35482,
                                 is.na(gno3_an_all) ~ 152312,
                                 TRUE ~ gno3_an_all)) %>% 
  replace_na(list(gno3_ac_all=0)) %>% 
  mutate(temp_ref_cohort = an - ac, temp_ref_gnomad = gno3_an_all - gno3_ac_all) %>% 
  filter(!ref_gene %in% c("OPN1LW", "OPN1MW"))

variant_ac <- input %>% select(ac, temp_ref_cohort, gno3_ac_all, temp_ref_gnomad)

df_fisher <- data.frame(t(apply(variant_ac, 1, row_fisher)))

variant_fisher <- cbind(input,df_fisher) %>% 
  mutate(OR = as.character(OR)) %>% 
  mutate(OR = ifelse(OR == "Inf","1000000000",OR)) %>% 
  mutate(P_value = as.numeric(as.character(P_value)),
         OR = as.numeric(as.character(OR))) %>% 
  select(-starts_with("temp_")) %>% 
  select(ref_gene:aaf, ac:an, caller:af_oglg,OR:P_value,everything()) %>% 
  filter(P_value < 0.00001, OR > 1, caller %in% c("fbDvg", "fbDv"))

openxlsx::write.xlsx(list("enriched" = variant_fisher), 
                     file = output_file, firstRow = TRUE, firstCol = TRUE)
