args <- commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(readxl)

#Rscript family_multiple_affected_filter_Gluacoma.R gemini_tsv_filtered/G04556.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04742.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04745.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04746.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04748.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04839.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04842.lp22-09.nano.gemini.filtered.tsv filePlaceholder gemini_tsv_filtered/G04743.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04747.lp22-09.nano.gemini.filtered.tsv
# args <- c("G00582.prasov23-1.23sample.gemini.filtered.tsv", 
#           "G00594.prasov23-1.23sample.gemini.filtered.tsv",
#           "G00595.prasov23-1.23sample.gemini.filtered.tsv", 
#           "G00614.prasov23-1.23sample.gemini.filtered.tsv",
#           "G00628.prasov23-1.23sample.gemini.filtered.tsv", 
#           "G00572.prasov23-1.23sample.gemini.filtered.tsv",
#           "G00598.prasov23-1.23sample.gemini.filtered.tsv", 
#           "G00602.prasov23-1.23sample.gemini.filtered.tsv",
#           "G00609.prasov23-1.23sample.gemini.filtered.tsv", 
#           "G00626.prasov23-1.23sample.gemini.filtered.tsv")

affected1_file <- args[1]
affected2_file <- args[2]
affected3_file <- args[3]
affected4_file <- args[4]
affected5_file <- args[5]
unaffected1_file <- args[6]
unaffected2_file <- args[7]
unaffected3_file <- args[8]
unaffected4_file <- args[9]
unaffected5_file <- args[10]

#F00058	G00572				Unaffected-pool
#F00058	G00582				Affected-pool
#F00058	G00594				Affected-pool
#F00058	G00595				Affected-pool
#F00058	G00598				Unaffected
#F00058	G00602				Unaffected-pool
#F00058	G00609				Unaffected
#F00058	G00614				Affected-pool
#F00058	G00626				Unaffected
#F00058	G00628				Affected

affected1 <- read_tsv(affected1_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  rename_all(list(~str_replace(., "\\.G0\\d{4}", "")))
affected2 <- read_tsv(affected2_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  rename_all(list(~str_replace(., "\\.G0\\d{4}", "")))
affected3 <- read_tsv(affected3_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  rename_all(list(~str_replace(., "\\.G0\\d{4}", "")))
affected4 <- read_tsv(affected4_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  rename_all(list(~str_replace(., "\\.G0\\d{4}", "")))
affected5 <- read_tsv(affected5_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
 rename_all(list(~str_replace(., "\\.G0\\d{4}", "")))


unaffected1 <- read_tsv(unaffected1_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  rename_all(list(~str_replace(., "\\.G0\\d{4}", "")))
unaffected2 <- read_tsv(unaffected2_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  rename_all(list(~str_replace(., "\\.G0\\d{4}", "")))
unaffected3 <- read_tsv(unaffected3_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  rename_all(list(~str_replace(., "\\.G0\\d{4}", "")))
unaffected4 <- read_tsv(unaffected4_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  rename_all(list(~str_replace(., "\\.G0\\d{4}", "")))
unaffected5 <- read_tsv(unaffected5_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  rename_all(list(~str_replace(., "\\.G0\\d{4}", "")))

un_affected <- bind_rows(unaffected1, unaffected2, unaffected3, unaffected4, unaffected5) %>% 
  select(chr_variant_id, sample) %>% group_by(chr_variant_id) %>% 
  summarise(samples_unaffected=paste(sample, collapse = ","), number_unaffected = n()) %>%
  ungroup()

variant_affected_gt3 <- bind_rows(affected1, affected2, affected3, affected4, affected5) %>% 
  select(chr_variant_id, sample) %>% group_by(chr_variant_id) %>% 
  summarise(samples_affected=paste(sample, collapse = ","), number_affected = n()) %>% 
  ungroup() %>% 
  filter(number_affected > 3) %>% 
  left_join(., un_affected, by = "chr_variant_id") %>% 
  replace_na(list(number_unaffected=0)) %>% 
  filter(number_unaffected < 2)

variant_selected <- variant_affected_gt3 %>% pull(chr_variant_id) 

all_column_affected_gt3 <- bind_rows(affected1, affected2, affected3, affected4, affected5) %>%
 select(ref_gene, chr_variant_id:grch37variant_id, panel_class:spliceaimasked50max,gno2x_filter:syn_z, ensembl_gene_id:an, clinpred_pred:primateai_score,revel_rankscore:motif_score_change,chrom,start_vcf) %>%
  filter(chr_variant_id %in% variant_selected) %>% distinct() %>% 
  left_join(., variant_affected_gt3, by = "chr_variant_id" ) %>% 
  type_convert() %>% 
  select(ref_gene:grch37variant_id,samples_affected:number_unaffected,everything()) %>% 
  filter(is.na(pmaxaf) | pmaxaf < 0.025, is.na(af_oglg) | af_oglg < 0.05, priority_score > 0 ) %>% 
  arrange(chrom, start_vcf)

openxlsx::write.xlsx(list("allow1mis" = all_column_affected_gt3), 
                     file = "family/F000058.allow1mis.all.variant.filtered.xlsx", firstRow = TRUE, firstCol = TRUE)
filter_affected <- function(dataframe, filename){
  affected_f <- filter(dataframe, chr_variant_id %in% variant_selected) %>% 
    left_join(., variant_affected_gt3, by = "chr_variant_id" ) %>% 
    filter(is.na(pmaxaf) | pmaxaf < 0.025, is.na(af_oglg) | af_oglg < 0.05, priority_score > 0 ) %>% 
    select(ref_gene:grch37variant_id,samples_affected:number_unaffected,everything()) %>% 
    arrange(chrom, start_vcf)
  openxlsx::write.xlsx(list("allow1mis" = affected_f), 
                       file = filename, firstRow = TRUE, firstCol = TRUE)
}

filter_affected(affected1, "family/G00582.allow1mis.filtered.xlsx")
filter_affected(affected2, "family/G00594.allow1mis.filtered.xlsx")
filter_affected(affected3, "family/G00595.allow1mis.filtered.xlsx")
filter_affected(affected4, "family/G00614.allow1mis.filtered.xlsx")
filter_affected(affected5, "family/G00628.allow1mis.filtered.xlsx")

# affected1f <- anti_join(affected1, un_affected, by = "chr_variant_id") %>%
#   left_join(., all_affected, by = "chr_variant_id") %>% filter(NoAffected > 5, is.na(pmaxaf) | pmaxaf < 0.01, is.na(af_oglg) | af_oglg < 0.05, priority_score > 0 ) %>% arrange(chrom, start_vcf) # at least in 4 of 5
# affected2f <- anti_join(affected2, un_affected, by = "chr_variant_id") %>%
#   left_join(., all_affected, by = "chr_variant_id") %>% filter(NoAffected > 2, is.na(pmaxaf) | pmaxaf < 0.01, is.na(af_oglg) | af_oglg < 0.05, priority_score > 0 ) %>% arrange(chrom, start_vcf)
# affected3f <- anti_join(affected3, un_affected, by = "chr_variant_id") %>%
#   left_join(., all_affected, by = "chr_variant_id") %>% filter(NoAffected > 2, is.na(pmaxaf) | pmaxaf < 0.01, is.na(af_oglg) | af_oglg < 0.05, priority_score > 0 ) %>% arrange(chrom, start_vcf)
# affected4f <- anti_join(affected4, un_affected, by = "chr_variant_id") %>%
#   left_join(., all_affected, by = "chr_variant_id") %>% filter(NoAffected > 2, is.na(pmaxaf) | pmaxaf < 0.01, is.na(af_oglg) | af_oglg < 0.05, priority_score > 0 ) %>% arrange(chrom, start_vcf)
#affected5f <- anti_join(affected5, un_affected, by = "chr_variant_id") %>%
#  left_join(., all_affected, by = "chr_variant_id") %>% filter(NoAffected > 2, is.na(pmaxaf) | pmaxaf < 0.01, is.na(af_oglg) | af_oglg < 0.05, priority_score > 0 ) %>% arrange(chrom, start_vcf)


# openxlsx::write.xlsx(list("G04742" = affected1f), 
#                      file = "G04742.4.affected.absent.in.unaffected.xlsx", firstRow = TRUE, firstCol = TRUE)
#openxlsx::write.xlsx(list("G04743" = affected2f), 
#                     file = "G04743.variant.absent.in.unaffected.xlsx", firstRow = TRUE, firstCol = TRUE)
# openxlsx::write.xlsx(list("G04746" = affected2f), 
#                      file = "G04746.4.affected.absent.in.unaffected.xlsx", firstRow = TRUE, firstCol = TRUE)
# openxlsx::write.xlsx(list("G04839" = affected3f), 
#                      file = "G04839.4.affected.absent.in.unaffected.xlsx", firstRow = TRUE, firstCol = TRUE)
# openxlsx::write.xlsx(list("G04842" = affected4f), 
#                      file = "G04842.4.affected.absent.in.unaffected.xlsx", firstRow = TRUE, firstCol = TRUE)

# w55_12 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.1.2.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
#   rename_all(funs(str_replace(., "W55.1.2", "")))
# w55_13 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.1.3.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
#   rename_all(funs(str_replace(., "W55.1.3", "")))
# w55_21 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.2.1.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
#   rename_all(funs(str_replace(., "W55.2.1", "")))
# w55_22 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.2.2.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
#   rename_all(funs(str_replace(., "W55.2.2", "")))
# w55_23 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.2.3.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
#   rename_all(funs(str_replace(., "W55.2.3", "")))
# w55_24 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.2.4.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
#   rename_all(funs(str_replace(., "W55.2.4", "")))
# 
# w55_all <- bind_rows(w55_12, w55_21, w55_22, w55_23, w55_24) %>% 
#   group_by(chr_variant_id) %>% summarize(NoAffected = n())
# 
# w55_12f <- anti_join(w55_12, w55_13, by = "chr_variant_id") %>% left_join(., w55_all, by = "chr_variant_id") %>% 
#   filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
# w55_21f <- anti_join(w55_21, w55_13, by = "chr_variant_id") %>% left_join(., w55_all, by = "chr_variant_id") %>% 
#   filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
# w55_22f <- anti_join(w55_22, w55_13, by = "chr_variant_id") %>% left_join(., w55_all, by = "chr_variant_id") %>% 
#   filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
# w55_23f <- anti_join(w55_23, w55_13, by = "chr_variant_id") %>% left_join(., w55_all, by = "chr_variant_id") %>% 
#   filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
# w55_24f <- anti_join(w55_24, w55_13, by = "chr_variant_id") %>% left_join(., w55_all, by = "chr_variant_id") %>% 
#   filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
# 
# 
# openxlsx::write.xlsx(list("W55.1.2" = w55_12f, "W55.2.1" = w55_21f, "W55.2.2" = w55_22f, "W55.2.3" = w55_23f, "W55.2.4" = w55_24f), 
#                      file = "Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.R.filtered.xlsx", firstRow = TRUE, firstCol = TRUE)
