args <- commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(readxl)

#Rscript family_multiple_affected_filter_Gluacoma.R gemini_tsv_filtered/G04556.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04742.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04745.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04746.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04748.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04839.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04842.lp22-09.nano.gemini.filtered.tsv filePlaceholder gemini_tsv_filtered/G04743.lp22-09.nano.gemini.filtered.tsv gemini_tsv_filtered/G04747.lp22-09.nano.gemini.filtered.tsv
affected1_file <- args[1]
affected2_file <- args[2]
affected3_file <- args[3]
affected4_file <- args[4]
affected5_file <- args[5]
affected6_file <- args[6]
affected7_file <- args[7]
affected8_file <- args[8]
unaffected1_file <- args[9]
unaffected2_file <- args[10]


affected1 <- read_tsv(affected1_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% rename_all(funs(str_replace(., "\\.G0\\d{4}", "")))
affected2 <- read_tsv(affected2_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% rename_all(funs(str_replace(., "\\.G0\\d{4}", "")))
affected3 <- read_tsv(affected3_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% rename_all(funs(str_replace(., "\\.G0\\d{4}", "")))
affected4 <- read_tsv(affected4_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% rename_all(funs(str_replace(., "\\.G0\\d{4}", "")))
affected5 <- read_tsv(affected5_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% rename_all(funs(str_replace(., "\\.G0\\d{4}", "")))
affected6 <- read_tsv(affected6_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% rename_all(funs(str_replace(., "\\.G0\\d{4}", "")))
#affected7 <- read_tsv(affected7_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
#  type_convert() %>% rename_all(funs(str_replace(., "\\.G0\\d{4}", "")))

# affected8 <- read_tsv(affected8_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
#   type_convert() %>% rename_all(funs(str_replace(., "\\.G0\\d{4}", "")))


unaffected1 <- read_tsv(unaffected1_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% rename_all(funs(str_replace(., "\\.G0\\d{4}", "")))
unaffected2 <- read_tsv(unaffected2_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% rename_all(funs(str_replace(., "\\.G0\\d{4}", "")))

un_affected <- bind_rows(unaffected1, unaffected2) %>% select(chr_variant_id) %>% distinct()

all_affected <- bind_rows(affected1, affected2, affected3, affected4, affected5, affected6) %>% 
  anti_join(., un_affected, by = "chr_variant_id")

affected_count <- all_affected %>% group_by(chr_variant_id) %>% summarize(NoAffected = n())
samples_affected <- all_affected %>% group_by(chr_variant_id) %>%
  mutate(Aff_samples = paste(sample, collapse = ",")) %>% 
  select(chr_variant_id, Aff_samples) %>% 
  distinct() %>% 
  ungroup() %>% 
  left_join(., affected_count, by = "chr_variant_id" )


present_6_affected_only <- anti_join(affected1, un_affected, by = "chr_variant_id") %>%
  left_join(., affected_count, by = "chr_variant_id") %>% filter(NoAffected == 6, is.na(pmaxaf) | pmaxaf < 0.025, is.na(af_oglg) | af_oglg < 0.05, priority_score > 0 ) %>% arrange(chrom, start_vcf)
openxlsx::write.xlsx(list("6affected" = present_6_affected_only), 
                     file = "F01926/F01926.in.6.myopia.affected.xlsx", firstRow = TRUE, firstCol = TRUE)
filter_affected <- function(dataframe, filename){
  affected_f <- anti_join(dataframe, un_affected, by = "chr_variant_id") %>%
    left_join(., samples_affected, by = "chr_variant_id") %>%
    filter(NoAffected > 4, is.na(pmaxaf) | pmaxaf < 0.025, is.na(af_oglg) | af_oglg < 0.05, priority_score > 0 ) %>% arrange(chrom, start_vcf) # at least in 4 of 5
  openxlsx::write.xlsx(list("In>4affected" = affected_f), 
                       file = filename, firstRow = TRUE, firstCol = TRUE)
}

filter_affected(affected1, "F01926/G04742.In5orMore.myopia.affected.xlsx")
filter_affected(affected2, "F01926/G04743.In5orMore.myopia.affected.xlsx")
filter_affected(affected3, "F01926/G04745.In5orMore.myopia.affected.xlsx")
filter_affected(affected4, "F01926/G04746.In5orMore.myopia.affected.xlsx")
filter_affected(affected5, "F01926/G04747.In5orMore.myopia.affected.xlsx")
filter_affected(affected6, "F01926/G04748.In5orMore.myopia.affected.xlsx")
#filter_affected(affected7, "F01926/G04842.In5orMore.affected.xlsx")


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
