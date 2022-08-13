args <- commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(readxl)

can31 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/can3.01.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
  rename_all(funs(str_replace(., "can3.01", "")))
can32 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/can3.02.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
  rename_all(funs(str_replace(., "can3.02", "")))
can33 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/can3.03.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
  rename_all(funs(str_replace(., "can3.03", "")))
can314 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/can3.14.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
  rename_all(funs(str_replace(., "can3.14", "")))
can316 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/can3.16.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
  rename_all(funs(str_replace(., "can3.16", "")))

all <- bind_rows(can31, can32, can33, can314, can316) %>% 
  group_by(chr_variant_id) %>% summarize(NoAffected = n())

can31f <- left_join(can31, all, by = "chr_variant_id") %>% filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
can32f <- left_join(can32, all, by = "chr_variant_id") %>% filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
can33f <- left_join(can33, all, by = "chr_variant_id") %>% filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
can314f <- left_join(can314, all, by = "chr_variant_id") %>% filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
can316f <- left_join(can316, all, by = "chr_variant_id") %>% filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)

openxlsx::write.xlsx(list("can3.01" = can31f, "can3.02" = can32f, "can3.03" = can33f, "can3.14" = can314f, "can3.16" = can316f), 
                     file = "Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/can.R.filtered.xlsx", firstRow = TRUE, firstCol = TRUE)


w55_12 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.1.2.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
  rename_all(funs(str_replace(., "W55.1.2", "")))
w55_13 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.1.3.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
  rename_all(funs(str_replace(., "W55.1.3", "")))
w55_21 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.2.1.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
  rename_all(funs(str_replace(., "W55.2.1", "")))
w55_22 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.2.2.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
  rename_all(funs(str_replace(., "W55.2.2", "")))
w55_23 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.2.3.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
  rename_all(funs(str_replace(., "W55.2.3", "")))
w55_24 <- read_xlsx("Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.2.4.FIHPt.GRCh37.fih.gemini.filtered.xlsx", sheet = "all", na = c("NA", "", "None", ".")) %>% 
  rename_all(funs(str_replace(., "W55.2.4", "")))

w55_all <- bind_rows(w55_12, w55_21, w55_22, w55_23, w55_24) %>% 
  group_by(chr_variant_id) %>% summarize(NoAffected = n())

w55_12f <- anti_join(w55_12, w55_13, by = "chr_variant_id") %>% left_join(., w55_all, by = "chr_variant_id") %>% 
  filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
w55_21f <- anti_join(w55_21, w55_13, by = "chr_variant_id") %>% left_join(., w55_all, by = "chr_variant_id") %>% 
  filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
w55_22f <- anti_join(w55_22, w55_13, by = "chr_variant_id") %>% left_join(., w55_all, by = "chr_variant_id") %>% 
  filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
w55_23f <- anti_join(w55_23, w55_13, by = "chr_variant_id") %>% left_join(., w55_all, by = "chr_variant_id") %>% 
  filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)
w55_24f <- anti_join(w55_24, w55_13, by = "chr_variant_id") %>% left_join(., w55_all, by = "chr_variant_id") %>% 
  filter(NoAffected > 3, pmaxaf < 0.001, aaf < 0.5, priority_score > 3 ) %>% arrange(chrom, start_vcf)


openxlsx::write.xlsx(list("W55.1.2" = w55_12f, "W55.2.1" = w55_21f, "W55.2.2" = w55_22f, "W55.2.3" = w55_23f, "W55.2.4" = w55_24f), 
                     file = "Z:/NextSeqAnalysis/FIH/prioritization/gemini_xlsx/W55.R.filtered.xlsx", firstRow = TRUE, firstCol = TRUE)
