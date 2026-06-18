library(tidyverse)
#library(readxl)

args <- commandArgs(trailingOnly=TRUE)
#setwd("W:/abca4/clinvar.hgmd")
#args <- c("Z:/geneSearch/exome.tsv", "586", "crossmap.hg19.gene.hgmd.clinvar__chr1.tsv", "test.gene.hgmd.clinvar__chr1.ps.tsv")

Input_file <- args[1]
probandNo <- args[2]
output_excel <- args[3]
#psOutput_file <- args[4]

row_fisher <- function(row, alt = 'two.sided', cnf = 0.95) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(OR = round(f$estimate[[1]], 2),
           Conf_Int = paste(round(f$conf.int, 2), collapse = ";"),
           P_value = signif(f$p.value, 3)))
}
#use: p <- data.frame(t(apply(test_df, 1, row_fisher)))
#cohort_an and ac only valid for autosome chromosomes

OGLanno <- read_tsv(Input_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  type_convert() %>% mutate(exon = sub("^", " ", exon), intron = sub("^", " ", intron)) %>% 
  mutate(cohort_an = 2*as.integer(probandNo)) %>% 
  replace_na(list(gno3_ac_all = 0, gno3_an_all = 152312)) %>% 
  mutate(temp_zygosity = ifelse(gt_types == "1", 1, 2))

count <- OGLanno %>%
  select(sample, chr_variant_id, temp_zygosity) %>% 
  separate(sample, c("familyID", "otherField"), sep = "_|x", remove = TRUE) %>% 
  select(-otherField) %>% 
  group_by(familyID, chr_variant_id) %>% 
  slice(which.max(temp_zygosity)) %>% ungroup() %>% 
  group_by(chr_variant_id) %>% summarise(AlleleCount = sum(temp_zygosity))

OGLanno_count <- left_join(OGLanno, count, by = "chr_variant_id")

variant_ac <- OGLanno_count %>% select(AlleleCount, cohort_an, gno3_ac_all, gno3_an_all) %>% 
  mutate(cohort_ref = cohort_an - AlleleCount,
         gno3_ref = gno3_an_all - gno3_ac_all) %>% 
  select(AlleleCount, cohort_ref, gno3_ac_all, gno3_ref)

df_fisher <- data.frame(t(apply(variant_ac, 1, row_fisher)))

OGLanno_count <- cbind(OGLanno_count, df_fisher) %>% 
  mutate(OR = as.character(OR)) %>% 
  mutate(OR = ifelse(OR == "Inf","1000000000",OR)) %>% 
  mutate(P_value = as.numeric(as.character(P_value)),
         OR = as.numeric(as.character(OR))) 

filtered <- OGLanno_count %>% filter(priority_score >=5 , 
                                     caller %in% c("dvFb", "fbDvg","fbDv", "clr3g"))  
recessive_count <- select(filtered, c(sample, temp_zygosity, AlleleCount)) %>% 
  group_by(sample) %>% summarize(recessive_cnt = sum(temp_zygosity), maxAC = max(AlleleCount)) 

filtered <- left_join(filtered, recessive_count, by = "sample") %>%
  select(sample, recessive_cnt, AlleleCount, everything()) %>% 
  arrange(desc(recessive_cnt),desc(maxAC), sample, gene, desc(AlleleCount), desc(priority_score), desc(prscore_intervar), chr_variant_id) %>% 
  select(-maxAC)
#write_tsv(OGLanno_count, output_file, na = ".")

readme <- data.frame("Item" = "Fisher", "Note" = "Allele No not accurate for chrX/Y/M") %>% 
  add_row("Item" = "Proband No", "Note" = paste0(probandNo, ", Proband No not accurate if sample name does not include family ID, i.e., G series in genome"))


openxlsx::write.xlsx(list("ps>5" = filtered, "all" = OGLanno_count, "readme" = readme), file = output_excel, firstRow = TRUE, firstCol = TRUE)

