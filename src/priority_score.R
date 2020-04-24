#updated 7/28/19
## dbscSNV's ada > 0.8 and rf>0.5, Score=+3 (equals weight of PM)
## |dpsi_max_tissue+dpsi_zscore| > 10, score=+3; > 5, score=+1 make a histogram of dpsi of WGS set and determine the cut-off
## genesplicer (H|M) or maxentscan_diff > 3, score=+3
## spliceai_rank >0.8, score=+8; >0.5, score=+6; >0.2, score=+3; >0.15, score=+1; calculated in snakemake file. make a histogram of splicai score and determine the cut-off
## splice_score = min(8, spliceai etc)
## if PVS == 1 or maxaf > 0.02, then splice score is not added to the  priority score.
## other_pred_score is not added to priority score if maxaf > 0.02
## if impact == "missense_variant" & mis_z >= 3.09 & SigmaAF_Missense_0001 < 0.005 & pmaxaf < 0.0005, priority socre += 2
## inframe indels are added by intervar (3 pts); (protein_altering in vep impact could include mean "inframe indel" if delins): priority score += 3 if pmaxaf < 0.01 & Priority_Score_intervar < 5 (PM)

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
#args <- c("Y:/NextSeqAnalysis/training_run/freebayesPrioritization0218/temp/20200215v1.freebayes__1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y.tsv", "OGLv1_panel_DxORcandidate.txt", "1368.sorted.tsv", "1368.filtered.tsv")

input_df <- read_tsv(args[1], col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert()
##The following line was meant to add 3 to Priority_Score when CSQ fields has truncating and PVS1 == 0 and pmaxaf < 0.01 & Priority_Score_intervar < 6
ps_df <-  input_df %>% mutate(truncating_vep = ifelse(grepl("frameshift_variant|splice_acceptor_variant|splice_donor_variant|start_lost|stop_gained|stop_lost", CSQ, ignore.case = TRUE), 1, 0)) %>% 
  mutate(temp_CSQ = sub(",.*", "", CSQ)) %>%
  separate(temp_CSQ, c('allele','consequence','codons','amino_acids','gene','symbol','feature','exon','hgvsc','hgvsp','max_af','max_af_pops','polyphen','sift','mpc','protein_position','biotype','canonical','domains','existing_variation','clin_sig','pick','pubmed','phenotypes','cadd_raw','cadd_phred','genesplicer','spliceregion','ada_score','rf_score','maxentscan_diff'), sep = "\\|", remove = TRUE, convert = TRUE) %>% 
  replace_na(list(gnomAD_exome_ALL_annovar=0, gnomAD_genome_ALL_annovar=0, PopFreqMax_annovar=0, max_af=0, mpc=0, gno_af_popmax=0, mis_z=0)) %>% 
  mutate(pmaxaf = pmax(gnomAD_exome_ALL_annovar, gnomAD_genome_ALL_annovar, na.rm = TRUE)) %>% #removed PopFreqMax_annovar, max_af, gno_af_popmax, 4/14/2020
  unite("temp_clinvar", clinvar_sig, clin_sig, Clinvar_intervar, sep = "-", remove = FALSE ) %>% 
  mutate(temp_clinvar = gsub("_interpretations_of_pathogenicity", "", temp_clinvar)) %>% 
  mutate(temp_existing_variant = gsub("rs\\d+", "", existing_variation)) %>% 
  mutate(temp_existing_variant = gsub("COSM\\d+", "", temp_existing_variant)) %>%
  mutate(temp_existing_variant = gsub("\\w+_\\w+", "", temp_existing_variant)) %>%
  mutate(clinvar_hgmd_score = ifelse(grepl("pathogenic", temp_clinvar, ignore.case = TRUE), 6, 0) - 
           ifelse(grepl("benign", temp_clinvar, ignore.case = TRUE) & grepl("pathogenic", temp_clinvar, ignore.case = TRUE), 3, 0) + 
           ifelse(grepl("[A-Z]", HGMD_Overlap) | grepl("[A-Z]", temp_existing_variant), 3, 0) ) %>% 
  mutate(clinvar_hgmd_score = ifelse(clinvar_hgmd_score > 6, 6, clinvar_hgmd_score)) %>% 
  mutate(other_predic_score = ifelse(is.na(ClinPred_Score), 0, ifelse(ClinPred_Score > 0.5, 0.5, 0)) + ifelse(grepl("deleterious", sift), 0.5, 0) +
           ifelse(grepl("damaging", polyphen), 0.5, 0) + ifelse(grepl("D", MetaSVM_pred), 0.5, 0) + 
           ifelse(is.na(PrimateDL), 0, ifelse(PrimateDL > 0.803, 0.5, 0)) +
           ifelse(is.na(mpc), 0, ifelse(mpc > 1.5 & pmaxaf < 0.02, 0.5, 0)) +
           ifelse(grepl("H|M", MutationAssessor_pred), 0.5, 0) + ifelse(grepl("D", MutationTaster_pred), 0.5, 0) + ifelse(grepl("D", PROVEAN_pred), 0.5, 0) +
           ifelse(is.na(Eigen_PC_raw), 0, ifelse(Eigen_PC_raw > 0 | Eigen_raw > 0, 0.5, 0)) + 
           ifelse(is.na(cadd_phred), 0, ifelse(cadd_phred > 15, 0.5, 0) ) +
           ifelse(is.na(phyloP_100way), 0, ifelse(phyloP_100way > 2, 0.5, 0)) +
           ifelse(is.na(GERP_RS_intervar), 0, ifelse(GERP_RS_intervar> 1, 0.5, 0))) %>% 
  replace_na(list(clinvar_hgmd_score=0, other_predic_score=0)) %>% 
  mutate(temp_genesplicer = case_when(grepl("High", genesplicer) ~ 3,
                                      grepl("Medium", genesplicer) ~ 1,
                                      TRUE ~ 0)) %>% 
  mutate(temp_maxentscan_diff =  case_when(abs(maxentscan_diff) > 6 ~ 6,
                                           abs(maxentscan_diff) > 3 ~ 3,
                                           TRUE ~ 0)) %>%
  mutate(temp_dpsi_max_tissue = dpsi_max_tissue_annovar) %>% 
  mutate(temp_dpsi_zscore = dpsi_zscore_annovar)%>%
  mutate(temp_dbscsnv_ada_score = dbscSNV_ADA_SCORE_intervar) %>% 
  mutate(temp_dbscsnv_rf_score = dbscSNV_RF_SCORE_intervar) %>% 
  replace_na(list(temp_dpsi_max_tissue = 0, temp_dpsi_zscore = 0, temp_dbscsnv_ada_score = 0, temp_dbscsnv_rf_score = 0 )) %>%
  mutate(temp_dpsi_score = case_when(pmaxaf < 0.02 & abs(temp_dpsi_max_tissue + temp_dpsi_zscore) > 10 ~ 3,
                                     pmaxaf < 0.02 & abs(temp_dpsi_max_tissue + temp_dpsi_zscore) > 5 ~ 1,
                                     TRUE ~ 0)) %>%
  mutate(temp_dbscSNV_score = ifelse((pmaxaf < 0.02 & temp_dbscsnv_ada_score>0.8 & temp_dbscsnv_rf_score>0.5), 3, 0)) %>% 
  mutate(splice_score = pmin(8, (spliceai_rank + temp_genesplicer + temp_maxentscan_diff + temp_dpsi_score + temp_dbscSNV_score)) ) %>% 
  mutate(priority_score = Priority_Score_intervar + ifelse(pmaxaf < 0.03, clinvar_hgmd_score, 0) + ifelse(PVS1 == 1 | pmaxaf >= 0.03, 0, splice_score) + 
           ifelse(PVS1 == 1 | pmaxaf >= 0.005 | Priority_Score_intervar > 6 | splice_score > 2, 0, truncating_vep*3) +
           ifelse(pmaxaf >= 0.02, 0, other_predic_score) +
           ifelse(grepl("protein_altering_variant", CSQ, ignore.case = TRUE) & pmaxaf < 0.01 & Priority_Score_intervar < 5, 3, 0) +
           ifelse(grepl("missense_variant", CSQ, ignore.case = TRUE) & mis_z >= 3.09 & SigmaAF_Missense_0001 < 0.005 & pmaxaf < 0.005, 2, 0)) %>% 
  mutate(priority_score = ifelse(priority_score < 5 & Ref_Gene_annovar == "ROM1" & truncating_vep == "1", 5, priority_score)) %>% 
  select(CHROM, POS, REF, ALT, priority_score, clinvar_hgmd_score, splice_score, other_predic_score, pmaxaf, truncating_vep)
#pmaxaf cutoff increased to 0.005 from 0.0005; 4/14/2020
#http://web.corral.tacc.utexas.edu/WGSAdownload/resources/dbNSFP/dbNSFP4.0b2c.readme.txt
#CADD: https://cadd.gs.washington.edu/info
#eigen:
#http://mutationassessor.org/r3/howitworks.php

#If all_AF<0.002, add 2 to priority score, if < 0.005, add 1 ???

#Read in OGLv1_panel gene class of either Dx or Candidate

#sampleData <- read.tsv(args[1], header = TRUE, check.names=FALSE ) #nrows = 5
# classes <- sapply(gemini, class)
# largeData <- read.csv("huge-file.csv", header = TRUE, colClasses = classes)
#write_tsv(ps_df, path = './test.tsv')
write_tsv(ps_df, path = args[2])
