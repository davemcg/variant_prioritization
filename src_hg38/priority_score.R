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
## max of (splice_score + in_silico score) = 8, max of in silico is 6: pmin(8, ifelse(PVS1 == 1 | pmaxaf >= 0.03, 0, splice_score) + ifelse(pmaxaf >= 0.02, 0, pmin(6, insilico_score)) ))
library(tidyverse)
#library(vroom)

args <- commandArgs(trailingOnly=TRUE)
#When testing, use the line below.
# setwd("Z:/genome/USUHS/prioritization/temp")
# args <- c("for.ps.tsv", "squirls.csv",
#           "pangolin.tsv", "crossmap.hg19.tsv",
#           "Z:/resources/gnomad/release-2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt",
#           "usuhs44__chr1:1-124478211.ps.tsv")

psInput_file <- args[1]
squirls_file <- args[2]
pangolin_file <- args[3]
crossmap_file <- args[4]
gnomad_metrics_file <- args[5]
psOutput_file <- args[6]

input_df <- read_tsv(psInput_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
  type_convert() %>% 
  mutate(CHROM = as.factor(CHROM)) %>% 
  mutate(across(where(is.character), ~na_if(., "1"))) %>% #bcftools query outputted "1" for blank for character column.
  type_convert()

gnomad_metrics <- read_tsv(gnomad_metrics_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
  type_convert() %>%
  select('gene','pLI','pNull','pRec','oe_lof_upper','syn_z','mis_z','oe_lof_upper_bin','classic_caf','max_af','p') %>% 
  rename(Ref_Gene = gene, LOEUF = oe_lof_upper, sum_lof_af_gnomad = classic_caf, max_lof_af_gnomad = max_af, proportion_pLoF_haplotypes = p)

input_df <- left_join(input_df, gnomad_metrics, by = "Ref_Gene")

squirls_annotation <-  read_csv(squirls_file,  col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "NaN"), col_types = "ficccdc") %>%
  rename(squirls_interpretation = INTERPRETATION, squirls_maxscore = MAX_SCORE, squirls_score = SCORES) %>%
  unite("variantkey", CHROM:ALT, sep = "-", remove = FALSE) %>% 
  group_by(variantkey) %>%
  replace_na(list(squirls_maxscore = 0)) %>%
  slice(which.max(squirls_maxscore)) %>% 
  ungroup(variantkey) 

pangolin <- read_tsv(pangolin_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = "ficcc") %>%
  mutate(pangolin = Pangolin) %>% 
  mutate(across(where(is.character), ~na_if(., "1"))) %>% #bcftools query outputted "1" for blank for character column.
  type_convert() %>% 
  unite("variantkey", CHROM:ALT, sep = "-", remove = TRUE) %>% 
  separate_rows(Pangolin, sep = "ENSG") %>%
  separate(Pangolin, c("temp_pangolin_gene", "temp_pangolin_increase", "temp_pangolin_decrease","temp_warning"), sep = "\\|", remove = FALSE) %>% 
  separate(temp_pangolin_increase, c("temp_pangolin_incr_pos", "temp_pangolin_incr_score"), sep = ":", remove = TRUE, convert = TRUE) %>% 
  separate(temp_pangolin_decrease, c("temp_pangolin_decr_pos", "temp_pangolin_decr_score"), sep = ":", remove = TRUE, convert = TRUE) %>%
  mutate(pangolin_max = pmax(abs(temp_pangolin_incr_score), abs(temp_pangolin_decr_score), na.rm = TRUE)) %>%
  replace_na(list(pangolin_max = 0)) %>% 
  group_by(variantkey) %>%
  slice(which.max(pangolin_max)) %>% 
  ungroup(variantkey) %>% 
  select(-Pangolin, -starts_with("temp_")) 

squirls_pangolin_annotation <- left_join(squirls_annotation, pangolin, by = "variantkey") %>% 
  select(-variantkey)

pickMaxScore <- function(x, y){
  x = as.character(x)
  scores <- as.list(strsplit(x, ";"))[[1]] 
  if (length(scores) > 1) {
    maxscore <- purrr::keep(scores, grepl(format(round(y, 6), nsmall = 6), scores)) 
    if (length(maxscore) == 0) {
      return(x)
    } else if (length(maxscore) == 1) { return(maxscore) }
    else { return(tail(maxscore, n=1)) }
  } else {
    return(x)
  }
}

squirls_pangolin_annotation$squirls_score <- mapply(pickMaxScore, squirls_pangolin_annotation$squirls_score, squirls_pangolin_annotation$squirls_maxscore)

crossmap <- read_tsv(crossmap_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
  type_convert() %>% 
  unite("grch37variant_id", CHROM, POS, REF, ALT, sep = "-")

ps_df_crossmap <- left_join(input_df, crossmap, by = "ID")

rm(input_df)
rm(crossmap)

##The following line was meant to add 3 to Priority_Score when CSQ fields has truncating and PVS1 == 0 and pmaxaf < 0.01 & Priority_Score_intervar < 6
##Empty fields showed as "" for character columns after separate(). Possibly use mutate(across(where(is.character), ~na_if(., "")))
ps_df <-  left_join(ps_df_crossmap, squirls_pangolin_annotation, by=c('CHROM', 'POS', 'REF', 'ALT')) %>%
  mutate(truncating_vep = ifelse(grepl("frameshift_variant|splice_acceptor_variant|splice_donor_variant|start_lost|stop_gained|stop_lost", CSQ, ignore.case = TRUE), 1, 0)) %>% 
  mutate(temp_CSQ = CSQ) %>% 
  separate_rows(temp_CSQ, sep = "\\,") %>%
  separate(temp_CSQ, c('allele','consequence','codons','amino_acids','gene','symbol','MANE_SELECT','feature','exon','intron','hgvsc','hgvsp','max_af','max_af_pops','protein_position','biotype','canonical','domains','existing_variation','clin_sig','pick','pubmed','phenotypes','sift','polyphen','cadd_raw','cadd_phred','genesplicer','spliceregion','MaxEntScan_alt','maxentscan_diff','MaxEntScan_ref','existing_inframe_oorfs','existing_outofframe_oorfs','existing_uorfs','five_prime_utr_variant_annotation','five_prime_utr_variant_consequence','MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','am_class','am_pathogenicity'), sep = "\\|", remove = TRUE, convert = TRUE) %>% 
  #gno2x_af_all,gno3_af_all,maxaf_annovar,gno2x_af_popmax,gno3_popmax,gno_gx_ratio,gno2x_an_all,gno3_an_all,gno2x_filter,gno3_filter AND max_af above from VEP.
  mutate(temp_genesplicer = case_when(grepl("High", genesplicer) ~ 3,
                                      grepl("Medium", genesplicer) ~ 1,
                                      TRUE ~ 0)) %>% 
  mutate(temp_maxentscan_diff =  case_when(abs(maxentscan_diff) > 6 ~ 6,
                                           abs(maxentscan_diff) > 3 ~ 3,
                                           TRUE ~ 0)) %>%
  replace_na(list(max_af=0)) %>%
  mutate(temp_csq_score = ifelse(grepl("deleterious", sift), 0.5, 0) +
           ifelse(grepl("damaging", polyphen), 0.5, 0) +
           ifelse(is.na(cadd_phred), 0, ifelse(cadd_phred > 15, 0.5, 0) ) +
           temp_genesplicer + temp_maxentscan_diff +
           ifelse(five_prime_utr_variant_consequence == "" | max_af > 0.001, 0, 1) +
           ifelse(am_class == "likely_pathogenic", 0.5, 0) ) %>% 
  group_by(ID) %>% slice(which.max(temp_csq_score)) %>% ungroup() %>% 
  mutate(gno2x_expected_an = case_when(CHROM %in% c("X", "chrX") & gno2x_nonpar == "1" ~ 183653,
                                       CHROM %in% c("Y", "chrY") & gno2x_nonpar == "1" ~ 67843,
                                       TRUE ~ 251496)) %>%
  mutate(gno3_expected_an = case_when(CHROM %in% c("X", "chrX") & gno3_nonpar == "1" ~ 116830,
                                       CHROM %in% c("Y", "chrY") & gno3_nonpar == "1" ~ 35482,
                                       TRUE ~ 152312)) %>%
  mutate(gno2x_filter = ifelse(gno2x_an_all > 0 & is.na(gno2x_filter) & gno2x_an_all < gno2x_expected_an/2, "Less50", gno2x_filter),
         gno3_filter = ifelse(gno3_an_all > 0 & is.na(gno3_filter) & gno3_an_all < gno3_expected_an/2, "Less50", gno3_filter) ) %>% 
  replace_na(list(gno2x_af_all=0, gno3_af_all=0, gno2x_af_popmax=0, gno3_maxaf = 0, esp6500siv2_all = 0, f1000g2015aug_all = 0)) %>% 
  mutate(pmaxaf = case_when(is.na(gno2x_filter) & !is.na(gno3_filter) ~ pmax(gno2x_af_all, gno2x_af_popmax, max_af, na.rm = TRUE),
                            !is.na(gno2x_filter) & is.na(gno3_filter) ~ pmax(gno3_af_all, gno3_maxaf, na.rm = TRUE),
                            !is.na(gno2x_filter) & !is.na(gno3_filter) ~ pmax(esp6500siv2_all,f1000g2015aug_all, na.rm = TRUE),
                            TRUE ~ pmax(gno2x_af_all, gno3_af_all, gno2x_af_popmax, gno3_maxaf, max_af, na.rm = TRUE))) %>% 
  replace_na(list(pmaxaf=0)) %>% 
  unite("temp_clinvar", CLNSIG, clin_sig, sep = "-", remove = FALSE ) %>% 
  mutate(temp_clinvar = gsub("_interpretations_of_pathogenicity", "", temp_clinvar)) %>% 
  mutate(temp_existing_variant = gsub("rs\\d+", "", existing_variation)) %>% 
  mutate(temp_existing_variant = gsub("COS\\w\\d+", "", temp_existing_variant)) %>% #VEP100: COSV, #VEP99: COSM
  mutate(temp_existing_variant = gsub("\\w+_\\w+", "", temp_existing_variant)) %>%
  mutate(temp_clinvar_score = ifelse(grepl("pathogenic", temp_clinvar, ignore.case = TRUE), 6, 0) - 
           ifelse(grepl("benign", temp_clinvar, ignore.case = TRUE) & grepl("pathogenic", temp_clinvar, ignore.case = TRUE), 3, 0) ) %>%
  mutate(temp_hgmd_score = case_when( hgmd_class == "DM" ~ 3,
                                      grepl("[A-Z]", temp_existing_variant) | !is.na(HGMD_Overlap4aa) & pmaxaf < 0.005 ~ 1,
                                      hgmd_class == "DM\\?" ~ 0.5, 
                                      TRUE ~ 0)) %>% 
  mutate(clinvar_hgmd_score = pmin(6, temp_clinvar_score + temp_hgmd_score) ) %>% 
  mutate(insilico_score = ifelse(is.na(ClinPred_score), 0, ifelse(ClinPred_score > 0.5, 0.5, 0)) + 
           ifelse(is.na(REVEL_score), 0, ifelse(REVEL_score > 0.55, 0.5, 0)) +
           ifelse(grepl("deleterious", sift), 0.5, 0) +
           ifelse(grepl("damaging", polyphen), 0.5, 0) + 
           #ifelse(grepl("D", SIFT_pred), 0.5, 0) + #T for tolerant, D for damaging
           #ifelse(grepl("D", Polyphen2_HVAR_pred), 0.5, 0) + #B for benign, P for possibly damaging, D for probably damgaing in dbSNP4
           ifelse(is.na(mutscore), 0, ifelse(mutscore > 0.7, 0.5, 0)) +
           ifelse(grepl("D", MetaSVM_pred), 0.5, 0) + 
           ifelse(is.na(PrimateAI_score), 0, ifelse(PrimateAI_score > 0.803, 0.5, 0)) +
           ifelse(is.na(MPC_score), 0, ifelse(MPC_score > 1.5 & pmaxaf < 0.02, 0.5, 0)) +
           ifelse(grepl("H|M", MutationAssessor_pred), 0.5, 0) + 
           ifelse(grepl("D", MutationTaster_pred), 0.5, 0) +
           ifelse(grepl("D", PROVEAN_pred), 0.5, 0) +
           ifelse(is.na(Eigen_PC_raw_coding), 0, ifelse(Eigen_PC_raw_coding > 0, 0.5, 0)) + 
           ifelse(is.na(cadd_phred), 0, ifelse(cadd_phred > 15, 0.5, 0) ) +
           ifelse(is.na(phyloP100way_vertebrate), 0, ifelse(phyloP100way_vertebrate > 2, 0.5, 0)) + 
           ifelse(is.na(GERPplus_RS), 0, ifelse(GERPplus_RS > 1, 0.5, 0)) +
           ifelse(is.na(ccr_pct) | is.na(ExonicFunc_refGeneWithVer), 0, ifelse(ccr_pct > 0.95 & grepl("nonframeshift|nonsynonymous", ExonicFunc_refGeneWithVer) & pmaxaf < 0.02, 0.5, 0)) +
           ifelse(is.na(remm), 0, ifelse(remm > 0.6 & pmaxaf < 0.02, 0.5, 0)) + #fine above
           ifelse(is.na(fathmm_XF_coding_score), 0, ifelse(fathmm_XF_coding_score > 0.6 & pmaxaf < 0.02, 0.5, 0)) +
           ifelse(is.na(fathmm_xf_noncoding), 0, ifelse(fathmm_xf_noncoding > 0.6 & pmaxaf < 0.02, 0.5, 0)) + 
           ifelse(is.na(hmc_score), 0, ifelse(hmc_score < 0.8, 0.5, 0)) +
           ifelse(is.na(gnomad_nc_constraint), 0, ifelse(gnomad_nc_constraint > 4 & pmaxaf < 0.01, 0.5, 0)) +
           ifelse(am_class == "likely_pathogenic", 0.5, 0)) %>% 
  replace_na(list(clinvar_hgmd_score=0, insilico_score=0)) %>% 
  mutate(temp_dpsi_max_tissue = dpsi_max_tissue) %>% 
  mutate(temp_dpsi_zscore = dpsi_zscore) %>%
  mutate(temp_dbscsnv_ada_score = dbscSNV_ADA_SCORE) %>% 
  mutate(temp_dbscsnv_rf_score = dbscSNV_RF_SCORE) %>% 
  replace_na(list(temp_dpsi_max_tissue = 0, temp_dpsi_zscore = 0, temp_dbscsnv_ada_score = 0, temp_dbscsnv_rf_score = 0 )) %>%
  mutate(temp_dpsi_score = case_when(pmaxaf < 0.02 & abs(temp_dpsi_max_tissue + temp_dpsi_zscore) > 10 ~ 3,
                                     pmaxaf < 0.02 & abs(temp_dpsi_max_tissue + temp_dpsi_zscore) > 5 ~ 1,
                                     TRUE ~ 0)) %>%
  mutate(temp_dbscSNV_score = ifelse((pmaxaf < 0.02 & temp_dbscsnv_ada_score>0.8 & temp_dbscsnv_rf_score>0.5), 3, 0)) %>% 
  mutate(temp_squirl_score = case_when(squirls_interpretation == "pathogenic" & squirls_maxscore > 0.8 ~ 3, #according to squirls paper Fig. S4.
                                       squirls_interpretation == "pathogenic" & squirls_maxscore > 0.2 ~ 0.5,
                                       TRUE ~ 0)) %>%
  mutate(temp_pangolin_score = case_when( pangolin_max > 0.5 ~ 6,
                                          pangolin_max > 0.2 ~ 3,
                                          TRUE ~ 0)) %>% 
  mutate(splice_score = pmin(8, (spliceai_rank + temp_genesplicer + temp_maxentscan_diff + temp_dpsi_score + temp_dbscSNV_score + temp_squirl_score + 
                                   temp_pangolin_score + ifelse(pmaxaf < 0.02 & !is.na(branchpoint_prob), 2, 0) + ifelse(pmaxaf < 0.02 & !is.na(labranchor_score), 2, 0) )) ) %>% 
  mutate(priority_score = PrScore_intervar + ifelse(pmaxaf < 0.03, clinvar_hgmd_score, 0) +  
           ifelse(PVS1 == 1 | pmaxaf >= 0.005 | PrScore_intervar > 6 | splice_score > 2, 0, truncating_vep*3) +
           pmin(8, ifelse(PVS1 == 1 | pmaxaf >= 0.03, 0, splice_score) + ifelse(pmaxaf >= 0.02, 0, pmin(6, insilico_score)) )) %>% 
  mutate(other_modification = ifelse(grepl("protein_altering_variant|inframe", CSQ, ignore.case = TRUE) & pmaxaf < 0.01 & PrScore_intervar < 5 & insilico_score < 3, 3, 0) +
           ifelse(grepl("missense_variant", CSQ, ignore.case = TRUE) & mis_z >= 3.09 & SigmaAF_Missense_0001 < 0.005 & pmaxaf < 0.005 & insilico_score < 4, 2, 0) +
           ifelse(grepl("upstream|downstream|UTR", Func_refGeneWithVer) & pmaxaf < 0.005, 0.5, 0) +
           ifelse(five_prime_utr_variant_consequence == "" | pmaxaf > 0.001, 0, 1) +
           ifelse(is.na(regsnp_disease) | regsnp_disease == "B" | pmaxaf > 0.001, 0, ifelse(regsnp_disease == "D", 1, 0.5)) +
           ifelse(grepl("0,255,0|255,0,0", atac_rpe_itemRgb) & pmaxaf < 0.001 & priority_score < 5, 1, 0) + 
           ifelse(is.na(ft_ret_rpe_score), 0, ifelse(pmaxaf < 0.001 & priority_score < 5, 0.5, 0)) +
           ifelse(is.na(cherry_sum_score) | cherry_sum_score > -5, 0, ifelse(pmaxaf < 0.001 & priority_score < 5, 1, 0)) ) %>% 
  mutate(priority_score = priority_score + other_modification) %>% 
  mutate(priority_score = ifelse(priority_score < 5 & Ref_Gene == "ROM1" & truncating_vep == "1", 5, priority_score)) %>%
  #mutate(priority_score = ifelse(priority_score < 5 & pmaxaf < 0.005 & is.na(Ref_Gene) & (!is.na(gene_gnomad) | !is.na(eyeIntegration_gene) | !is.na(omim_Gene)), priority_score + 2, priority_score)) %>% 
  mutate(priority_score = case_when(priority_score < 5 & pmaxaf < 0.001 & Ref_Gene %in% c("MIR184","MIR204") ~ 5,
                                    priority_score < 4.5 & pmaxaf < 0.001 & grepl("^MIR", Ref_Gene, ignore.case = TRUE) ~ 4,
                                    TRUE ~ priority_score)) %>% 
  select(CHROM, POS, REF, ALT, priority_score, clinvar_hgmd_score, splice_score, insilico_score, pmaxaf, truncating_vep, squirls_interpretation, squirls_maxscore, squirls_score, pangolin, grch37variant_id,
         'pLI','pNull','pRec','LOEUF','syn_z','mis_z','oe_lof_upper_bin','sum_lof_af_gnomad','max_lof_af_gnomad','proportion_pLoF_haplotypes') %>% 
  arrange(CHROM, POS)

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
#write_tsv(ps_df, file = psOutput_file)
write_tsv(ps_df, file.path('.',  psOutput_file), na=".")
