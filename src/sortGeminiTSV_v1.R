#updated 7/28/19
## dbscSNV's ada > 0.8 and rf>0.5, Score=+3 (equals weight of PM)
## |dpsi_max_tissue+dpsi_zscore| > 5, score=+3; > 2.5, score=+1 make a histogram of dpsi of WGS set and determine the cut-off
## genesplicer (H|M) or maxentscan_diff > 3, score=+3
## spliceai_rank >0.8, score=+8; >0.5, score=+6; >0.2, score=+3; >0.15, score=+1; calculated in snakemake file. make a histogram of splicai score and determine the cut-off
## splice_score = min(8, spliceai etc)
## if PVS == 1 or maxaf > 0.02, then splice score is not added to the  priority score.
## other_pred_score is not added to priority score if maxaf > 0.02

args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
#args <- c("1368_D332_001.gemini.tsv", "OGLv1_panel_DxORcandidate.txt", "1368.sorted.tsv", "1368.filtered.tsv")

library(tidyverse)

gemini_input <- read_tsv(args[1], col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert() %>% mutate(exon = sub("^", " ", exon))

gemini <-  gemini_input %>% mutate( temp_start_vcf = start + 1 ) %>% 
  unite("chr_variant_id", chrom, temp_start_vcf, ref, alt, sep = "-", remove = FALSE ) %>%
  replace_na(list(gnomad_exome_all_annovar=0, gnomad_genome_all_annovar=0, popfreqmax_annovar=0, max_af=0, gno_af_popmax=0, max_aaf_all=0)) %>% 
  mutate(maxaf_postgemini = pmax(gnomad_exome_all_annovar, gnomad_genome_all_annovar, popfreqmax_annovar, max_af, gno_af_popmax, max_aaf_all, na.rm = TRUE)) %>% 
  unite("temp_clinvar", clinvar_sig, clin_sig, clinvar_intervar, sep = "-", remove = FALSE ) %>% 
  mutate(temp_clinvar = gsub("_interpretations_of_pathogenicity", "", temp_clinvar)) %>% 
  mutate(temp_existing_variant = gsub("rs\\d+", "", existing_variation)) %>% 
  mutate(temp_existing_variant = gsub("COSM\\d+", "", temp_existing_variant)) %>%
  mutate(temp_existing_variant = gsub("\\w+_\\w+", "", temp_existing_variant)) %>%
  mutate(clinvar_hgmd_score = ifelse(grepl("pathogenic", temp_clinvar, ignore.case = TRUE), 6, 0) - 
           ifelse(grepl("benign", temp_clinvar, ignore.case = TRUE) & grepl("pathogenic", temp_clinvar, ignore.case = TRUE), 3, 0) + 
           ifelse(grepl("[A-Z]", hgmd_overlap) | grepl("[A-Z]", temp_existing_variant), 3, 0) ) %>% 
  mutate(clinvar_hgmd_score = ifelse(clinvar_hgmd_score > 6, 6, clinvar_hgmd_score)) %>% 
  mutate(clinvar_hgmd_score = ifelse(maxaf_postgemini < 0.02, clinvar_hgmd_score, 0)) %>% 
  mutate(other_predic_score = ifelse(is.na(clinpred_score), 0, ifelse(clinpred_score > 0.5, 0.5, 0)) + ifelse(grepl("deleterious", sift_pred), 0.5, 0) +
           ifelse(grepl("damaging", polyphen_pred), 0.5, 0) + ifelse(grepl("D", metasvm_pred), 0.5, 0) + 
           ifelse(is.na(primatedl), 0, ifelse(primatedl > 0.803, 0.5, 0)) +
           ifelse(is.na(mpc), 0, ifelse(mpc > 1.5 & maxaf_postgemini < 0.02, 0.5, 0)) +
           ifelse(grepl("H|M", mutationassessor_pred), 0.5, 0) + ifelse(grepl("D", mutationtaster_pred), 0.5, 0) + ifelse(grepl("D", provean_pred), 0.5, 0) +
           ifelse(is.na(eigen_pc_raw), 0, ifelse(eigen_pc_raw > 0 | eigen_raw > 0, 0.5, 0)) + 
           ifelse(is.na(cadd_phred), 0, ifelse(cadd_phred > 15, 0.5, 0) ) +
           ifelse(is.na(phylop_100way), 0, ifelse(phylop_100way > 2, 0.5, 0)) +
           ifelse(is.na(gerp_rs_intervar), 0, ifelse(gerp_rs_intervar> 1, 0.5, 0))) %>% 
  replace_na(list(clinvar_hgmd_score=0, other_predic_score=0)) %>% 
  mutate(temp_genesplicer_maxent_score = ifelse((!is.na(genesplicer) | (!is.na(maxentscan_diff) & abs(maxentscan_diff) > 3)) & maxaf_postgemini < 0.02, 3, 0)) %>% 
  mutate(temp_dpsi_max_tissue = dpsi_max_tissue_annovar) %>% 
  mutate(temp_dpsi_zscore = dpsi_zscore_annovar)%>%
  mutate(temp_dbscsnv_ada_score = dbscsnv_ada_score_intervar) %>% 
  mutate(temp_dbscsnv_rf_score = dbscsnv_rf_score_intervar) %>% 
  replace_na(list(temp_dpsi_max_tissue = 0, temp_dpsi_zscore = 0, temp_dbscsnv_ada_score = 0, temp_dbscsnv_rf_score = 0 )) %>%
  mutate(temp_dpsi_score = ifelse(maxaf_postgemini < 0.02 & abs(temp_dpsi_max_tissue + temp_dpsi_zscore) > 5, 3, (ifelse((maxaf_postgemini < 0.02 & abs(temp_dpsi_max_tissue + temp_dpsi_zscore) > 2.5), 1, 0)))  ) %>%
  mutate(temp_dbscSNV_score = ifelse((maxaf_postgemini < 0.02 & temp_dbscsnv_ada_score>0.8 & temp_dbscsnv_rf_score>0.5), 3, 0)) %>% 
  mutate(splice_score = pmin(8, (spliceai_rank + temp_genesplicer_maxent_score + temp_dpsi_score + temp_dbscSNV_score)) ) %>% 
  mutate(priority_score = priority_score_intervar + clinvar_hgmd_score + ifelse(pvs1 == 1 | maxaf_postgemini >= 0.02, 0, splice_score) + 
           ifelse(maxaf_postgemini >= 0.02, 0, other_predic_score) )


#http://web.corral.tacc.utexas.edu/WGSAdownload/resources/dbNSFP/dbNSFP4.0b2c.readme.txt
#CADD: https://cadd.gs.washington.edu/info
#eigen:
#http://mutationassessor.org/r3/howitworks.php

#If all_AF<0.002, add 2 to priority score, if < 0.005, add 1 ???

#Read in OGLv1_panel gene class of either Dx or Candidate

#sampleData <- read.tsv(args[1], header = TRUE, check.names=FALSE ) #nrows = 5
# classes <- sapply(gemini, class)
# largeData <- read.csv("huge-file.csv", header = TRUE, colClasses = classes)


OGLv1_gene_class <- read.delim(args[2], sep = "\t", header = T, colClasses = c("character","character") )

# get max_priority_score for each gene
max_priority_score <- select(gemini, c(ref_gene_annovar, priority_score)) %>% group_by(ref_gene_annovar) %>% summarize(maxpriorityscore = max(priority_score)) 

#arrange by max_priority_score, then by gene, and priority score. None gene region?
gemini_max_priority_score <- left_join(gemini, max_priority_score, by=c("ref_gene_annovar"))

#VEP hg19 version's gene names are the same as in the IDT ordering design sheets. This is what used for left_join
gemini_sorted <- left_join(gemini_max_priority_score, OGLv1_gene_class, by=c("gene")) %>% 
  arrange(desc(maxpriorityscore), ref_gene_annovar, desc(priority_score)) %>% 
  select(-starts_with('temp_')) %>% 
  select('chr_variant_id','chr_annovar', 'start_annovar', 'ref_annovar', 'alt_annovar', 'qual', 'filter', starts_with('gts'), starts_with('gt_'), 
         'panel_class', 'priority_score', 'priority_score_intervar', 'clinvar_hgmd_score', 'splice_score', 'other_predic_score', 'maxaf_postgemini', 'exac_num_hom_alt', 'gno_hom',  'gene_refgenewithver_annovar', 
         'genedetail_refgenewithver_annovar', 'exonicfunc_refgenewithver_annovar', 'aachange_refgenewithver_annovar', 'hgvsc', 'hgvsp', 'omim_genesymbol', 'omim_inheritance', 'omim_phenotypes', 'gene', 'exon', 'aa_length', 'ref_gene_annovar', 'func_refgene_annovar', 'hgmd_overlap', 'existing_variation', 'clinvar_intervar', 'intervar_and_evidence', 
         'clinvar_id', 'clinvar_pathogenic', 'clinvar_sig', 'clin_sig', 'gnomad_exome_all_annovar', 'gnomad_genome_all_annovar', 'popfreqmax_annovar', 'max_af', 'max_af_pops',  
         'spliceai', 'spliceai_maxscore', 'spliceai_filtered', 'dbscsnv_ada_score_intervar', 'dbscsnv_rf_score_intervar', 'dpsi_max_tissue_annovar', 'dpsi_zscore_annovar', 'genesplicer', 'maxentscan_diff', 'branchpoint_u2_binding_energy', 'branchpoint_prob', 'branchpoint_to_3prime', 
         'sift_pred', 'polyphen_pred', 'mutationassessor_pred', 'mutationtaster_pred', 'metasvm_pred','metasvm_score_intervar', 'clinpred_score', 'mpc', 'cadd_raw', 'cadd_phred', 'eigen_pc_raw', 'eigen_raw', 'gerp_rs_intervar', 'phylop46way_placental_intervar', 'phylop_100way', 'primatedl', 'func_refgenewithver_annovar', 'exonicfunc_refgene_intervar', 'avsnp150_annovar', 'interpro_domain_intervar', 
         'pfam_domain', 'tfbs', 'pli', 'lof_z', 'mis_z', 'sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01', 'atac_rpe_itemrgb', 'atac_rpe_score', 
         'eyeintegration_rnaseq_tpm_rpe_adulttissue', 'eyeintegration_rnaseq_tpm_rpe_cellline', 'eyeintegration_rnaseq_tpm_rpe_fetaltissue', 'eyeintegration_rnaseq_tpm_rpe_stemcellline', 'eyeintegration_rnaseq_tpm_retina_adulttissue', 'eyeintegration_rnaseq_tpm_retina_stemcellline', 'eyeintegration_rnaseq_tpm_wholeblood', 
         'chrom', 'start', 'end', 'ref', 'alt', 'freq_esp6500siv2_all_annovar', 'freq_1000g2015aug_all_annovar', 'aaf_esp_all', 'aaf_1kg_all', 'af_exac_all', 'pubmed', everything() )

write_tsv(gemini_sorted, file.path('.', args[3]))

gemini_filtered <- gemini_sorted %>% mutate(temp_group = ifelse(priority_score >= 3, 3, ifelse(priority_score >= -2, -2, -3))) %>% 
  filter(temp_group >= -2, maxaf_postgemini < 0.2) %>% arrange(desc(temp_group), desc(maxpriorityscore), ref_gene_annovar, desc(priority_score)) %>% 
  select(-maxpriorityscore, -temp_group)


write_tsv(gemini_filtered, file.path('.', args[4]))
