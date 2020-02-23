#updated 7/28/19
## dbscSNV's ada > 0.8 and rf>0.5, Score=+3 (equals weight of PM)
## |dpsi_max_tissue+dpsi_zscore| > 6, score=+3; > 3, score=+1 make a histogram of dpsi of WGS set and determine the cut-off
## genesplicer (H|M) or maxentscan_diff > 3, score=+3
## spliceai_rank >0.8, score=+8; >0.5, score=+6; >0.2, score=+3; >0.15, score=+1; calculated in snakemake file. make a histogram of splicai score and determine the cut-off
## splice_score = min(8, spliceai etc)
## if PVS == 1 or maxaf > 0.02, then splice score is not added to the  priority score.
## other_pred_score is not added to priority score if maxaf > 0.02
## if impact == "missense_variant" & mis_z >= 3.09 & sigmaaf_missense_0001 < 0.005 & pmaxaf < 0.0005, priority socre += 2

args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
# args <- c("Y:/NextSeqAnalysis/trainingRun/freebayesPrioritization/gemini_tsv/G05068.20200215v1.freebayes.training_RPGR_nano.gemini.tsv",
#           "Z:/OGL_NGS/variant_prioritization/data/OGLv1_panel_DxORcandidate.txt", "sorted.tsv", "filtered.tsv", "G05068", "filtered.xlsx")

library(tidyverse)
library(readxl)

gemini_input <- read_tsv(args[1], col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert() %>% mutate(exon = sub("^", " ", exon))

gemini <-  gemini_input %>% mutate( temp_start_vcf = start + 1 ) %>% 
  unite("chr_variant_id", chrom, temp_start_vcf, ref, alt, sep = "-", remove = FALSE ) %>% 
  mutate(gene = toupper(gene))
  

#http://web.corral.tacc.utexas.edu/WGSAdownload/resources/dbNSFP/dbNSFP4.0b2c.readme.txt
#CADD: https://cadd.gs.washington.edu/info
#eigen:
#http://mutationassessor.org/r3/howitworks.php

#If all_AF<0.002, add 2 to priority score, if < 0.005, add 1 ???

#Read in OGLv1_panel gene class of either Dx or Candidate

#sampleData <- read.tsv(args[1], header = TRUE, check.names=FALSE ) #nrows = 5
# classes <- sapply(gemini, class)
# largeData <- read.csv("huge-file.csv", header = TRUE, colClasses = classes)


OGLv1_gene_class <- read.delim(args[2], sep = "\t", header = T, colClasses = c("character","character","character") ) %>% 
  select('gene', 'panel_class')

# get max_priority_score for each gene
max_priority_score <- select(gemini, c(ref_gene_annovar, priority_score)) %>% group_by(ref_gene_annovar) %>% summarize(maxpriorityscore = max(priority_score)) 

#arrange by max_priority_score, then by gene, and priority score. None gene region?
gemini_max_priority_score <- left_join(gemini, max_priority_score, by=c("ref_gene_annovar"))

#VEP hg19 version's gene names are the same as in the IDT ordering design sheets. This is what used for left_join
gemini_sorted <- left_join(gemini_max_priority_score, OGLv1_gene_class, by = c("gene")) %>% 
  mutate(note = "NA") %>% 
  arrange(desc(maxpriorityscore), ref_gene_annovar, desc(priority_score)) %>% 
  select(-starts_with('temp_')) %>% 
  mutate(note = "") %>% 
  select('chr_variant_id','chr_annovar', 'start_annovar', 'ref_annovar', 'alt_annovar', 'qual', 'filter', starts_with('gts'), starts_with('gt_'), 'aaf',
         'panel_class', 'priority_score', 'priority_score_intervar', 'clinvar_hgmd_score', 'splice_score', 'other_predic_score', 'pmaxaf', 'exac_num_hom_alt', 'gno_hom', 'pvs1', 'truncating_vep', 'gene_refgenewithver_annovar', 'note', 
         'genedetail_refgenewithver_annovar', 'exonicfunc_refgenewithver_annovar', 'aachange_refgenewithver_annovar', 'hgvsc', 'hgvsp', 'omim_genesymbol', 'omim_inheritance', 'omim_phenotypes', 'gene', 'exon', 'aa_length', 'ref_gene_annovar', 'func_refgene_annovar', 'hgmd_overlap', 'existing_variation', 'clinvar_intervar', 'intervar_and_evidence', 
         'clinvar_id', 'clinvar_pathogenic', 'clinvar_sig', 'clin_sig', 'gnomad_exome_all_annovar', 'gnomad_genome_all_annovar', 'popfreqmax_annovar', 'max_af', 'max_af_pops',  
         'spliceai', 'spliceai_maxscore', 'spliceai_filtered', 'dbscsnv_ada_score_intervar', 'dbscsnv_rf_score_intervar', 'dpsi_max_tissue_annovar', 'dpsi_zscore_annovar', 'genesplicer', 'maxentscan_diff', 'branchpoint_u2_binding_energy', 'branchpoint_prob', 'branchpoint_to_3prime', 
         'sift_pred', 'polyphen_pred', 'mutationassessor_pred', 'mutationtaster_pred', 'metasvm_pred','metasvm_score_intervar', 'clinpred_score', 'mpc', 'cadd_raw', 'cadd_phred', 'eigen_pc_raw', 'eigen_raw', 'gerp_rs_intervar', 'phylop46way_placental_intervar', 'phylop_100way', 'primatedl', 'func_refgenewithver_annovar', 'exonicfunc_refgene_intervar', 'avsnp150_annovar', 'interpro_domain_intervar', 
         'pfam_domain', 'tfbs', 'pli', 'lof_z', 'mis_z', 'sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01', 'atac_rpe_itemrgb', 'atac_rpe_score', 
         'eyeintegration_rnaseq_tpm_rpe_adulttissue', 'eyeintegration_rnaseq_tpm_rpe_cellline', 'eyeintegration_rnaseq_tpm_rpe_fetaltissue', 'eyeintegration_rnaseq_tpm_rpe_stemcellline', 'eyeintegration_rnaseq_tpm_retina_adulttissue', 'eyeintegration_rnaseq_tpm_retina_stemcellline', 'eyeintegration_rnaseq_tpm_wholeblood', 
         'chrom', 'start', 'end', 'ref', 'alt', 'freq_esp6500siv2_all_annovar', 'freq_1000g2015aug_all_annovar', 'aaf_esp_all', 'aaf_1kg_all', 'af_exac_all', 'pubmed', everything() )

write_tsv(gemini_sorted, path = args[3])

gemini_filtered <- gemini_sorted %>% mutate(temp_group = ifelse(priority_score >= 3, 3, ifelse(priority_score >= -2, -2, -3))) %>% 
  filter(temp_group >= -2, pmaxaf < 0.2, aaf < args[7]) %>% arrange(desc(temp_group), desc(maxpriorityscore), ref_gene_annovar, desc(priority_score)) %>% 
  select(-maxpriorityscore, -temp_group)


write_tsv(gemini_filtered, path = args[4])

gemini_filtered2 <- gemini_sorted %>% filter(priority_score >= 4, pmaxaf < 0.2, aaf < args[7]) %>% rename_all(funs(str_replace(., args[5], ""))) 

recessive_count <- select(gemini_filtered2, c(ref_gene_annovar, priority_score, gt_types.)) %>%
  filter(priority_score >= 5) %>% select(-priority_score) %>% 
  group_by(ref_gene_annovar) %>% summarize(recessive_cnt = sum(gt_types.)) 

gemini_filtered3 <- left_join(gemini_filtered2, recessive_count, by=c("ref_gene_annovar")) %>% 
  replace_na(list(recessive_cnt=0)) %>% 
  mutate(recessive_cnt = as.integer(recessive_cnt))

xR <- gemini_filtered3 %>% filter(chr_annovar == "X", recessive_cnt >= 2) %>% select(-recessive_cnt)
xD <- gemini_filtered3 %>% filter(chr_annovar == "X", recessive_cnt == 1, pmaxaf < 0.002) %>% select(-recessive_cnt)

#ar are those genes with homozygous or compound hets variants of ps >= 5. However, ps = 4 variants were also listed if there are 2 ps>=5.
#ar gene with 1 hit will not be here.
ar <- gemini_filtered3 %>% filter(!chr_annovar %in% c("X", "Y"), !omim_inheritance %in% c("AD"), recessive_cnt >= 2) %>% 
  mutate(knownAR = ifelse(grepl("AR", omim_inheritance), 1, 0)) %>% 
  arrange(desc(knownAR), desc(maxpriorityscore), ref_gene_annovar, desc(priority_score)) %>% 
  select(-maxpriorityscore, -knownAR, -recessive_cnt)

#ad are those omim unknown or AD inheritance. 
ad <- gemini_filtered3 %>% filter(!chr_annovar %in% c("X", "Y"), 
                                  recessive_cnt == 1 & !omim_inheritance %in% c("AR") | recessive_cnt >=2 & grepl("AD", omim_inheritance),
                                  pmaxaf < 0.002, priority_score >= 5) %>% 
  arrange(desc(maxpriorityscore), ref_gene_annovar, desc(priority_score)) %>% 
  select(-maxpriorityscore, -recessive_cnt)
#grepl("AD", omim_inheritance) | is.na(omim_inheritance)

#AD: score > 4 , AR: score > 4, all: score > 3

acmg_genes = c('ACTA2','ACTC1','APC','APOB','ATP7B','BMPR1A','BRCA1','BRCA2',
               'CACNA1S','COL3A1','DSC2','DSG2','DSP','FBN1','GLA','KCNH2','KCNQ1',
               'LDLR','LMNA','MEN1','MLH1','MSH2','MSH6','MUTYH','MYBPC3','MYH11',
               'MYH7','MYL2','MYL3','NF2','OTC','PCSK9','PKP2','PMS2','PRKAG2',
               'PTEN','RB1','RET','RYR1','RYR2','SCN5A','SDHAF2','SDHB','SDHC',
               'SDHD','SMAD3','SMAD4','STK11','TGFBR1','TGFBR2','TMEM43','TNNI3',
               'TNNT2','TP53','TPM1','TSC1','TSC2','VHL','WT1')
acmg <- gemini_filtered3 %>% filter(gene %in% acmg_genes, priority_score > 4)
openxlsx::write.xlsx(list("AD" = ad, "AR" = ar, "XR" = xR, "XD" = xD, "ACMG59" = acmg, "all" = gemini_filtered2), file = args[6])

