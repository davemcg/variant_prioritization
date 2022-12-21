#updated 7/28/19
## dbscSNV's ada > 0.8 and rf>0.5, Score=+3 (equals weight of PM)
## |dpsi_max_tissue+dpsi_zscore| > 6, score=+3; > 3, score=+1 make a histogram of dpsi of WGS set and determine the cut-off
## genesplicer (H|M) or maxentscan_diff > 3, score=+3
## spliceai_rank >0.8, score=+8; >0.5, score=+6; >0.2, score=+3; >0.15, score=+1; calculated in snakemake file. make a histogram of splicai score and determine the cut-off
## splice_score = min(8, spliceai etc)
## if PVS == 1 or maxaf > 0.02, then splice score is not added to the  priority score.
## other_pred_score is not added to priority score if maxaf > 0.02
## if impact == "missense_variant" & mis_z >= 3.09 & sigmaaf_missense_0001 < 0.005 & pmaxaf < 0.0005, priority socre += 2

#When testing, switch comment lines below.

args <- commandArgs(trailingOnly=TRUE)

#args <- c("gemini_tsv/D1280_01_BP70015.BP.blueprint_pedigree.gemini.tsv",
#           "gemini_tsv/D1280_01_BP70015.BP.blueprint_pedigree.gemini.ref.tsv", "Z:/resources/OGLpanelGeneDxORcandidate.xlsx", "rearranged.tsv", "filtered.tsv", "108976P", "filtered.xlsx", "0.5", "W:/ddl_nisc_custom_capture/042020/CoNVaDING/CNV_hiSens/108976P.b37.aligned.only.best.score.shortlist.txt")

gemini_file <- args[1]
gemini_ref_var_file <- args[2]
geneCategory_file <- args[3]
rearrangedGemini_file <- args[4]
filteredGemini_tsv_file <- args[5]
sampleName <- args[6]
gemini_xlsx_file <- args[7]
aafCutoff <- args[8]
manta_file <- args[9]
manta_freq_file <- args[10]
roh_file <- args[11]
scramble_mei_file <- args[12]
scramble_del_file <- args[13]
config_file <- args[14]
convading_file <- args[15]
convading_LAF <- args[16]

library(tidyverse)
library(readxl)
library(RColorBrewer)

gemini_input <- read_tsv(gemini_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  mutate(atac_rpe_score = gsub(",", "_", atac_rpe_score)) %>% 
  type_convert() %>% mutate(exon = sub("^", " ", exon), intron = sub("^", " ", intron))
print("###gemini tsv loaded### 10%")
gemini <-  gemini_input %>% mutate( start_vcf = start + 1 ) %>% 
  unite("chr_variant_id", chrom, start_vcf, ref, alt, sep = "-", remove = FALSE ) %>% 
  mutate(sample = sampleName) %>%
  mutate(ref_gene = ifelse(is.na(ref_gene), gene, ref_gene)) 
rm(gemini_input)

#mutate(temp_genes_bed = pmap_chr(list(eyeintegration_gene, gene_gnomad, omim_gene, gene, gene_refgenewithver), ~toString(unique(na.omit(c(...)))) )) %>%
#  mutate(temp_genes_bed = na_if(temp_genes_bed, "") ) %>% 
#above removed so that it doesn't interfere with scoring system. May try in cohort analyses. 10/1/2022

#  mutate(temp_gene = ifelse(grepl(",", gene_refgenewithver), gene, gene_refgenewithver)) 
#InterVar seperate multiple genes for a variant to mulitple lines, then InterVar.R picks the gene with higher priority score. thus this might be safer
#use InterVar gene annotation should be fine, thus remove this line.  ref_gene is from intervar

#http://web.corral.tacc.utexas.edu/WGSAdownload/resources/dbNSFP/dbNSFP4.0b2c.readme.txt
#CADD: https://cadd.gs.washington.edu/info
#eigen:
#http://mutationassessor.org/r3/howitworks.php

#If all_AF<0.002, add 2 to priority score, if < 0.005, add 1 ???

#Read in OGLv1_panel gene class of either Dx or Candidate

#sampleData <- read.tsv(args[1], header = TRUE, check.names=FALSE ) #nrows = 5
# classes <- sapply(gemini, class)
# largeData <- read.csv("huge-file.csv", header = TRUE, colClasses = classes)


#panelGene <- read.delim(args[2], sep = "\t", header = T, colClasses = c("character","character","character") ) %>% 
#  select('gene', 'panel_class') %>% rename(ref_gene = gene)
panelGene <- read_xlsx(geneCategory_file, sheet = "analysis", na = c("NA", "", "None", "NONE", ".")) %>%
  #mutate(ref_gene = toupper(gene)) %>%
  rename(ref_gene = gene) %>% 
  select(ref_gene, panel_class) %>% distinct()
blacklistGene <- read_xlsx(geneCategory_file, sheet = "IVA", na = c("NA", "", "None", "NONE", "."))  %>% filter(Blacklist == "Excluded") %>% pull(Gene)

# get max_priority_score for each gene
max_priority_score <- select(gemini, c(ref_gene, priority_score)) %>% group_by(ref_gene) %>% summarize(maxpriorityscore = max(priority_score)) 
print("###max priority score### 20%")
#arrange by max_priority_score, then by gene, and priority score. None gene region?
gemini_max_priority_score <- left_join(gemini, max_priority_score, by=c("ref_gene"))
rm(gemini)
#VEP hg19 version's gene names are the same as in the IDT ordering design sheets. This is what used for left_join
gemini_rearrangeCol <- left_join(gemini_max_priority_score, panelGene, by = c("ref_gene")) %>% 
  mutate(note = "") %>% separate(vcf_id, c('caller', 'hg38_id'), sep = "_") %>% 
  mutate(hg38_pos = sub("[ACGT]*>[ACGT]*", "", hg38_id)) %>% 
  mutate(gno2e3g_hom = ifelse(is.na(gno2x_hom) & is.na(gno3_nhomalt), NA, ifelse(is.na(gno2x_hom), 0, gno2x_hom) + ifelse(is.na(gno3_nhomalt), 0, gno3_nhomalt) ), 
         gno2e3g_ac = ifelse(is.na(gno2x_ac_all) & is.na(gno3_ac_all), NA, ifelse(is.na(gno2x_ac_all), 0, gno2x_ac_all) + ifelse(is.na(gno3_ac_all), 0, gno3_ac_all) ), 
         gno2e3g_an = ifelse(is.na(gno2x_an_all) & is.na(gno3_an_all), NA, ifelse(is.na(gno2x_an_all), 0, gno2x_an_all) + ifelse(is.na(gno3_an_all), 0, gno3_an_all) )) %>% 
  mutate(gno2e3g_af = gno2e3g_ac/gno2e3g_an) %>% 
  unite("gno2e3g_acan", gno2e3g_ac, gno2e3g_an, sep = "/", remove = TRUE) %>% 
  mutate(gno2x_expected_an = case_when(chrom %in% c("X", "chrX") & gno2x_nonpar == "1" ~ 183653,
                                       chrom %in% c("Y", "chrY") & gno2x_nonpar == "1" ~ 67843,
                                       TRUE ~ 251496)) %>%
  mutate(gno3_expected_an = case_when(chrom %in% c("X", "chrX") & gno3_nonpar == "1" ~ 116830,
                                      chrom %in% c("Y", "chrY") & gno3_nonpar == "1" ~ 35482,
                                      TRUE ~ 152312)) %>%
  mutate(gno2x_filter = ifelse(gno2x_an_all > 0 & is.na(gno2x_filter) & gno2x_an_all < gno2x_expected_an/2, "Less50", gno2x_filter),
         gno3_filter = ifelse(gno3_an_all > 0 & is.na(gno3_filter) & gno3_an_all < gno3_expected_an/2, "Less50", gno3_filter) ) %>%
  mutate(eyeGene = case_when(panel_class == "Dx" ~ 2,
                             panel_class == "Candidate" ~ 1,
                             TRUE ~ 0)) %>% 
  select(-gno2x_expected_an, -gno3_expected_an) %>% 
  select('ref_gene','sample', 'chr_variant_id','grch37variant_id','chrom', 'start_vcf', 'qual', 'filter', starts_with('gts'), starts_with('gt_'), 'aaf', 'caller','hg38_pos',
         'panel_class', 'priority_score', 'prscore_intervar', 'clinvar_hgmd_score', 'splice_score', 'insilico_score', 'gno2e3g_af', 'gno2e3g_acan', 'pmaxaf',gno2x_af_all,gno2x_filter,gno3_af_all,gno3_filter,'max_af', 'max_af_pops', 'gno2e3g_hom', 'note', func_refgenewithver, exonicfunc_refgenewithver, 
         'refgenewithver', 'gene', mane_select, 'hgvsc', 'hgvsp', 'exon', 'intron', 'aa_length', 'omim_gene', 'omim_inheritance', 'omim_phen', 'pvs1', 'truncating_vep', 'hgmd_id', 'hgmd_class', 'hgmd_phen', hgmd_overlap4aa, 'existing_variation',clnid, clnalleleid,clnsig,'clin_sig',clnsigconf, clnrevstat, clndn, clndisdb, 
         af_oglx, ac_oglx, ac_hom_oglx, an_oglx, af_oglg, ac_oglg, ac_hom_oglg, an_oglg, intervar_and_evidence, 'interpro_domain', 'pfam_domain', 'rmsk', 'existing_inframe_oorfs','existing_outofframe_oorfs','existing_uorfs','five_prime_utr_variant_annotation','five_prime_utr_variant_consequence',
         'spliceai', 'spliceai_maxscore', 'spliceaimasked50', 'spliceaimasked50max', 'squirls_interpretation', 'squirls_maxscore', 'squirls_score', 'dbscsnv_ada_score', 'dbscsnv_rf_score', 'regsnp_fpr','regsnp_disease','regsnp_splicing_site','dpsi_max_tissue', 'dpsi_zscore', 'genesplicer', 'maxentscan_diff', 'branchpoint_prob', 'regsnp_fpr','regsnp_disease','regsnp_splicing_site',  
         'sift_pred', 'polyphen_pred', 'mutscore', 'mutationassessor_pred', 'mutationtaster_pred', 'metasvm_pred','metasvm_score', 'clinpred_score', 'primateai_rankscore', 'revel_score', hmc_score, 'ccr_pct','mpc_score', 'mtr_score', 'mtr_fdr', 'mtr_pct', 'cadd_raw', 'cadd_phred','remm', 'fathmm_xf_coding_score','fathmm_xf_noncoding','eigen_pc_raw_coding', 'gerpplus_rs', 'phylop100way_vertebrate', 
         'atac_rpe_score','atac_rpe_itemrgb', 'ft_ret_rpe_score', cherry_sum_score, 'gene_refgenewithver', 'avsnp150', 'tfbs', 'pli','pnull', 'prec', 'mis_z', 'sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01', 'atac_rpe_itemrgb', 'atac_rpe_score', 
         'eyeintegration_rpe_adulttissue', 'eyeintegration_rpe_cellline', 'eyeintegration_rpe_fetaltissue', 'eyeintegration_rpe_stemcellline', 'eyeintegration_retina_adulttissue', 'eyeintegration_retina_stemcellline', 'eyeintegration_wholeblood', 
         'pubmed','sift_score', 'polyphen_score', 'metasvm_rankscore', 'metasvm_score', 'provean_score', 'provean_converted_rankscore',	'provean_pred',
         'f1000g2015aug_all','esp6500siv2_all',gno2_xg_ratio:gno3_popmax, 'syn_z', everything() ) 
#4/12/20: removed 'chr_annovar', 'start_annovar', 'ref_annovar', 'alt_annovar',
print("###gemini_rearranged###")
write_tsv(gemini_rearrangeCol, file.path('.', rearrangedGemini_file), na="")
print("###rearranged file written### 30%")
# mastermind_counts       mastermind_mmid3 are in the VEP output
#nrow(df) == 0 or dim(df)[1] == 0
gemini_ref_var_input <- read_tsv(gemini_ref_var_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) 

if (nrow(gemini_ref_var_input) == 0) {
  gemini_ref_var_rearrangeCol <- data.frame("sample" = sampleName, "note" = "Empty rare reference allele search")
} else {
  gemini_ref_var_input <- gemini_ref_var_input %>%
    mutate(atac_rpe_score = sub(",", "_", atac_rpe_score)) %>% 
    type_convert() %>% mutate(exon = sub("^", " ", exon), intron = sub("^", " ", intron)) %>% 
    mutate( start_vcf = start + 1 ) %>% 
    unite("chr_variant_id", chrom, start_vcf, ref, alt, sep = "-", remove = FALSE ) %>%
    mutate(sample = sampleName) %>%
    mutate(ref_gene = ifelse(is.na(ref_gene), gene, ref_gene)) #%>% 
    #mutate(temp_gene = toupper(ref_gene)) %>% 
    #select(-temp_genes_bed) 
  
  gemini_ref_var_rearrangeCol <- left_join(gemini_ref_var_input, panelGene, by = c("ref_gene")) %>% 
    mutate(note = "") %>% separate(vcf_id, c('caller', 'hg38_id'), sep = "_") %>% 
    mutate(hg38_pos = sub("[ACGT]*>[ACGT]*", "", hg38_id)) %>% 
    mutate(gno2e3g_hom = ifelse(is.na(gno2x_hom) & is.na(gno3_nhomalt), NA, ifelse(is.na(gno2x_hom), 0, gno2x_hom) + ifelse(is.na(gno3_nhomalt), 0, gno3_nhomalt) ), 
           gno2e3g_ac = ifelse(is.na(gno2x_ac_all) & is.na(gno3_ac_all), NA, ifelse(is.na(gno2x_ac_all), 0, gno2x_ac_all) + ifelse(is.na(gno3_ac_all), 0, gno3_ac_all) ), 
           gno2e3g_an = ifelse(is.na(gno2x_an_all) & is.na(gno3_an_all), NA, ifelse(is.na(gno2x_an_all), 0, gno2x_an_all) + ifelse(is.na(gno3_an_all), 0, gno3_an_all) )) %>% 
    mutate(gno2e3g_af = gno2e3g_ac/gno2e3g_an) %>% 
    filter(!ref_gene %in% blacklistGene, gno2e3g_af > 0.98) %>% 
    unite("gno2e3g_acan", gno2e3g_ac, gno2e3g_an, sep = "/", remove = TRUE) %>% 
    mutate(gno2x_expected_an = case_when(chrom %in% c("X", "chrX") & gno2x_nonpar == "1" ~ 183653,
                                         chrom %in% c("Y", "chrY") & gno2x_nonpar == "1" ~ 67843,
                                         TRUE ~ 251496)) %>%
    mutate(gno3_expected_an = case_when(chrom %in% c("X", "chrX") & gno3_nonpar == "1" ~ 116830,
                                        chrom %in% c("Y", "chrY") & gno3_nonpar == "1" ~ 35482,
                                        TRUE ~ 152312)) %>%
    mutate(gno2x_filter = ifelse(gno2x_an_all > 0 & is.na(gno2x_filter) & gno2x_an_all < gno2x_expected_an/2, "Less50", gno2x_filter),
           gno3_filter = ifelse(gno3_an_all > 0 & is.na(gno3_filter) & gno3_an_all < gno3_expected_an/2, "Less50", gno3_filter) ) %>%
    mutate(eyeGene = case_when(panel_class == "Dx" ~ 2,
                               panel_class == "Candidate" ~ 1,
                               TRUE ~ 0)) %>% 
    select(-gno2x_expected_an, -gno3_expected_an) %>% 
    select('ref_gene','sample', 'chr_variant_id','grch37variant_id','chrom', 'start_vcf', 'qual', 'filter', starts_with('gts'), starts_with('gt_'), 'aaf', 'caller','hg38_pos',
           'panel_class', 'priority_score', 'prscore_intervar', 'clinvar_hgmd_score', 'splice_score', 'insilico_score', 'gno2e3g_af', 'gno2e3g_acan', 'pmaxaf',gno2x_af_all,gno2x_filter,gno3_af_all,gno3_filter,'max_af', 'max_af_pops', 'gno2e3g_hom', 'note', func_refgenewithver, exonicfunc_refgenewithver, 
           'refgenewithver', 'gene', mane_select, 'hgvsc', 'hgvsp', 'exon', 'intron', 'aa_length', 'omim_gene', 'omim_inheritance', 'omim_phen', 'pvs1', 'truncating_vep', 'hgmd_id', 'hgmd_class', 'hgmd_phen', hgmd_overlap4aa, 'existing_variation', clnalleleid,clnsig,'clin_sig', clnrevstat, clndn, clndisdb, 
           af_oglx, ac_oglx, ac_hom_oglx, an_oglx, af_oglg, ac_oglg, ac_hom_oglg, an_oglg, intervar_and_evidence, 'interpro_domain', 'pfam_domain', 'rmsk', 'existing_inframe_oorfs','existing_outofframe_oorfs','existing_uorfs','five_prime_utr_variant_annotation','five_prime_utr_variant_consequence',
           'spliceai', 'spliceai_maxscore', 'spliceaimasked50', 'spliceaimasked50max', 'squirls_interpretation', 'squirls_maxscore', 'squirls_score', 'dbscsnv_ada_score', 'dbscsnv_rf_score', 'regsnp_fpr','regsnp_disease','regsnp_splicing_site','dpsi_max_tissue', 'dpsi_zscore', 'genesplicer', 'maxentscan_diff', 'branchpoint_prob', 'regsnp_fpr','regsnp_disease','regsnp_splicing_site',  
           'sift_pred', 'polyphen_pred', 'mutscore', 'mutationassessor_pred', 'mutationtaster_pred', 'metasvm_pred','metasvm_score', 'clinpred_score', 'primateai_rankscore', 'revel_score', hmc_score, 'ccr_pct','mpc_score', 'mtr_score', 'mtr_fdr', 'mtr_pct', 'cadd_raw', 'cadd_phred','remm', 'fathmm_xf_coding_score','fathmm_xf_noncoding','eigen_pc_raw_coding', 'gerpplus_rs', 'phylop100way_vertebrate', 
           'atac_rpe_score','atac_rpe_itemrgb', 'ft_ret_rpe_score', cherry_sum_score, 'gene_refgenewithver', 'avsnp150', 'tfbs', 'pli','pnull', 'prec', 'mis_z', 'sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01', 'atac_rpe_itemrgb', 'atac_rpe_score', 
           'eyeintegration_rpe_adulttissue', 'eyeintegration_rpe_cellline', 'eyeintegration_rpe_fetaltissue', 'eyeintegration_rpe_stemcellline', 'eyeintegration_retina_adulttissue', 'eyeintegration_retina_stemcellline', 'eyeintegration_wholeblood', 
           'pubmed','sift_score', 'polyphen_score', 'metasvm_rankscore', 'metasvm_score', 'provean_score', 'provean_converted_rankscore',	'provean_pred',
           'f1000g2015aug_all','esp6500siv2_all',gno2_xg_ratio:gno3_popmax, 'syn_z', everything()) %>%  
    arrange(desc(eyeGene), ref_gene)
}

# mutate(temp_genes_bed = pmap_chr(list(eyeintegration_gene, gene_gnomad, omim_gene, gene, gene_refgenewithver), ~toString(unique(na.omit(c(...)))) )) %>%
#   mutate(temp_genes_bed = na_if(temp_genes_bed, "") ) %>% 
#   mutate(ref_gene = ifelse(is.na(ref_gene), temp_genes_bed, ref_gene)) %>%  

gemini_filtered <- gemini_rearrangeCol %>% mutate(temp_group = ifelse(priority_score >= 3, 3, ifelse(priority_score >= -3, -3, -4))) %>% # checked OPA1 non-coding regions, the AF for some of variants with score -3 are around 0.01.
  filter(!ref_gene %in% blacklistGene, priority_score >= 15 | (temp_group >= -3 & pmaxaf < 0.1 & aaf < aafCutoff & af_oglg < 0.05 & af_oglx < 0.05) ) %>%
  arrange(desc(eyeGene), desc(temp_group), desc(maxpriorityscore), ref_gene, desc(priority_score)) 

gemini_filtered0 <- gemini_filtered %>% select(-maxpriorityscore) #filtered0: temp_group >= -3 & pmaxaf < 0.2 & aaf < aafCutoff) | priority_score >= 10 

#AnnotSV does not have GD_POPMAX_AF column as of 3/2/2021, if GD_POPMAX_AF in columnames, then use it in the next version.
write_tsv(gemini_filtered0, file.path('.', filteredGemini_tsv_file), na="")

#found out whether blank is 0 after importing, use 1197 for filtering
#consider adding manta to the main gemini df for sorting/filtering after knowing the specificity of the manta calls. To better sort AR, AD, and ACMG 2nd.

manta_freq <- read_tsv(manta_freq_file, col_names = TRUE, na = c("NA", "full=NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
  type_convert()
  
if (file.size(manta_file) == 0) {
  manta_sort <- data.frame("sample" = sampleName, "note" = "Empty manta")
} else {
  manta_original <- read_tsv(manta_file, col_names = TRUE, na = c("NA", "full=NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
    filter(FILTER == 'PASS') %>% 
    mutate(ACMG_class = sub("full=", "", ACMG_class)) %>% 
    type_convert() %>%
    mutate(temp_SV_start = round(SV_start, -3), temp_SV_end = round(SV_end, -3)) %>%
    unite("variant", SV_chrom, temp_SV_start, temp_SV_end, SV_type, sep = "-", remove = FALSE) 
  manta <- left_join(manta_original, manta_freq, by = c("variant") ) %>% 
    filter(is.na(CohortFreq) | CohortFreq < 0.025) %>% 
    filter( Annotation_mode == "split" | Gene_count == 0) %>%
    mutate(ACMG_class = case_when(!is.na(ACMG_class) & SV_type == "DEL" & !is.na(B_loss_source) ~ ACMG_class - 2,
                                  !is.na(ACMG_class) & SV_type == "DUP" & !is.na(B_gain_source) ~ ACMG_class - 1,
                                  SV_type == "BND" & ( ( SV_chrom %in% c("X", "chrX") & between(SV_start, 140420272, 140421278) | ( grepl("chrX|X", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 140420272, 140421278) ) )) ~ 5,
                                  ( SV_chrom %in% c("17", "chr17") & between(SV_start, 59105674, 59683460) ) | ( grepl("chr17|17", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 59105674, 59683460) )  ~ 5,
                                  TRUE ~ ACMG_class ) ) %>% 
    filter(ACMG_class > 1 | is.na(ACMG_class), is.na(SV_length) | abs(SV_length) < 1000000) %>% #added is.na(ACMG_class) 11/17/2021, may need to remove this part if too many lines in the results
    separate(Location, c('temp_location1', 'temp_location2'), sep = "-", remove = FALSE, convert = FALSE) %>% 
    filter(!(grepl("intron", temp_location1) & temp_location1 == temp_location2 & (!is.na(B_gain_source) | !is.na(B_loss_source) | !is.na(B_ins_source)) )) %>% 
    select(-starts_with('temp_'), -variant) %>% 
    mutate(ref_gene = toupper(Gene_name))
  manta_sort <- left_join(manta, panelGene, by = c("ref_gene")) %>% 
    mutate(note = "") %>% 
    mutate(panel_class == ifelse(SV_type == "BND" & ( ( SV_chrom %in% c("X", "chrX") & between(SV_start, 140420272, 140421278) | ( grepl("chrX|X", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 140420272, 140421278) ) )), "Dx", panel_class)) %>% 
    mutate(eyeGene = case_when(panel_class == "Dx" ~ 2,
                               panel_class == "Candidate" ~ 1,
                               SV_type == "BND" & ( ( SV_chrom %in% c("X", "chrX") & between(SV_start, 140420272, 140421278) | ( grepl("chrX|X", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 140420272, 140421278) ) )) ~ 3,
                               ( SV_chrom %in% c("17", "chr17") & between(SV_start, 59105674, 59683460) ) | ( grepl("chr17|17", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 59105674, 59683460) )  ~ 2.5,
                               TRUE ~ 0)) %>% # GRCh38 coordinates
    arrange(desc(eyeGene), desc(ACMG_class)) %>% 
    mutate(note = ifelse(SV_type == "BND" & ( ( SV_chrom %in% c("X", "chrX") & between(SV_start, 140420272, 140421278) | ( grepl("chrX|X", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 140420272, 140421278) ) )), "XLFD", note)) %>% 
    mutate(note = ifelse(( SV_chrom %in% c("17", "chr17") & between(SV_start, 59105674, 59683460) ) | ( grepl("chr17|17", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 59105674, 59683460) ), "RP17", note)) %>% 
    select(AnnotSV_ID:Gene_name,panel_class,ACMG_class,note,CohortFreq,`NaltP/NtotalP`,OMIM_phenotype:GnomAD_pLI,Frameshift,everything() )
  
}

if ( roh_file == "filePlaceholder") {
  roh <- data.frame("sample" = sampleName, "note" = "Roh not analyzed.")
} else if ( file.size(roh_file) == 0) { roh <- data.frame("sample" = sampleName, "note" = "Empty roh") 
} else {
  roh <- read_tsv(roh_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
    type_convert() %>% 
    mutate(autosome = ifelse(`#Chr` %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
                                           "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22" ), "autosome", "X")) %>% 
    group_by(autosome) %>% 
    mutate(pct = sum(`Size(Mb)`)/2900) %>% 
    select(-autosome)
}

if (file.size(scramble_mei_file) == 0) {
  scramble_mei <- data.frame("sample" = sampleName, "note" = "Empty scramble mei")
} else {
  scramble_mei <- read_xlsx(scramble_mei_file, na = c("NA", "", "None", "NONE", "."))
}

if (scramble_del_file == "filePlaceholder") {
  scramble_del_sort <- data.frame("sample" = sampleName, "note" = "Scramble del not performed")
} else if (file.size(scramble_del_file) == 0) {
  scramble_del_sort <- data.frame("sample" = sampleName, "note" = "Empty scramble del")
} else {
  scramble_del <- read_tsv(scramble_del_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
    filter( is.na(Gene_name) | Annotation_mode == 'split' & !is.na(Gene_name) ) %>%
    filter(is.na(CohortFreq) | CohortFreq < 0.025) %>% 
    mutate(ACMG_class = sub("full=", "", ACMG_class)) %>% 
    type_convert() %>% 
    mutate(ACMG_class = case_when(SV_type == "DEL" & !is.na(B_loss_source) ~ ACMG_class - 2,
                                  SV_type == "DUP" & !is.na(B_gain_source) ~ ACMG_class - 1,
                                  TRUE ~ ACMG_class )) %>% 
    filter(ACMG_class > 1, is.na(SV_length) | abs(SV_length) < 1000000) %>% 
    separate(Location, c('temp_location1', 'temp_location2'), sep = "-", remove = FALSE, convert = FALSE) %>% 
    filter(!(grepl("intron", Location) & temp_location1 == temp_location2)) %>% 
    select(-starts_with('temp_')) %>% 
    rename(ref_gene = "Gene_name")
  scramble_del_sort <- left_join(scramble_del, panelGene, by = c("ref_gene")) %>% 
    mutate(note = "") %>% 
    mutate(eyeGene = case_when(panel_class == "Dx" ~ 2,
                               panel_class == "Candidate" ~ 1,
                               TRUE ~ 0)) %>% 
    arrange(desc(eyeGene), desc(ACMG_class)) %>% 
    select(AnnotSV_ID:Annotation_mode,OMIM_phenotype:ACMG_class,panel_class,note,everything() ) 
  if (dim(scramble_del_sort)[1] == 0) {
    scramble_del_sort <- scramble_del_sort %>% add_row(note = "no scramble del candidate after filtering with scramble db")
  } 
}
  
#deleted sampleName, after FORMAT, add this back for production

#manta_file "Z:/NextSeqAnalysis/test2/manta/manta.1197.annotated.tsv"
gemini_filtered1 <- gemini_filtered %>% filter(priority_score >= 3) %>% select(-maxpriorityscore) %>% 
  arrange(desc(eyeGene), desc(priority_score))
  
gemini_filtered2 <- gemini_filtered %>% filter(priority_score >= 4) %>% 
  rename_all(funs(str_replace(., sampleName, ""))) %>% 
  mutate(temp_ref_gene = case_when(ref_gene %in% c("ROM1", "PRPH2") ~ "PRPH2-ROM1",
                                   ref_gene %in% c("PCDH15", "CDH23") ~ "PCDH15-CDH23",
                                   ref_gene %in% c("CNGA3", "CNGB3") ~ "CNGA3-CNGB3",
                                   ref_gene %in% c("CNGA1", "CNGB1") ~ "CNGA1-CNGB1",
                                   TRUE ~ ref_gene)) # digenic recessive

recessive_count <- select(gemini_filtered2, c(temp_ref_gene, priority_score, gt_types.)) %>%
  filter(priority_score >= 5) %>% select(-priority_score) %>% 
  group_by(temp_ref_gene) %>% summarize(recessive_cnt = sum(gt_types.)) 

gemini_filtered3 <- left_join(gemini_filtered2, recessive_count, by=c("temp_ref_gene")) %>% 
  replace_na(list(recessive_cnt=0)) %>% 
  mutate(recessive_cnt = as.integer(recessive_cnt)) %>% 
  select(-temp_ref_gene)

xR <- gemini_filtered3 %>% filter(chrom %in% c("X", "chrX"), priority_score >= 5, recessive_cnt >= 2) %>% select(-maxpriorityscore, -recessive_cnt)
xD <- gemini_filtered3 %>% filter(chrom %in% c("X", "chrX"), recessive_cnt == 1, pmaxaf < 0.002) %>% select(-maxpriorityscore, -recessive_cnt)

#ar are those genes with homozygous or compound hets variants of ps >= 5. However, ps = 4 variants were also listed if there are 2 ps>=5.
#ar gene with 1 hit will not be here.
ar <- gemini_filtered3 %>% filter(!chrom %in% c("X", "Y", "chrX", "chrY"), !omim_inheritance %in% c("AD"), priority_score >= 5, recessive_cnt >= 2) %>% 
  mutate(knownAR = ifelse(grepl("AR", omim_inheritance), 1, 0)) %>% 
  arrange(desc(eyeGene), desc(knownAR), desc(maxpriorityscore), ref_gene, desc(priority_score)) %>% 
  mutate(note = ifelse(ref_gene %in% c("PRPH2", "ROM1", "PCDH15", "CDH23", "CNGA1", "CNGB1", "CNGA3", "CNGB3"), "Digenic?", note) ) %>%
  select(-maxpriorityscore, -knownAR, -recessive_cnt) # digenic recessive

#ad are those omim unknown or AD inheritance. 
ad <- gemini_filtered3 %>% filter(!chrom %in% c("X", "Y", "chrX", "chrY"), 
                                  recessive_cnt == 1 & !omim_inheritance %in% c("AR") | recessive_cnt >=2 & grepl("AD", omim_inheritance),
                                  pmaxaf < 0.002, priority_score >= 5) %>% 
  arrange(desc(eyeGene), desc(maxpriorityscore), ref_gene, desc(priority_score)) %>% 
  select(-maxpriorityscore, -recessive_cnt)
print("###inheritance search done### 50%")
#grepl("AD", omim_inheritance) | is.na(omim_inheritance)

#AD: score > 4 , AR: score > 4, all: score >= 3

acmg_genes <- read_xlsx(geneCategory_file, sheet = "ACMG", na = c("NA", "", "None", "NONE", ".")) %>% pull(Gene) %>% unique()

# acmg_genes = c('ACTA2','ACTC1','APC','APOB','ATP7B','BMPR1A','BRCA1','BRCA2',
#                'CACNA1S','COL3A1','DSC2','DSG2','DSP','FBN1','GLA','KCNH2','KCNQ1',
#                'LDLR','LMNA','MEN1','MLH1','MSH2','MSH6','MUTYH','MYBPC3','MYH11',
#                'MYH7','MYL2','MYL3','NF2','OTC','PALB2','PCSK9','PKP2','PMS2','PRKAG2',
#                'PTEN','RB1','RET','RYR1','RYR2','SCN5A','SDHAF2','SDHB','SDHC',
#                'SDHD','SMAD3','SMAD4','STK11','TGFBR1','TGFBR2','TMEM43','TNNI3',
#                'TNNT2','TP53','TPM1','TSC1','TSC2','VHL','WT1')
acmg <- gemini_filtered3 %>% filter(ref_gene %in% acmg_genes, priority_score > 4) %>% select(-maxpriorityscore, -recessive_cnt)
print("###acmg done### 70%")

config <- read_tsv(config_file, col_names = FALSE, na = c("NA", ""), col_types = cols(.default = col_character())) %>% 
  separate("X1", c("tool", "version", "note"), sep = "\\:|\\#", remove = TRUE)

summaryInfo <- data.frame("sample" = sampleName, "PatientDxPhenotype" = NA, "DxOutcome"= NA, "variant" = NA, "reviewer" = NA, "date" = NA) %>% 
  add_row("sample" = sampleName, "PatientDxPhenotype" = NA, "DxOutcome"= NA, "variant" = NA, "reviewer" = NA, "date" = NA)

if (scramble_del_file == "filePlaceholder") {
  openxlsx::write.xlsx(list("AR" = ar, "AD" = ad, "XR" = xR, "XD" = xD, "ACMG" = acmg, "all" = gemini_filtered1, "rareRef" = gemini_ref_var_rearrangeCol, "manta" = manta_sort, "scramble_mei" = scramble_mei, "roh" = roh, "config" = config, "summary" = summaryInfo), file = gemini_xlsx_file, firstRow = TRUE, firstCol = TRUE)
} else if (is.na(convading_file)) {
  openxlsx::write.xlsx(list("AR" = ar, "AD" = ad, "XR" = xR, "XD" = xD, "ACMG" = acmg, "all" = gemini_filtered1, "rareRef" = gemini_ref_var_rearrangeCol, "manta" = manta_sort, "scramble_mei" = scramble_mei, "scramble_del" = scramble_del_sort, "roh" = roh, "config" = config, "summary" = summaryInfo), file = gemini_xlsx_file, firstRow = TRUE, firstCol = TRUE)
} else {
    cnv <- read_tsv(convading_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
      type_convert() 
    if (dim(cnv)[1] == 0) {
      openxlsx::write.xlsx(list("AR" = ar, "AD" = ad, "XR" = xR, "XD" = xD, "ACMG" = acmg, "all" = gemini_filtered1, "rareRef" = gemini_ref_var_rearrangeCol, "manta" = manta_sort, "scramble_mei" = scramble_mei, "scramble_del" = scramble_del_sort, "config" = config, "summary" = summaryInfo), file = gemini_xlsx_file, firstRow = TRUE, firstCol = TRUE)
    } else {
      cnv_gene <- as.list(distinct(cnv, GENE)[[1]])
      #cnv_gene <- dplyr::pull(cnv, GENE) #pull column as a vector
      cnv_variant <- gemini_rearrangeCol %>% 
        rename_all(funs(str_replace(., sampleName, ""))) %>% 
        filter(ref_gene %in% cnv_gene) %>% 
        mutate(LAF = ifelse(gt_alt_freqs. > 0.5, 1 - gt_alt_freqs., gt_alt_freqs.)) %>% 
        select(-gt_alt_freqs.) %>% 
        select('chr_variant_id', 'chrom', 'start', 'qual', 'filter', starts_with('gts'), starts_with('gt_'), 'LAF', 'ref_gene', 'exon', 'ref_gene',  
              'refgenewithver', 'exonicfunc_refgenewithver', 'hgvsc', 'hgvsp', 'type') %>% 
        rename(gene = ref_gene)
      openxlsx::write.xlsx(list("AR" = ar, "AD" = ad, "XR" = xR, "XD" = xD, "ACMG" = acmg, "all" = gemini_filtered1, "rareRef" = gemini_ref_var_rearrangeCol, "CoNVaDING" = cnv, "CNV-variant" = cnv_variant, "manta" = manta_sort, "scramble_mei" = scramble_mei, "scramble_del" = scramble_del_sort, "config" = config, "summary" = summaryInfo), file = gemini_xlsx_file, firstRow = TRUE, firstCol = TRUE)
      cnv_edit <- cnv %>% 
        mutate(START = START - 100, STOP = STOP + 100, type = "snp") %>% #padding of 100 nt
        gather(START:STOP, key = "datatype", value = "position") %>% 
        mutate(LAF = 0.54, gt_depths. = ifelse(datatype == "START", 40, 200) ) %>% 
        rename(chrom = "CHR", gene = "GENE") %>% 
        select(chrom, gene, datatype, position, LAF, gt_depths., type)
      
      cnv_variant_edit <- cnv_variant %>% 
        select(chrom, start, gt_depths.,	LAF, gene, type) %>% 
        rename(position = "start") %>% 
        mutate(datatype="LAF") %>% 
        filter(gt_depths. >= 15)
      
      variantForPlot <- rbind(cnv_edit, cnv_variant_edit) %>% mutate(gene = as.factor(gene)) %>% 
        arrange(chrom, position) %>% 
        mutate(variant_no = row_number()) %>% 
        mutate(DepthGroup = case_when(gt_depths. >= 100 ~ "DP>=100",
                                      gt_depths. >= 30 ~ "DP30-99",
                                      TRUE ~ "DP15-29")) %>% 
        mutate(type = factor(type, levels = c("snp", "indel")))
      
      variantForPlot$DepthGroup = factor(variantForPlot$DepthGroup,levels = c("DP>=100", "DP30-99", "DP15-29"))
      
      print(length(cnv_gene))
      if (length(cnv_gene) <= 6) {
        plot_pdf <- ggplot(variantForPlot, aes(x= variant_no, y = LAF, color = DepthGroup, shape = type)) + 
          scale_color_brewer(palette = "Set1") +
          coord_cartesian(ylim = c(0, 0.55)) +
          labs(title = sampleName, x= 'Variants', y = 'Lesser allele frequency') +
          facet_wrap(~ gene, ncol = 1, scales = "free_x") +
          geom_point(size = 0.5) +
          geom_segment(data = subset(variantForPlot, datatype %in% c('START', 'STOP')), aes(x= variant_no, xend=variant_no, y=0, yend=LAF), linetype="dotted") +
          theme_bw() +
          scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
          theme(axis.text.x  = element_blank(), axis.text.y  = element_text(size=16)) +
          theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16)) +
          theme(legend.position = "right") 
        ggsave(convading_LAF, plot = plot_pdf, width = 16, height = 8 * length(cnv_gene), units = "cm")
      } else if (length(cnv_gene) <= 12) {
        plot_pdf <- ggplot(variantForPlot, aes(x= variant_no, y = LAF, color = DepthGroup, shape = type)) + 
          scale_color_brewer(palette = "Set1") +
          coord_cartesian(ylim = c(0, 0.55)) +
          labs(title = sampleName, x= 'Variants', y = 'Lesser allele frequency') +
          facet_wrap(~ gene, ncol = 2, scales = "free_x") +
          geom_point(size = 0.8) +
          geom_segment(data = subset(variantForPlot, datatype %in% c('START', 'STOP')), aes(x= variant_no, xend=variant_no, y=0, yend=LAF), linetype="dotted") +
          theme_bw() +
          scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
          theme(axis.text.x  = element_blank(), axis.text.y  = element_text(size=8)) +
          theme(axis.title.x = element_text(size=8), axis.title.y = element_text(size=8)) +
          theme(legend.title = element_blank(), legend.position = 'none') 
        ggsave(convading_LAF, plot = plot_pdf, width = 32, height = 4 * length(cnv_gene), units = "cm")
      } else {
        plot_pdf <- ggplot(variantForPlot, aes(x= variant_no, y = LAF, color = DepthGroup, shape = type)) + 
          scale_color_brewer(palette = "Set1") +
          coord_cartesian(ylim = c(0, 0.55)) +
          labs(title = sampleName, x= 'Variants', y = 'Lesser allele frequency') +
          facet_wrap(~ gene, ncol = 6, scales = "free_x") +
          geom_point(size = 0.8) +
          geom_segment(data = subset(variantForPlot, datatype %in% c('START', 'STOP')), aes(x= variant_no, xend=variant_no, y=0, yend=LAF), linetype="dotted") +
          theme_bw() +
          scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
         # theme(axis.text.x  = element_text(size=8), axis.text.y  = element_text(size=8)) +
         # theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16)) +
          theme(legend.title = element_blank(), legend.position = 'none') 
        ggsave(convading_LAF, plot = plot_pdf, width = 48, height = 4/6 * length(cnv_gene), units = "cm", limitsize = FALSE)
      } 
       
    }
}
print("### 100%")

