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
## needs to use if (filePlaceholder) condition and/or file.size condition so that OGL NGS_calling/or vcf input only are used.

args <- commandArgs(trailingOnly=TRUE)

# setwd("Z:/exome/x/prioritization")
# args <- c("gemini_tsv/J01x3.YH1.gt3.J01.gemini.tsv.gz",
#           "gemini_tsv/J01x3.YH1.gt3.J01.gemini.ref.tsv.gz",
#           "Z:/resources/OGLpanelGeneDxORcandidate.xlsx",
#           "gemini_tsv/J01x3.YH1.gt3.J01.gemini.rearranged.tsv.gz",
#           "gemini_tsv_filtered/J01x3.YH1.gt3.J01.gemini.filtered.tsv.gz",
#           "J01x3",
#           "gemini_xlsx/J01x3.YH1.gt3.J01.gemini.filtered.xlsx",
#           "1.1",
#           "../manta/manta.J01x3.annotated.tsv",
#           "Z:/resources/manta/manta.OGL.freq.2022-09.tsv",
#           "filePlaceholder",
#           "Z:/resources/jaxCNV/jaxCNV.OGL.freq.2025-01.tsv",
#           "../AutoMap/J01x3/J01x3.HomRegions.annot.tsv",
#           "../scramble_anno/J01x3.scramble.xlsx",
#           "../scramble_anno/J01x3.scramble.del.tsv",
#           "config_variant_prioritization_J01.yaml",
#           "filePlaceholder")

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
jaxcnv_file <- args[11]
jaxcnv_freq_file <- args[12]
roh_file <- args[13]
scramble_mei_file <- args[14]
scramble_del_file <- args[15]
config_file <- args[16]
clinsv_file <- args[17]
mutserve_file <- args[18]
convading_file <- args[19]
convading_LAF <- args[20]

library(tidyverse)
library(readxl)
library(RColorBrewer)

gemini_input <- read_tsv(gemini_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  mutate(atac_rpe_score = gsub(",", "_", atac_rpe_score)) %>% 
  type_convert() %>% mutate(exon = sub("^", " ", exon), intron = sub("^", " ", intron)) %>% 
  rename_with(., ~ sub(paste0("\\.", sampleName, "$"), "", .x))
print("###gemini tsv loaded### 10%")
gemini <-  gemini_input %>% mutate( start_vcf = start + 1 ) %>% 
  unite("chr_variant_id", chrom, start_vcf, ref, alt, sep = "-", remove = FALSE ) %>% 
  mutate(sample = sampleName) %>%
  mutate(ref_gene = ifelse(is.na(ref_gene), gene, ref_gene)) %>% 
  replace_na(list(max_af = -1)) 
  
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
  select(ref_gene, panel_class, GenePhenotypeCategory, Phenotypes) %>% distinct()
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
  mutate(gno2e3g_hemi = ifelse(gno2x_nonpar + gno3_nonpar == 0, NA, ifelse(gno2x_nonpar == 0, 0, gno2x_ac_xy) + ifelse(gno3_nonpar == 0, 0, round(gno3_an_xy * gno3_af_xy, 0) ) )) %>%
  mutate(gno2e3g_hom = ifelse(is.na(gno2x_hom) & is.na(gno3_nhomalt), NA, ifelse(is.na(gno2x_hom), 0, gno2x_hom) + ifelse(is.na(gno3_nhomalt), 0, gno3_nhomalt) ), 
         gno2e3g_ac = ifelse(is.na(gno2x_ac_all) & is.na(gno3_ac_all), NA, ifelse(is.na(gno2x_ac_all), 0, gno2x_ac_all) + ifelse(is.na(gno3_ac_all), 0, gno3_ac_all) ), 
         gno2e3g_an = ifelse(is.na(gno2x_an_all) & is.na(gno3_an_all), NA, ifelse(is.na(gno2x_an_all), 0, gno2x_an_all) + ifelse(is.na(gno3_an_all), 0, gno3_an_all) )) %>% 
  mutate(gno2e3g_af = gno2e3g_ac/gno2e3g_an) %>% 
  replace_na(list(gno2e3g_af = -1 )) %>% 
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
                             panel_class == "Candidate-High" ~ 1.6,
                             panel_class == "Candidate-refutedDxGene" ~ 1.4,
                             panel_class == "Candidate" ~ 1,
                             panel_class == "Dx-modifier-rare" ~ 0.8,
                             panel_class == "Dx-modifier-common" ~ 0.6,
                             TRUE ~ 0)) %>% 
  unite("temp_panel_class", panel_class, GenePhenotypeCategory, sep = "|", remove = TRUE, na.rm = TRUE) %>%
  rename(panel_class = temp_panel_class) %>%
  unite("temp_omim_phen", omim_phen, Phenotypes, sep = "::", remove = TRUE, na.rm = TRUE) %>%
  rename(omim_phen = temp_omim_phen) %>%
  select(-gno2x_expected_an, -gno3_expected_an) %>% 
  mutate(gt_alt_freqs = round(gt_alt_freqs,2),
         aaf = round(aaf, 3),
         gno2x_af_all = round(gno2x_af_all, 6),
         gno3_af_all = round(gno3_af_all, 6),
         max_af = round(max_af, 6),
         af_oglx = round(af_oglx, 5),
         af_oglg = round(af_oglg, 5),
         pli = round(pli, 3),
         loeuf = round(loeuf, 2),
         mis_z = round(mis_z, 2),
         gno2e3g_af = round(gno2e3g_af,6),
         gnomad_nc_constraint = round(gnomad_nc_constraint,2) ) %>%
  #unite("hgmd", hgmd_class, hgmd_id, na.rm = TRUE, remove = FALSE) %>%
  mutate(
    temp_hgvsc = if_else(is.na(hgvsc), NA_character_, str_replace(hgvsc, "^[^:]*:", "")),
    temp_prot  = if_else(is.na(hgvsp), NA_character_, str_replace(hgvsp, "^[^:]*:p.", "")),
    temp_prot = sub("%3D", "=", temp_prot),
    hgvs_trscrpt = str_extract(hgvsc, "^[^:]+"),
    vep_gene = if_else(is.na(hgvs_trscrpt), gene, str_c(gene, " ", hgvs_trscrpt, sep = "")),
    vep_hgvs = if_else(is.na(temp_prot), temp_hgvsc, str_c(temp_hgvsc, " p.(", temp_prot, ")", sep = "")),
  ) %>%
  mutate(gnomad_af = if_else(is.na(gno2x_af_all) | gno2x_af_all == -1, gno3_af_all, gno2x_af_all)) %>% 
  select(-temp_hgvsc, -temp_prot, -hgvs_trscrpt) %>%
  #mutate(zygosity = ifelse(gt_types == 1, "Heterozygous", ifelse(chrom %in% c("chrX", "X"), "Hemizygous", "Homozygous"))) %>% 
  select('ref_gene','sample','chr_variant_id',gts,'hgmd_class',clnsig,clnsigconf,'refgenewithver',vep_gene,vep_hgvs,gt_types,hgmd_id,clnid,gnomad_af,'gno2e3g_acan', 'gno2e3g_hom',gno2e3g_hemi,
         'panel_class','caller', 'priority_score', 'prscore_intervar', 'clinvar_hgmd_score', 'splice_score', 'insilico_score', 
         'omim_inheritance','omim_phen','hg38_pos', 'note','aaf','oe_lof_upper_bin','pli','loeuf','mis_z',exonicfunc_refgenewithver, 'exon','aa_length','intron',
         'grch37variant_id',gt_alt_freqs,gt_depths,gt_quals, 
         'max_af', 'max_af_pops',pmaxaf, af_oglx,af_oglg,gno2x_af_all,gno3_af_all,'interpro_domain', 'pfam_domain','pnull','prec', 'rmsk', 
         promoterai, 'am_pathogenicity','am_class','spliceai','spliceai_maxscore','spliceaimasked50','spliceaimasked50max',
         'qual',gno2x_filter,gno3_filter,func_refgenewithver,'gene',mane_select,'hgvsc','hgvsp','omim_gene', 'hgmd_phen', hgmd_overlap4aa, 'existing_variation', clnalleleid,'clin_sig', clnrevstat, clndn, clndisdb, 
         intervar_and_evidence, 'pvs1', 'truncating_vep', 'gno2e3g_af', ac_oglx, ac_hom_oglx, an_oglx, ac_oglg, ac_hom_oglg, an_oglg, 'existing_inframe_oorfs','existing_outofframe_oorfs','existing_uorfs','five_prime_utr_variant_annotation','five_prime_utr_variant_consequence',
         'squirls_interpretation', 'squirls_maxscore', 'squirls_score', 'dbscsnv_ada_score', 'dbscsnv_rf_score', 'dpsi_max_tissue', 'dpsi_zscore', 'genesplicer', 'maxentscan_diff', 'branchpoint_prob', 'labranchor_score', 'regsnp_fpr','regsnp_disease','regsnp_splicing_site',  
         'sift_pred', 'polyphen_pred', 'mutscore', 'mutationassessor_pred', 'mutationtaster_pred', 'metasvm_pred','clinpred_score', 'primateai_rankscore', 'revel_score', hmc_score, 'ccr_pct','mpc_score', 'mtr_score', 'mtr_fdr', 'mtr_pct', 'cadd_raw', 'cadd_phred','remm', 'fathmm_xf_coding_score','fathmm_xf_noncoding','eigen_pc_raw_coding', 'gerpplus_rs', 'phylop100way_vertebrate', 
         gnomad_nc_constraint, 'atac_rpe_score','atac_rpe_itemrgb', 'ft_ret_rpe_score', cherry_sum_score, 'gene_refgenewithver', 'avsnp150', 'tfbs',  'sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01',  
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
    rename_with(., ~ sub(paste0("\\.", sampleName, "$"), "", .x)) %>% 
    mutate( start_vcf = start + 1 ) %>% 
    unite("chr_variant_id", chrom, start_vcf, ref, alt, sep = "-", remove = FALSE ) %>%
    mutate(sample = sampleName) %>%
    mutate(ref_gene = ifelse(is.na(ref_gene), gene, ref_gene)) #%>% 
    #mutate(temp_gene = toupper(ref_gene)) %>% 
    #select(-temp_genes_bed) 
  
  gemini_ref_var_rearrangeCol <- left_join(gemini_ref_var_input, panelGene, by = c("ref_gene")) %>% 
    mutate(note = "") %>% separate(vcf_id, c('caller', 'hg38_id'), sep = "_") %>% 
    mutate(hg38_pos = sub("[ACGT]*>[ACGT]*", "", hg38_id)) %>%
    mutate(gno2e3g_hemi = ifelse(gno2x_nonpar + gno3_nonpar == 0, NA, ifelse(gno2x_nonpar == 0, 0, gno2x_ac_xy) + ifelse(gno3_nonpar == 0, 0, round(gno3_an_xy * gno3_af_xy, 0) ) )) %>% 
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
                               panel_class == "Candidate-High" ~ 1.6,
                               panel_class == "Candidate-refutedDxGene" ~ 1.4,
                               panel_class == "Candidate" ~ 1,
                               panel_class == "Dx-modifier-rare" ~ 0.8,
                               panel_class == "Dx-modifier-common" ~ 0.6,
                               TRUE ~ 0)) %>%
    unite("temp_panel_class", panel_class, GenePhenotypeCategory, sep = "|", remove = TRUE, na.rm = TRUE) %>%
    rename(panel_class = temp_panel_class) %>%
    select(-gno2x_expected_an, -gno3_expected_an) %>% 
    mutate(gt_alt_freqs = round(gt_alt_freqs,2),
           aaf = round(aaf, 3),
           gno2x_af_all = round(gno2x_af_all, 5),
           gno3_af_all = round(gno3_af_all, 5),
           max_af = round(max_af, 5),
           af_oglx = round(af_oglx, 5),
           af_oglg = round(af_oglg, 5),
           pli = round(pli, 3),
           loeuf = round(loeuf, 2),
           mis_z = round(mis_z, 2),
           gno2e3g_af = round(gno2e3g_af,5),
           gnomad_nc_constraint = round(gnomad_nc_constraint,2) ) %>% 
    #unite("hgmd", hgmd_class, hgmd_id, na.rm = TRUE, remove = FALSE) %>%
    mutate(
      temp_hgvsc = if_else(is.na(hgvsc), NA_character_, str_replace(hgvsc, "^[^:]*:", "")),
      temp_prot  = if_else(is.na(hgvsp), NA_character_, str_replace(hgvsp, "^[^:]*:p.", "")),
      temp_prot = sub("%3D", "=", temp_prot),
      hgvs_trscrpt = str_extract(hgvsc, "^[^:]+"),
      vep_gene = if_else(is.na(hgvs_trscrpt), gene, str_c(gene, " ", hgvs_trscrpt, sep = "")),
      vep_hgvs = if_else(is.na(temp_prot), temp_hgvsc, str_c(temp_hgvsc, " p.(", temp_prot, ")", sep = "")),
    ) %>%
    mutate(gnomad_af = if_else(is.na(gno2x_af_all) | gno2x_af_all == -1, gno3_af_all, gno2x_af_all)) %>% 
    select(-temp_hgvsc, -temp_prot, -hgvs_trscrpt) %>%
    #mutate(zygosity = ifelse(gt_types == 1, "Heterozygous", ifelse(chrom %in% c("chrX", "X"), "Hemizygous", "Homozygous"))) %>% 
    select('ref_gene','sample','chr_variant_id',gts,'hgmd_class',clnsig,clnsigconf,'refgenewithver',vep_gene,vep_hgvs,gt_types,hgmd_id,clnid,gnomad_af,'gno2e3g_acan', 'gno2e3g_hom',gno2e3g_hemi,
           'panel_class','caller', 'priority_score', 'prscore_intervar', 'clinvar_hgmd_score', 'splice_score', 'insilico_score', 
           'omim_inheritance','omim_phen','hg38_pos', 'note','aaf','oe_lof_upper_bin','pli','loeuf','mis_z',exonicfunc_refgenewithver, 'exon','aa_length','intron',
           'grch37variant_id',gt_alt_freqs,gt_depths,gt_quals, 
           'max_af', 'max_af_pops',pmaxaf, af_oglx,af_oglg,gno2x_af_all,gno3_af_all,'interpro_domain', 'pfam_domain','pnull','prec', 'rmsk', 
           promoterai, 'am_pathogenicity','am_class','spliceai','spliceai_maxscore','spliceaimasked50','spliceaimasked50max',
           'qual',gno2x_filter,gno3_filter,func_refgenewithver,'gene',mane_select,'hgvsc','hgvsp','omim_gene', 'hgmd_phen', hgmd_overlap4aa, 'existing_variation', clnalleleid,'clin_sig', clnrevstat, clndn, clndisdb, 
           intervar_and_evidence, 'pvs1', 'truncating_vep', 'gno2e3g_af', ac_oglx, ac_hom_oglx, an_oglx, ac_oglg, ac_hom_oglg, an_oglg, 'existing_inframe_oorfs','existing_outofframe_oorfs','existing_uorfs','five_prime_utr_variant_annotation','five_prime_utr_variant_consequence',
           'squirls_interpretation', 'squirls_maxscore', 'squirls_score', 'dbscsnv_ada_score', 'dbscsnv_rf_score', 'dpsi_max_tissue', 'dpsi_zscore', 'genesplicer', 'maxentscan_diff', 'branchpoint_prob', 'labranchor_score', 'regsnp_fpr','regsnp_disease','regsnp_splicing_site',  
           'sift_pred', 'polyphen_pred', 'mutscore', 'mutationassessor_pred', 'mutationtaster_pred', 'metasvm_pred','clinpred_score', 'primateai_rankscore', 'revel_score', hmc_score, 'ccr_pct','mpc_score', 'mtr_score', 'mtr_fdr', 'mtr_pct', 'cadd_raw', 'cadd_phred','remm', 'fathmm_xf_coding_score','fathmm_xf_noncoding','eigen_pc_raw_coding', 'gerpplus_rs', 'phylop100way_vertebrate', 
           gnomad_nc_constraint, 'atac_rpe_score','atac_rpe_itemrgb', 'ft_ret_rpe_score', cherry_sum_score, 'gene_refgenewithver', 'avsnp150', 'tfbs',  'sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01',  
           'eyeintegration_rpe_adulttissue', 'eyeintegration_rpe_cellline', 'eyeintegration_rpe_fetaltissue', 'eyeintegration_rpe_stemcellline', 'eyeintegration_retina_adulttissue', 'eyeintegration_retina_stemcellline', 'eyeintegration_wholeblood', 
           'pubmed','sift_score', 'polyphen_score', 'metasvm_rankscore', 'metasvm_score', 'provean_score', 'provean_converted_rankscore',	'provean_pred',
           'f1000g2015aug_all','esp6500siv2_all',gno2_xg_ratio:gno3_popmax, 'syn_z', everything() ) %>% 
    arrange(desc(eyeGene), ref_gene)
}

# select('ref_gene','sample','chr_variant_id','hg38_pos', gts,gt_types,gt_alt_freqs, 'aaf', 'caller',
#          'panel_class', 'priority_score', 'prscore_intervar', 'clinvar_hgmd_score', 'splice_score', 'insilico_score', 
#          exonicfunc_refgenewithver, 'refgenewithver',mane_select,'gene','hgvsc','hgvsp','exon','aa_length','intron','omim_inheritance','omim_phen',hgmd,clnsig,clnsigconf,'note', 
#          gno2x_af_all,gno3_af_all,'gno2e3g_acan', 'gno2e3g_hom',gno2e3g_hemi,'max_af', 'max_af_pops',af_oglx,af_oglg,'interpro_domain', 'pfam_domain', 
#          'oe_lof_upper_bin','pli','loeuf','mis_z','pnull','prec', 'rmsk', 
#          promoterai, 'am_pathogenicity','am_class','spliceai','spliceai_maxscore','spliceaimasked50','spliceaimasked50max',
#          'grch37variant_id','qual',gt_depths,gt_quals,gno2x_filter,gno3_filter,func_refgenewithver,'omim_gene','hgmd_id','hgmd_class', 'hgmd_phen', hgmd_overlap4aa, 'existing_variation',clnid, clnalleleid,'clin_sig', clnrevstat, clndn, clndisdb, 
#          intervar_and_evidence, 'pvs1', 'truncating_vep', 'gno2e3g_af','pmaxaf', ac_oglx, ac_hom_oglx, an_oglx, ac_oglg, ac_hom_oglg, an_oglg, 'existing_inframe_oorfs','existing_outofframe_oorfs','existing_uorfs','five_prime_utr_variant_annotation','five_prime_utr_variant_consequence',
#          'squirls_interpretation', 'squirls_maxscore', 'squirls_score', 'dbscsnv_ada_score', 'dbscsnv_rf_score', 'regsnp_fpr','regsnp_disease','regsnp_splicing_site','dpsi_max_tissue', 'dpsi_zscore', 'genesplicer', 'maxentscan_diff', 'branchpoint_prob', 'labranchor_score', 'regsnp_fpr','regsnp_disease','regsnp_splicing_site',  
#          'sift_pred', 'polyphen_pred', 'mutscore', 'mutationassessor_pred', 'mutationtaster_pred', 'metasvm_pred','metasvm_score', 'clinpred_score', 'primateai_rankscore', 'revel_score', hmc_score, 'ccr_pct','mpc_score', 'mtr_score', 'mtr_fdr', 'mtr_pct', 'cadd_raw', 'cadd_phred','remm', 'fathmm_xf_coding_score','fathmm_xf_noncoding','eigen_pc_raw_coding', 'gerpplus_rs', 'phylop100way_vertebrate', 
#          gnomad_nc_constraint, 'atac_rpe_score','atac_rpe_itemrgb', 'ft_ret_rpe_score', cherry_sum_score, 'gene_refgenewithver', 'avsnp150', 'tfbs',  'sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01',  
#          'eyeintegration_rpe_adulttissue', 'eyeintegration_rpe_cellline', 'eyeintegration_rpe_fetaltissue', 'eyeintegration_rpe_stemcellline', 'eyeintegration_retina_adulttissue', 'eyeintegration_retina_stemcellline', 'eyeintegration_wholeblood', 
#          'pubmed','sift_score', 'polyphen_score', 'metasvm_rankscore', 'metasvm_score', 'provean_score', 'provean_converted_rankscore',	'provean_pred',
#          'f1000g2015aug_all','esp6500siv2_all',gno2_xg_ratio:gno3_popmax, 'syn_z', everything() )
# mutate(temp_genes_bed = pmap_chr(list(eyeintegration_gene, gene_gnomad, omim_gene, gene, gene_refgenewithver), ~toString(unique(na.omit(c(...)))) )) %>%
#   mutate(temp_genes_bed = na_if(temp_genes_bed, "") ) %>% 
#   mutate(ref_gene = ifelse(is.na(ref_gene), temp_genes_bed, ref_gene)) %>%  

gemini_filtered <- gemini_rearrangeCol %>% mutate(temp_group = ifelse(priority_score >= 3, 3, ifelse(priority_score >= -3, -3, -4))) %>% # checked OPA1 non-coding regions, the AF for some of variants with score -3 are around 0.01.
  filter(!ref_gene %in% blacklistGene, priority_score >= 15 | (temp_group >= -3 & pmaxaf < 0.1 & aaf < aafCutoff & af_oglg < 0.05 & af_oglx < 0.05) ) %>%
  arrange(desc(eyeGene), desc(temp_group), desc(maxpriorityscore), ref_gene, desc(priority_score)) %>% 
  filter(priority_score >= 3)

gemini_filtered0 <- gemini_filtered %>% select(-maxpriorityscore) #filtered0: temp_group >= -3 & pmaxaf < 0.2 & aaf < aafCutoff) | priority_score >= 10 

#AnnotSV does not have GD_POPMAX_AF column as of 3/2/2021, if GD_POPMAX_AF in columnames, then use it in the next version.
write_tsv(gemini_filtered0, file.path('.', filteredGemini_tsv_file), na="") #gemini_filtered tsv file main filter is priority_score >= 3

#found out whether blank is 0 after importing, use 1197 for filtering
#consider adding manta to the main gemini df for sorting/filtering after knowing the specificity of the manta calls. To better sort AR, AD, and ACMG 2nd.

manta_freq <- read_tsv(manta_freq_file, col_names = TRUE, na = c("NA", "full=NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
  type_convert()

if ( !file.exists(manta_file)) {
  manta_sort <- data.frame("sample" = sampleName, "note" = "No manta calling")
} else if (file.size(manta_file) == 0) {
  manta_sort <- data.frame("sample" = sampleName, "note" = "manta empty or not performed")
} else {
  manta_original <- read_tsv(manta_file, col_names = TRUE, na = c("NA", "full=NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
    filter(FILTER == 'PASS') %>% 
    mutate(ACMG_class = sub("full=", "", ACMG_class)) %>% 
    type_convert() 
  if (nrow(manta_original) == 0) {
    manta_sort <- data.frame("sample" = sampleName, "note" = "No PASS manta call") 
    } else { manta1 <- manta_original %>% mutate(temp_SV_start = round(SV_start, -3), temp_SV_end = round(SV_end, -3)) %>%
      unite("variant", SV_chrom, temp_SV_start, temp_SV_end, SV_type, sep = "-", remove = FALSE) %>% 
      left_join(., manta_freq, by = c("variant") ) %>% 
      filter(is.na(CohortFreq) | CohortFreq < 0.025) %>% 
      filter( Annotation_mode == "split" | Gene_count == 0)
    if (nrow(manta1) == 0) {
      manta_sort <- data.frame("sample" = sampleName, "note" = "No rare manta call") 
    } else {
      manta <- manta1 %>%
        mutate(ACMG_class = case_when(!is.na(ACMG_class) & SV_type == "DEL" & !is.na(B_loss_source) ~ ACMG_class - 2,
                                      !is.na(ACMG_class) & SV_type == "DUP" & !is.na(B_gain_source) ~ ACMG_class - 1,
                                      SV_type == "BND" & ( ( SV_chrom %in% c("X", "chrX") & between(SV_start, 140420272, 140421278) | ( grepl("chrX|X", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 140420272, 140421278) ) )) ~ 5,
                                      ( SV_chrom %in% c("17", "chr17") & between(SV_start, 59105674, 59683460) ) | ( grepl("chr17|17", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 59105674, 59683460) )  ~ 5,
                                      TRUE ~ ACMG_class ) ) %>% 
        filter(ACMG_class > 1 | is.na(ACMG_class), is.na(SV_length) | abs(SV_length) < 1000000) %>% #added is.na(ACMG_class) 11/17/2021, may need to remove this part if too many lines in the results
        separate(Location, c('temp_location1', 'temp_location2'), sep = "-", remove = FALSE, convert = FALSE) %>% 
        filter(!(grepl("intron", temp_location1) & temp_location1 == temp_location2 & (!is.na(B_gain_source) | !is.na(B_loss_source) | !is.na(B_ins_source)) )) %>% #remove del,dup,ins,bnd if only in one intron and present in benign databases
        filter(!(grepl("intron", temp_location1) & temp_location1 == temp_location2 & SV_type %in% c("DEL", "INS", "DUP") )) %>% #remove del,dup,ins if only in one intron even not found in benign database
        select(-starts_with('temp_'), -variant) 
       manta_sort <- left_join(manta, panelGene, by = c("Gene_name" = "ref_gene")) %>%
        mutate(note = "") %>% 
        mutate(panel_class == ifelse(SV_type == "BND" & ( ( SV_chrom %in% c("X", "chrX") & between(SV_start, 140420272, 140421278) | ( grepl("chrX|X", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 140420272, 140421278) ) )), "Dx", panel_class)) %>% 
        mutate(eyeGene = case_when(panel_class == "Dx" ~ 2,
                                   panel_class == "Candidate-High" ~ 1.6,
                                   panel_class == "Candidate-refutedDxGene" ~ 1.4,
                                   panel_class == "Candidate" ~ 1,
                                   panel_class == "Dx-modifier-rare" ~ 0.8,
                                   panel_class == "Dx-modifier-common" ~ 0.6,
                                   SV_type == "BND" & ( ( SV_chrom %in% c("X", "chrX") & between(SV_start, 140420272, 140421278) | ( grepl("chrX|X", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 140420272, 140421278) ) )) ~ 3,
                                   ( SV_chrom %in% c("17", "chr17") & between(SV_start, 59105674, 59683460) ) | ( grepl("chr17|17", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 59105674, 59683460) )  ~ 2.5,
                                   TRUE ~ 0)) %>% # GRCh38 coordinates
        arrange(desc(eyeGene), desc(ACMG_class)) %>% 
        unite("temp_panel_class", panel_class, GenePhenotypeCategory, sep = "|", remove = TRUE, na.rm = TRUE) %>%
        rename(panel_class = temp_panel_class) %>% 
        mutate(note = ifelse(SV_type == "BND" & ( ( SV_chrom %in% c("X", "chrX") & between(SV_start, 140420272, 140421278) | ( grepl("chrX|X", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 140420272, 140421278) ) )), "XLFD", note)) %>% 
        mutate(note = ifelse(( SV_chrom %in% c("17", "chr17") & between(SV_start, 59105674, 59683460) ) | ( grepl("chr17|17", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 59105674, 59683460) ), "RP17", note)) %>% 
        select(AnnotSV_ID:Gene_name,panel_class,ACMG_class,note,CohortFreq,`NaltP/NtotalP`,OMIM_phenotype:GnomAD_pLI,Frameshift,everything() )
    }
      }
}
#SV_type TRA?? in manta call.

jaxcnv_freq <- read_tsv(jaxcnv_freq_file, col_names = TRUE, na = c("NA", "full=NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
  type_convert()

if ( !file.exists(jaxcnv_file)) {
  jaxcnv_sort <- data.frame("sample" = sampleName, "note" = "No jaxcnv calling")
} else if (file.size(jaxcnv_file) == 0) {
  jaxcnv_sort <- data.frame("sample" = sampleName, "note" = "Empty jaxcnv")
} else {
  jaxcnv_original <- read_tsv(jaxcnv_file, col_names = TRUE, na = c("NA", "full=NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
    mutate(ACMG_class = sub("full=", "", ACMG_class)) %>% 
    type_convert()
  trim_cell <- function(x) {
    if (is.character(x)) substr(x, 1L, 32767L) else x
  }
  jaxcnv_original[] <- lapply(jaxcnv_original, trim_cell)
  if (nrow(jaxcnv_original) == 0) {
    jaxcnv_sort <- data.frame("sample" = sampleName, "note" = "No PASS jaxcnv call") 
  } else { jaxcnv1 <- jaxcnv_original %>% mutate(temp_SV_start = round(SV_start, -3), temp_SV_end = round(SV_end, -3)) %>%
    unite("variant", SV_chrom, temp_SV_start, temp_SV_end, SVtype, sep = "-", remove = FALSE) %>% 
    left_join(., jaxcnv_freq, by = c("variant") ) %>% 
    filter(is.na(CohortFreq) | CohortFreq < 0.025) #%>% 
   # filter( Annotation_mode == "split" | Gene_count == 0)
  if (nrow(jaxcnv1) == 0) {
    jaxcnv_sort <- data.frame("sample" = sampleName, "note" = "No rare jaxcnv call") 
  } else {
    jaxcnv <- jaxcnv1 %>%
      mutate(ACMG_class = case_when(!is.na(ACMG_class) & SVtype == "DEL" & !is.na(B_loss_source) ~ ACMG_class - 2,
                                    !is.na(ACMG_class) & SVtype == "DUP" & !is.na(B_gain_source) ~ ACMG_class - 1,
                                    SV_chrom %in% c("17", "chr17") & ( between(SV_start, 59105674, 59683460) | between(SV_end, 59105674, 59683460) ) ~ 5,
                                    TRUE ~ ACMG_class ) )
    jaxcnv_sort <- left_join(jaxcnv, panelGene, by = c("Gene_name" = "ref_gene")) %>% 
      mutate(note = "") %>% 
      mutate(panel_class == ifelse(SV_chrom %in% c("17", "chr17") & ( between(SV_start, 59105674, 59683460) | between(SV_end, 59105674, 59683460) ), "Dx", panel_class)) %>% 
      mutate(eyeGene = case_when(panel_class == "Dx" ~ 2,
                                 panel_class == "Candidate-High" ~ 1.6,
                                 panel_class == "Candidate-refutedDxGene" ~ 1.4,
                                 panel_class == "Candidate" ~ 1,
                                 panel_class == "Dx-modifier-rare" ~ 0.8,
                                 panel_class == "Dx-modifier-common" ~ 0.6,
                                 SV_chrom %in% c("17", "chr17") & ( between(SV_start, 59105674, 59683460) | between(SV_end, 59105674, 59683460) ) ~ 3,
                                 TRUE ~ 0)) %>% # GRCh38 coordinates
      arrange(desc(eyeGene), desc(ACMG_class)) %>%
      unite("temp_panel_class", panel_class, GenePhenotypeCategory, sep = "|", remove = TRUE, na.rm = TRUE) %>%
      rename(panel_class = temp_panel_class) %>% 
      mutate(note = ifelse(SV_chrom %in% c("17", "chr17") & ( between(SV_start, 59105674, 59683460) | between(SV_end, 59105674, 59683460) ), "RP17", note)) %>% 
      select(AnnotSV_ID:Gene_name,panel_class,ACMG_class,note,CohortFreq,`NaltP/NtotalP`,OGLsamples,OMIM_phenotype:GnomAD_pLI,Frameshift,everything() )
  }
  }
}

if ( roh_file == "filePlaceholder") {
  roh <- data.frame("sample" = sampleName, "note" = "Roh not analyzed.")
} else if (file.size(roh_file) == 0) { roh <- data.frame("sample" = sampleName, "note" = "Empty roh") 
} else {
  roh <- read_tsv(roh_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
    type_convert() %>% 
    mutate(autosome = ifelse(`#Chr` %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
                                           "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22" ), "autosome", "X")) %>% 
    group_by(autosome) %>% 
    mutate(pct = sum(`Size(Mb)`)/2900) %>% 
    ungroup() %>% 
    mutate(pct = round(pct, digits = 3)) %>% 
    select(-autosome)
}

if ( scramble_mei_file == "filePlaceholder") {
  scramble_mei <- data.frame("sample" = sampleName, "note" = "Scramble mei not analyzed.")
} else if (file.size(scramble_mei_file) == 0) {
  scramble_mei <- data.frame("sample" = sampleName, "note" = "Empty scramble mei or not performed")
} else {
  scramble_mei <- read_xlsx(scramble_mei_file, na = c("NA", "", "None", "NONE", "."))
}

if (scramble_del_file == "filePlaceholder") {
  scramble_del_sort <- data.frame("sample" = sampleName, "note" = "Scramble del not performed")
} else if (file.size(scramble_del_file) == 0) {
  scramble_del_sort <- data.frame("sample" = sampleName, "note" = "Empty scramble del")
} else {
  scramble_del <- read_tsv(scramble_del_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
    mutate(ACMG_class = sub("full=", "", ACMG_class)) %>% 
    type_convert() %>% 
    filter( is.na(Gene_name) | Annotation_mode == 'split' & !is.na(Gene_name) ) %>%
    filter(is.na(CohortFreq) | CohortFreq < 0.025) %>% 
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
    mutate(eyeGene = case_when(panel_class == "Dx" ~ 9,
                               panel_class == "Candidate-High" ~ 8,
                               panel_class == "Candidate-refutedDxGene" ~ 7,
                               panel_class == "Candidate" ~ 6,
                               panel_class == "Dx-modifier-rare" ~ 3,
                               panel_class == "Dx-modifier-common" ~ 2,
                               TRUE ~ 0)) %>% 
    arrange(desc(eyeGene), desc(ACMG_class)) %>% 
    unite("temp_panel_class", panel_class, GenePhenotypeCategory, sep = "|", remove = TRUE, na.rm = TRUE) %>%
    rename(panel_class = temp_panel_class) %>% 
    select(AnnotSV_ID:Annotation_mode,OMIM_phenotype:ACMG_class,panel_class,CohortFreq,`NaltP/NtotalP`,note,everything() ) 
  if (dim(scramble_del_sort)[1] == 0) {
    scramble_del_sort <- scramble_del_sort %>% add_row(note = "no scramble del candidate after filtering with scramble db")
  } 
}
  
#deleted sampleName, after FORMAT, add this back for production

#manta_file "Z:/NextSeqAnalysis/test2/manta/manta.1197.annotated.tsv"
##Add clinSV for genome, the position is "filePlaceholder" for other analysis type.


gemini_filtered <- gemini_filtered %>% 
  mutate(temp_ref_gene = case_when(ref_gene %in% c("ROM1", "PRPH2") ~ "PRPH2-ROM1",
                                   ref_gene %in% c("PCDH15", "CDH23") ~ "PCDH15-CDH23",
                                   ref_gene %in% c("CNGA3", "CNGB3") ~ "CNGA3-CNGB3",
                                   ref_gene %in% c("CNGA1", "CNGB1") ~ "CNGA1-CNGB1",
                                   TRUE ~ ref_gene)) %>% # digenic recessive
 filter(gt_quals > 10) #added 2/6/2026
gemini_filtered1 <- gemini_filtered %>% filter( priority_score > 4.5 | 
                                                ( between(priority_score, 3, 4.5) & (clinvar_hgmd_score >= 3 | splice_score >= 4 | insilico_score >= 3))
                                              ) %>% # changed from score of 3 to 5 on 9/11/2025 
  select(-maxpriorityscore) %>% 
  arrange(desc(eyeGene), desc(priority_score)) #this goes to "all" sheet after updating clinSV genes later.
  
gemini_filtered2 <- gemini_filtered %>% filter( priority_score > 5.5 | 
                                                ( between(priority_score, 3, 5.5) & (clinvar_hgmd_score >= 3 | splice_score >= 6 | insilico_score >= 3))
) # changed from score of 4 to 5 on 9/11/2025 

recessive_count <- select(gemini_filtered2, c(temp_ref_gene, gt_types)) %>%
  group_by(temp_ref_gene) %>% summarize(recessive_cnt = sum(gt_types)) #filtered further than gemini_filtered2 prior v3.4.4

eyeGeneList <- read_xlsx(geneCategory_file, sheet = "analysis", na = c("NA", "", "None", ".")) %>% 
  select(gene) %>% distinct() %>% 
  pull(gene)
acmg_genes <- read_xlsx(geneCategory_file, sheet = "ACMG", na = c("NA", "", "None", "NONE", ".")) %>% pull(Gene) %>% unique()
acmg_ARgenes <- read_xlsx(geneCategory_file, sheet = "ACMG", na = c("NA", "", "None", "NONE", ".")) %>% 
  filter(Inheritance == "AR") %>% pull(Gene) %>% unique()

selectGene <- function(x, y){
  x = as.character(x)
  geneNames <- as.list(strsplit(x, ","))[[1]] 
  if (length(geneNames) > 0) {
    eyeGene <- purrr::keep(geneNames, geneNames %in% y) 
    if (length(eyeGene) == 0) {
      return(NA)
    } else if (length(eyeGene) == 1) { return(eyeGene) }
    else { return(paste(eyeGene, collapse = ",")) }
  } else {
    return(NA)
  }
}

#clinSV cds gene added to reccessive_count
if ( clinsv_file == "filePlaceholder") {
  clinsv <- data.frame("sample" = sampleName, "note" = "ClinSV not analyzed.")
} else if (file.size(clinsv_file) == 0) {
  clinsv <- data.frame("sample" = sampleName, "note" = "Empty clinsv or not performed")
} else {
  clinsv <- read_xlsx(clinsv_file, sheet = "clinSV", na = c("NA", "", "None", "NONE", ".")) %>% 
    filter(is.na(VariantID) | VariantID != "ID") #temporary fix for extra column name line, to be updated/edited
  clinsv$ACMG2nd <- sapply(1:nrow(clinsv), function(x) {selectGene(clinsv[x, "Genes"], acmg_genes)})
  clinsv$eyeGene <- sapply(1:nrow(clinsv), function(x) {selectGene(clinsv[x, "Genes"], eyeGeneList)})
  clinsv <- select(clinsv, Sample:SEGD, eyeGene, ACMG2nd, everything()) %>%
    filter(PopAF_MGRB < 0.05, PopAF1k < 0.05) %>% 
    arrange(eyeGene, ACMG2nd) %>% 
    left_join(., panelGene, by = c("eyeGene" = "ref_gene")) %>% 
    unite("temp_panel_class", panel_class, GenePhenotypeCategory, sep = "|", remove = TRUE, na.rm = TRUE) %>%
    rename(panel_class = temp_panel_class) %>% 
    select(Sample:SEGD, eyeGene, panel_class, everything())
  clinsv_cds <- filter(clinsv, SVtype == "BND" | GeneFeature != "intron") %>% 
    select(Genes, GT) %>% separate_rows(Genes, sep=",") %>% unique() %>% 
    separate(GT, c("GT1", "GT2"), sep = "\\/", convert = TRUE) %>% 
    mutate(clinsvGT = GT1 + GT2) %>% 
    select(Genes, clinsvGT) %>% 
    mutate(temp_ref_gene = case_when(Genes %in% c("ROM1", "PRPH2") ~ "PRPH2-ROM1",
                                     Genes %in% c("PCDH15", "CDH23") ~ "PCDH15-CDH23",
                                     Genes %in% c("CNGA3", "CNGB3") ~ "CNGA3-CNGB3",
                                     Genes %in% c("CNGA1", "CNGB1") ~ "CNGA1-CNGB1",
                                     TRUE ~ Genes)) # digenic recessive
  clinsv_cds_geneList <- select(clinsv_cds, temp_ref_gene) %>% 
    distinct() %>% pull(temp_ref_gene)
  recessive_count <- left_join(recessive_count, clinsv_cds, by = c("temp_ref_gene")) %>% 
    replace_na(list(clinsvGT = 0)) %>% 
    mutate(recessive_cnt = recessive_cnt + clinsvGT) %>% 
    select(temp_ref_gene, recessive_cnt)
  gemini_filtered1 <- gemini_filtered1 %>% 
    mutate(note = ifelse(temp_ref_gene %in% clinsv_cds_geneList, ifelse(note == "", "clinSV", paste0(note, "; clinSV")), note)) %>% 
    select(-temp_ref_gene)
  gemini_filtered2 <- gemini_filtered2 %>% 
    mutate(note = ifelse(temp_ref_gene %in% clinsv_cds_geneList, ifelse(note == "", "clinSV", paste0(note, "; clinSV")), note))
}

gemini_filtered3 <- left_join(gemini_filtered2, recessive_count, by=c("temp_ref_gene")) %>% 
  replace_na(list(recessive_cnt=0)) %>% 
  mutate(recessive_cnt = as.integer(recessive_cnt)) %>% 
  select(-temp_ref_gene)

acmg <- gemini_filtered3 %>% filter(ref_gene %in% acmg_genes, priority_score > 4) %>%
  filter(ref_gene != "HFE" | (chr_variant_id == "chr6-26092913-G-A" & recessive_cnt > 1) ) %>% #HFE report hmz C282Y only
  filter(recessive_cnt > 1 | !ref_gene %in% acmg_ARgenes | 
           grepl("RPE65",refgenewithver) & grepl("D477G|E519K",refgenewithver)) %>% #AR gene requires two variants, RPE65 has 2 dominant variants
  filter(ref_gene != "ABCD1" | recessive_cnt > 1 ) %>% #ABCD1 (XLR) requires hemi or hmz or 2 het
  filter(ref_gene != "TTN" | pvs1 == 1 | truncating_vep == 1 ) %>%  #Truncating only for TTN
  select(-maxpriorityscore, -recessive_cnt)
 
print("###acmg done### 70%")

xR <- gemini_filtered3 %>% filter(chrom %in% c("X", "chrX"), recessive_cnt >= 2) %>% select(-maxpriorityscore, -recessive_cnt) #removed priority_score >= 5.5, as already filtered above
xD <- gemini_filtered3 %>% filter(chrom %in% c("X", "chrX"), recessive_cnt == 1, pmaxaf < 0.001) %>% select(-maxpriorityscore, -recessive_cnt)

#ar are those genes with homozygous or compound hets variants of ps >= 5. However, ps = 4 variants were also listed if there are 2 ps>=5.
#ar gene with 1 hit will not be here.
ar <- gemini_filtered3 %>% filter(!chrom %in% c("X", "Y", "chrX", "chrY"), recessive_cnt >= 2) %>% #removed priority_score >= 5.5, as already filtered above
  mutate(knownAR = ifelse(grepl("AR", omim_inheritance), 1, 0)) %>% 
  arrange(desc(eyeGene),  desc(maxpriorityscore), ref_gene, desc(knownAR), desc(priority_score)) %>% 
  mutate(note = ifelse(ref_gene %in% c("PRPH2", "ROM1", "PCDH15", "CDH23", "CNGA1", "CNGB1", "CNGA3", "CNGB3"), 
                      ifelse(note == "", "Digenic?", paste0(note, "; Digenic?")), note) ) %>%
  select(-maxpriorityscore, -knownAR, -recessive_cnt) # digenic recessive

#ad are those omim unknown or AD inheritance. 
ad <- gemini_filtered3 %>% filter(!chrom %in% c("X", "Y", "chrX", "chrY"), 
                                  recessive_cnt == 1 & !omim_inheritance %in% c("AR") | recessive_cnt >=2 & grepl("AD", omim_inheritance),
                                  pmaxaf < 0.001) %>% #removed priority_score >= 5.5, as already filtered above, reduced pmaxaf cutoff from 0.002 to 0.001
  arrange(desc(eyeGene), desc(maxpriorityscore), ref_gene, desc(priority_score)) %>% 
  select(-maxpriorityscore, -recessive_cnt)
print("###inheritance search done### 80%")
#grepl("AD", omim_inheritance) | is.na(omim_inheritance)

#AD: score > 4 , AR: score > 4, all: score >= 3


config <- read_tsv(config_file, col_names = FALSE, na = c("NA", ""), col_types = cols(.default = col_character())) %>% 
  separate("X1", c("tool", "version", "note"), sep = "\\:|\\#", remove = TRUE)

if (is.na(mutserve_file) | mutserve_file == "filePlaceholder") {
  mutserve <- data.frame("sample" = sampleName, "note" = "Mitochondria not analyzed or empty.")
} else {
  mutserve <- read_tsv(mutserve_file, col_names = TRUE, na = c("NA", "full=NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
    type_convert() %>% 
    replace_na(list(Helix_vaf_hom = -1, Helix_vaf_het = -1)) %>% 
    mutate(VariantGroup = ifelse(Helix_vaf_hom > 0.01 | Helix_vaf_het > 0.01 | VariantLevel < 0.01, -1, 1) ) %>% 
    mutate(VariantGroup = case_when(grepl("P]", mitomapDiseaseStatus) ~ VariantGroup + 4,
                                    grepl("B]", mitomapDiseaseStatus) ~ VariantGroup - 4,
                                    !is.na(mitomapDiseaseStatus) ~ VariantGroup + 1,
                                    TRUE ~ VariantGroup)) %>% 
    arrange(desc(VariantGroup))
}
summaryInfo <- data.frame("sample" = sampleName, "PatientDxPhenotype" = NA, "DxOutcome"= NA, "variant" = NA, "reviewer" = NA, "date" = NA) %>% 
  add_row("sample" = sampleName, "PatientDxPhenotype" = NA, "DxOutcome"= NA, "variant" = NA, "reviewer" = NA, "date" = NA)

#scramble_del_file == "filePlaceholder", genome as scramble_del is not run for genome.
#is.na(convading_file): exome 
if (scramble_del_file == "filePlaceholder") {
  openxlsx::write.xlsx(list("AR" = ar, "AD" = ad, "XR" = xR, "XD" = xD, "MT" = mutserve, "ACMG" = acmg, "all" = gemini_filtered1, "rareRef" = gemini_ref_var_rearrangeCol, "clinSV" = clinsv, "manta" = manta_sort, "jaxCNV" = jaxcnv_sort, "scramble_mei" = scramble_mei, "roh" = roh, "config" = config, "summary" = summaryInfo), file = gemini_xlsx_file, firstRow = TRUE, firstCol = TRUE, keepNA = FALSE)
} else if (is.na(convading_file)) {
  openxlsx::write.xlsx(list("AR" = ar, "AD" = ad, "XR" = xR, "XD" = xD, "ACMG" = acmg, "all" = gemini_filtered1, "rareRef" = gemini_ref_var_rearrangeCol, "manta" = manta_sort, "scramble_mei" = scramble_mei, "scramble_del" = scramble_del_sort, "roh" = roh, "config" = config, "summary" = summaryInfo), file = gemini_xlsx_file, firstRow = TRUE, firstCol = TRUE, keepNA = FALSE)
} else {
    cnv <- read_tsv(convading_file, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
      type_convert() 
    if (dim(cnv)[1] == 0) {
      openxlsx::write.xlsx(list("AR" = ar, "AD" = ad, "XR" = xR, "XD" = xD, "ACMG" = acmg, "all" = gemini_filtered1, "rareRef" = gemini_ref_var_rearrangeCol, "manta" = manta_sort, "scramble_mei" = scramble_mei, "scramble_del" = scramble_del_sort, "config" = config, "summary" = summaryInfo), file = gemini_xlsx_file, firstRow = TRUE, firstCol = TRUE, keepNA = FALSE)
    } else {
      cnv_gene <- as.list(distinct(cnv, GENE)[[1]])
      #cnv_gene <- dplyr::pull(cnv, GENE) #pull column as a vector
      cnv_variant <- gemini_rearrangeCol %>% 
        filter(ref_gene %in% cnv_gene) %>% 
        mutate(LAF = ifelse(gt_alt_freqs > 0.5, 1 - gt_alt_freqs, gt_alt_freqs)) %>% 
        select(-gt_alt_freqs) %>% 
        select('chr_variant_id', 'chrom', 'start', 'qual', 'filter', starts_with('gts'), starts_with('gt_'), 'LAF', 'ref_gene', 'exon', 'ref_gene',  
              'refgenewithver', 'exonicfunc_refgenewithver', 'hgvsc', 'hgvsp', 'type') %>% 
        rename(gene = ref_gene)
      openxlsx::write.xlsx(list("AR" = ar, "AD" = ad, "XR" = xR, "XD" = xD, "ACMG" = acmg, "all" = gemini_filtered1, "rareRef" = gemini_ref_var_rearrangeCol, "CoNVaDING" = cnv, "CNV-variant" = cnv_variant, "manta" = manta_sort, "scramble_mei" = scramble_mei, "scramble_del" = scramble_del_sort, "config" = config, "summary" = summaryInfo), file = gemini_xlsx_file, firstRow = TRUE, firstCol = TRUE, keepNA = FALSE)
      cnv_edit <- cnv %>% 
        mutate(START = START - 100, STOP = STOP + 100, type = "snp") %>% #padding of 100 nt
        gather(START:STOP, key = "datatype", value = "position") %>% 
        mutate(LAF = 0.54, gt_depths = ifelse(datatype == "START", 40, 200) ) %>% 
        rename(chrom = "CHR", gene = "GENE") %>% 
        select(chrom, gene, datatype, position, LAF, gt_depths, type)
      
      cnv_variant_edit <- cnv_variant %>% 
        select(chrom, start, gt_depths,	LAF, gene, type) %>% 
        rename(position = "start") %>% 
        mutate(datatype="LAF") %>% 
        filter(gt_depths >= 15)
      
      variantForPlot <- rbind(cnv_edit, cnv_variant_edit) %>% mutate(gene = as.factor(gene)) %>% 
        arrange(chrom, position) %>% 
        mutate(variant_no = row_number()) %>% 
        mutate(DepthGroup = case_when(gt_depths >= 100 ~ "DP>=100",
                                      gt_depths >= 30 ~ "DP30-99",
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

