args = commandArgs(trailingOnly=TRUE)

# argument handling
geneCategory_file <- args[1]
aaf_freq <- args[2]
output_xlsx <- args[3]
family_name <- args[4]
denovo_file <- args[5]
ad_file <- args[6]
ar_file <- args[7]
comphets_file <- args[8]
xdenovo_file <- args[9]
xd_file <- args[10]
xr_file <- args[11]
mendel_errors_file <- args[12]

library(tidyverse)
library(readxl)

sortFilterGemini <- function(fileName) {
  InheritanceTest <- read_tsv(fileName, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
    mutate(atac_rpe_score = sub(",", "_", atac_rpe_score)) %>% 
    type_convert() %>% mutate(exon = sub("^", " ", exon)) %>% mutate( start_vcf = start + 1 ) %>% 
    unite("chr_variant_id", chrom, start_vcf, ref, alt, sep = "-", remove = FALSE ) %>% 
    mutate(gene = toupper(gene)) %>% 
    filter(pmaxaf < 0.2, aaf < aaf_freq, !ref_gene %in% blacklistGene, qual > 10)
  # get max_priority_score for each gene
  InheritanceTest_max_priority_score <- select(InheritanceTest, c(ref_gene, priority_score)) %>% group_by(ref_gene) %>% summarize(maxpriorityscore = max(priority_score)) 
  #arrange by max_priority_score, then by gene, and priority score. None gene region?
  InheritanceTest_max_priority_score1 <- left_join(InheritanceTest, InheritanceTest_max_priority_score, by=c("ref_gene"))
  #VEP hg19 version's gene names are the same as in the IDT ordering design sheets. This is what used for left_join
  InheritanceTest_rearrangeCol <- left_join(InheritanceTest_max_priority_score1, panelGene, by = c("ref_gene")) %>% 
    mutate(note = "") %>% separate(vcf_id, c('caller', 'hg38_id'), sep = "_") %>% 
    mutate(hg38_pos = sub("[ACGT]*>[ACGT]*", "", hg38_id)) %>% 
    mutate(gnomad_hom = ifelse(is.na(gno_hom) & is.na(gnog_hom), NA, ifelse(is.na(gno_hom), 0, gno_hom) + ifelse(is.na(gnog_hom), 0, gnog_hom) ), 
           gnomad_ac = ifelse(is.na(gno_ac_all) & is.na(gnog_ac_all), NA, ifelse(is.na(gno_ac_all), 0, gno_ac_all) + ifelse(is.na(gnog_ac_all), 0, gnog_ac_all) ), 
           gnomad_an = ifelse(is.na(gno_an_all) & is.na(gnog_an_all), NA, ifelse(is.na(gno_an_all), 0, gno_an_all) + ifelse(is.na(gnog_an_all), 0, gnog_an_all) )) %>% 
    mutate(gnomad_af = gnomad_ac/gnomad_an) %>% 
    unite("gnomad_acan", gnomad_ac, gnomad_an, sep = "/", remove = TRUE) %>% 
    mutate(eyeGene = case_when(panel_class == "Dx" ~ 2,
                               panel_class == "Candidate" ~ 1,
                               TRUE ~ 0)) %>% 
    select('chr_variant_id', 'chrom', 'start_vcf', 'qual', 'filter', 'family_id','family_members', 'family_genotypes', 'samples', 'aaf', 'caller', 'hg38_pos',
           'panel_class', 'priority_score', 'priority_score_intervar', 'clinvar_hgmd_score', 'splice_score', 'other_predic_score', 'gnomad_af', 'gnomad_acan', 'pmaxaf', 'max_af', 'max_af_pops', 'gnomad_hom', 'ref_gene', 'note', 
           'exonicfunc_ensgene', 'refgenewithver', 'hgvsc', 'hgvsp', 'gene', 'exon', 'aa_length', 'omim_gene', 'omim_inheritance', 'omim_phenotype', 'pvs1', 'truncating_vep', 'hgmd_overlap', 'existing_variation', 'clinvar_intervar', 'intervar_and_evidence', 
           'clinvar_id', 'clinvar_pathogenic', 'clinvar_sig', 'clin_sig', 'existing_inframe_oorfs','existing_outofframe_oorfs','existing_uorfs','five_prime_utr_variant_annotation','five_prime_utr_variant_consequence',
           'spliceai', 'spliceai_maxscore', 'spliceaimasked50', 'spliceaimasked50max', 'squirls_interpretation', 'squirls_maxscore', 'squirls_score', 'dbscsnv_ada_score_intervar', 'dbscsnv_rf_score_intervar', 'dpsi_max_tissue_annovar', 'dpsi_zscore_annovar', 'genesplicer', 'maxentscan_diff', 'branchpoint_u2_binding_energy', 'branchpoint_prob', 'branchpoint_to_3prime', 
           'sift_pred', 'polyphen_pred', 'mutationassessor_pred', 'mutationtaster_pred', 'metasvm_pred','metasvm_score_intervar', 'clinpred_score', 'primatedl', 'revel', 'ccr_pct','mpc', 'cadd_raw', 'cadd_phred','remm', 'fathmm_xf_coding','fathmm_xf_noncoding','eigen_pc_raw', 'eigen_raw', 'gerp_rs_intervar', 'phylop46way_placental_intervar', 'phylop_100way', 
           'atac_rpe_score','atac_rpe_itemrgb', 'ft_ret_rpe_score', 'genedetail_ensgene', 'aachange_ensgene', 'gene_refgenewithver', 'func_refgene', 'func_refgenewithver', 'exonicfunc_refgenewithver', 'exonicfunc_refgene', 'avsnp150_annovar', 'interpro_domain_intervar', 
           'pfam_domain', 'tfbs', 'pli', 'lof_z', 'mis_z', 'sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01', 'atac_rpe_itemrgb', 'atac_rpe_score', 
           'eyeintegration_rnaseq_tpm_rpe_adulttissue', 'eyeintegration_rnaseq_tpm_rpe_cellline', 'eyeintegration_rnaseq_tpm_rpe_fetaltissue', 'eyeintegration_rnaseq_tpm_rpe_stemcellline', 'eyeintegration_rnaseq_tpm_retina_adulttissue', 'eyeintegration_rnaseq_tpm_retina_stemcellline', 'eyeintegration_rnaseq_tpm_wholeblood', 
           'pubmed', 'polyphen_score', 'sift_score', 'eigen_pc_phred', 'eigen_pc_raw_rankscore', 'eigen_phred',	'eigen_coding_or_noncoding', 'fathmm_converted_rankscore',	'fathmm_pred',	'fathmm_score',	'gerp',	'genocanyon_score',	'genocanyon_score_rankscore',	'linsight', 'lrt_omega', 'lrt_converted_rankscore', 'lrt_pred',	'lrt_score',
           'm_cap_pred', 'm_cap_rankscore',	'm_cap_score', 'metalr_pred',	'metalr_rankscore',	'metalr_score',	'metasvm_rankscore', 'metasvm_score', 'mutationassessor_uniprotid',	'mutationassessor_score',	'mutationassessor_score_rankscore',	'mutationtaster_converted_rankscore',	'mutationtaster_model',	'mutationtaster_score',	'provean_converted_rankscore',	'provean_pred',	'provean_score', 'vest3_rankscore',	'vest3_score',
           'cpg_island', 'gno_an_all', 'gno_af_all',	'gno_hom', 'gnog_an_all', 'gnog_af_all',	'gnog_hom', 'gno_af_asj', 'gnog_af_asj', 'gwas_pubmed_trait', 'pnull', 'precessive', 'rmsk', 'syn_z', 'eyeGene', 'maxpriorityscore') %>%
    arrange(desc(eyeGene), desc(maxpriorityscore), ref_gene, desc(priority_score)) %>% select(-maxpriorityscore, -eyeGene) 
  return(InheritanceTest_rearrangeCol)
}

panelGene <- read_xlsx(geneCategory_file, sheet = "analysis", na = c("NA", "", "None", ".")) %>% select(gene, panel_class) %>% rename(ref_gene = gene)
blacklistGene <- read_xlsx(geneCategory_file, sheet = "IVA", na = c("NA", "", "None", "."))  %>% filter(Blacklist == "Excluded") %>% pull(Gene)

all <- data.frame()
if (file.size(denovo_file) == 0) {
  denovo <- data.frame("family_id" = family_name, "note" = "Empty denovo query")
} else {
  denovo <- sortFilterGemini(denovo_file)
  all <- rbind(all, denovo)
}

if (file.size(ad_file) == 0) {
  ad <- data.frame("family_id" = family_name, "note" = "Empty ad query")
} else {
  ad <- sortFilterGemini(ad_file) %>% filter(!chrom %in% c("X", "chrX"))
  all <- rbind(all, ad)
}

if (file.size(ar_file) == 0) {
  ar <- data.frame("family_id" = family_name, "note" = "Empty ar query")
} else {
  ar <- sortFilterGemini(ar_file) %>% filter(!chrom %in% c("X", "chrX"))
  all <- rbind(all, ar)
}

if (file.size(comphets_file) == 0) {
  comphets <- data.frame("family_id" = family_name, "note" = "Empty comphets query")
} else {
  comphets <- sortFilterGemini(comphets_file) %>% 
    filter(!is.na(gene)) %>% 
    distinct(chr_variant_id, .keep_all = TRUE)
  all <- rbind(all, comphets)
}

if (file.size(xdenovo_file) == 0) {
  xdenovo <- data.frame("family_id" = family_name, "note" = "Empty xdenovo query")
} else {
  xdenovo <- sortFilterGemini(xdenovo_file)
  all <- rbind(all, xdenovo)
}

if (file.size(xd_file) == 0) {
  xd <- data.frame("family_id" = family_name, "note" = "Empty xd query")
} else {
  xd <- sortFilterGemini(xd_file)
  all <- rbind(all, xd)
}

if (file.size(xr_file) == 0) {
  xr <- data.frame("family_id" = family_name, "note" = "Empty xr query")
} else {
  xr <- sortFilterGemini(xr_file)
  all <- rbind(all, xr)
}

if (file.size(mendel_errors_file) == 0) {
  mendel_errors <- data.frame("family_id" = family_name, "note" = "Empty mendel_errors query")
} else {
  mendel_errors <- sortFilterGemini(mendel_errors_file)
  all <- rbind(all, mendel_errors)
}

summaryInfo <- data.frame("family_id" = family_name, "DxOutcome"= NA, "variant" = NA, "reviewer" = NA, "date" = NA, "SecondReviewer" = NA, "SecondReviewDate" = NA)

acmg_genes <- read_xlsx(geneCategory_file, sheet = "ACMG3", na = c("NA", "", "None", ".")) %>% pull(Gene) %>% unique()
acmg <- all %>% filter(ref_gene %in% acmg_genes, priority_score > 4) %>% distinct(chr_variant_id, .keep_all = TRUE)
#in the next version, consider TTN truncating, HFE C..Y hmz etc.

openxlsx::write.xlsx(list( "de_novo" = denovo, 
                           "AD" = ad, 
                           "AR" = ar, 
                           "comp_hets" = comphets, 
                           "X_denovo" = xdenovo, 
                           "XD" = xd, 
                           "XR" = xr, 
                           "mendel_errors" = mendel_errors, 
                           "ACMG3" = acmg,
                           "summary" = summaryInfo), 
                     file = output_xlsx, firstRow = TRUE, firstCol = TRUE)
