args = commandArgs(trailingOnly=TRUE)
#args = c("Z:/resources/OGLpanelGeneDxORcandidate.xlsx", "0.8", "gemini_xlsx/D703.20220609.panelVeri.lenientYes.xlsx", "D703","temp/D703.20220609.panelVeri.denovo.tsv","temp/D703.20220609.panelVeri.ad.tsv","temp/D703.20220609.panelVeri.ar.tsv", "temp/D703.20220609.panelVeri.comphets.tsv", "temp/D703.20220609.panelVeri.xdenovo.tsv","temp/D703.20220609.panelVeri.xd.tsv", "temp/D703.20220609.panelVeri.xr.tsv", "temp/D703.20220609.panelVeri.mendel_errors.tsv", "config_variant_prioritization.yaml")
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
config_file <- args[13]

library(tidyverse)
library(readxl)

sortFilterGemini <- function(fileName) {
  InheritanceTest <- read_tsv(fileName, col_names = TRUE, na = c("NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
    mutate(atac_rpe_score = gsub(",", "_", atac_rpe_score)) %>% 
    type_convert() %>% mutate(exon = sub("^", " ", exon), intron = sub("^", " ", intron)) %>%
    mutate( start_vcf = start + 1 ) %>% 
    unite("chr_variant_id", chrom, start_vcf, ref, alt, sep = "-", remove = FALSE ) %>% 
    #mutate(gene = toupper(gene)) %>% 
    filter(pmaxaf < 0.2, aaf < aaf_freq, !ref_gene %in% blacklistGene, qual > 10)
  # get max_priority_score for each gene
  if (nrow(InheritanceTest) == 0) {
    return(data.frame()) 
  }
  else {
    InheritanceTest_max_priority_score <- select(InheritanceTest, c(ref_gene, priority_score)) %>% group_by(ref_gene) %>% summarize(maxpriorityscore = max(priority_score)) 
    #arrange by max_priority_score, then by gene, and priority score. None gene region?
    InheritanceTest_max_priority_score1 <- left_join(InheritanceTest, InheritanceTest_max_priority_score, by=c("ref_gene"))
    #VEP hg19 version's gene names are the same as in the IDT ordering design sheets. This is what used for left_join
    InheritanceTest_rearrangeCol <- left_join(InheritanceTest_max_priority_score1, panelGene, by = c("ref_gene")) %>% 
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
      mutate(gno2x_filter = ifelse(gno2x_an_all > 0 & is.na(gno2x_filter) & gno2x_an_all < gno2x_expected_an/2, "AN<half", gno2x_filter),
             gno3_filter = ifelse(gno3_an_all > 0 & is.na(gno3_filter) & gno3_an_all < gno3_expected_an/2, "AN<half", gno3_filter) ) %>%
      mutate(eyeGene = case_when(panel_class == "Dx" ~ 2,
                                 panel_class == "Candidate" ~ 1,
                                 TRUE ~ 0)) %>% 
      select(-gno2x_expected_an, -gno3_expected_an) %>%
      mutate(aaf = round(aaf, 3),
             gno2x_af_all = round(gno2x_af_all, 5),
             gno3_af_all = round(gno3_af_all, 5),
             max_af = round(max_af, 5),
             af_oglx = round(af_oglx, 5),
             af_oglg = round(af_oglg, 5),
             pli = round(pli, 3),
             loeuf = round(loeuf, 2),
             mis_z = round(mis_z, 2),
             gno2e3g_af = round(gno2e3g_af,5),
             gnomad_nc_constraint = round(gnomad_nc_constraint,2)) %>%
      select('ref_gene', 'chr_variant_id','grch37variant_id', 'family_id','family_members', 'family_genotypes', 'samples','aaf', 'caller',
             'panel_class', 'priority_score', 'prscore_intervar', 'clinvar_hgmd_score', 'splice_score', 'insilico_score', gno2x_af_all,gno3_af_all,'gno2e3g_acan', 'gno2e3g_hom','max_af', 'max_af_pops',af_oglx,af_oglg,'refgenewithver',
             'exon','aa_length','intron','omim_inheritance', 'omim_phen','hgmd_id',clnsig,clnsigconf,'oe_lof_upper_bin','pli','loeuf','mis_z', 'interpro_domain', 'pfam_domain', 'rmsk','note', 'spliceai', 'spliceai_maxscore', 'spliceaimasked50', 'spliceaimasked50max','hg38_pos', 'qual',
             gno2x_filter,gno3_filter,'gene', mane_select, 'hgvsc', 'hgvsp',func_refgenewithver, exonicfunc_refgenewithver, 'pnull','prec', 'omim_gene',  
             'hgmd_class', 'hgmd_phen', hgmd_overlap4aa, 'existing_variation',clnid, clnalleleid,'clin_sig', clnrevstat, clndn, clndisdb, 
             intervar_and_evidence, 'pvs1', 'truncating_vep', 'gno2e3g_af','pmaxaf', ac_oglx, ac_hom_oglx, an_oglx, ac_oglg, ac_hom_oglg, an_oglg, 'existing_inframe_oorfs','existing_outofframe_oorfs','existing_uorfs','five_prime_utr_variant_annotation','five_prime_utr_variant_consequence',
             'squirls_interpretation', 'squirls_maxscore', 'squirls_score', 'dbscsnv_ada_score', 'dbscsnv_rf_score', 'regsnp_fpr','regsnp_disease','regsnp_splicing_site','dpsi_max_tissue', 'dpsi_zscore', 'genesplicer', 'maxentscan_diff', 'branchpoint_prob', 'regsnp_fpr','regsnp_disease','regsnp_splicing_site',  
             'sift_pred', 'polyphen_pred', 'mutscore', 'mutationassessor_pred', 'mutationtaster_pred', 'metasvm_pred','metasvm_score', 'clinpred_score', 'primateai_rankscore', 'revel_score', hmc_score, 'ccr_pct','mpc_score', 'mtr_score', 'mtr_fdr', 'mtr_pct', 'cadd_raw', 'cadd_phred','remm', 'fathmm_xf_coding_score','fathmm_xf_noncoding','eigen_pc_raw_coding', 'gerpplus_rs', 'phylop100way_vertebrate', 
             gnomad_nc_constraint, 'atac_rpe_score','atac_rpe_itemrgb', 'ft_ret_rpe_score', cherry_sum_score, 'gene_refgenewithver', 'avsnp150', 'tfbs',  'sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01',
             'eyeintegration_rpe_adulttissue', 'eyeintegration_rpe_cellline', 'eyeintegration_rpe_fetaltissue', 'eyeintegration_rpe_stemcellline', 'eyeintegration_retina_adulttissue', 'eyeintegration_retina_stemcellline', 'eyeintegration_wholeblood', 
             'pubmed','sift_score', 'polyphen_score', 'metasvm_rankscore', 'metasvm_score', 'provean_score', 'provean_converted_rankscore',	'provean_pred',
             'f1000g2015aug_all','esp6500siv2_all',gno2_xg_ratio:gno3_popmax, 'syn_z', everything() ) %>%
      arrange(desc(eyeGene), desc(maxpriorityscore), ref_gene, desc(priority_score)) %>% select(-maxpriorityscore, -eyeGene) 
    return(InheritanceTest_rearrangeCol)
  }
}


#mutate(temp_genes_bed = pmap_chr(list(eyeintegration_gene, gene_gnomad, omim_gene, gene, gene_refgenewithver), ~toString(unique(na.omit(c(...)))) )) %>%
#mutate(temp_genes_bed = na_if(temp_genes_bed, "") ) %>% 
#  mutate(ref_gene = ifelse(is.na(ref_gene), temp_genes_bed, ref_gene)) %>%
#  select(-temp_genes_bed) %>%
  
panelGene <- read_xlsx(geneCategory_file, sheet = "analysis", na = c("NA", "", "None", "NONE", ".")) %>% select(gene, panel_class) %>% rename(ref_gene = gene)
blacklistGene <- read_xlsx(geneCategory_file, sheet = "IVA", na = c("NA", "", "None", "NONE", "."))  %>% filter(Blacklist == "Excluded") %>% pull(Gene)

all <- data.frame()
if (file.size(denovo_file) == 0) {
  denovo <- data.frame("family_id" = family_name, "note" = "Empty denovo query")
} else {
  denovo <- sortFilterGemini(denovo_file)
  if (dim(denovo)[1] > 0) { all <- rbind(all, denovo) }
}

if (file.size(ad_file) == 0) {
  ad <- data.frame("family_id" = family_name, "note" = "Empty ad query")
} else {
  ad <- sortFilterGemini(ad_file) 
  if (nrow(ad) > 0) { ad <- filter(ad, !chrom %in% c("X", "chrX"))
    all <- rbind(all, ad) } 
}

if (file.size(ar_file) == 0) {
  ar <- data.frame("family_id" = family_name, "note" = "Empty ar query")
} else {
  ar <- sortFilterGemini(ar_file) 
  if (nrow(ar) > 0) { ar <- filter(ar, !chrom %in% c("X", "chrX"))
    all <- rbind(all, ar) } 
}

if (file.size(comphets_file) == 0) {
  comphets <- data.frame("family_id" = family_name, "note" = "Empty comphets query")
} else {
  comphets <- sortFilterGemini(comphets_file) 
  if (nrow(comphets) > 0) { comphets <- filter(comphets, !is.na(gene)) %>% 
    distinct(chr_variant_id, .keep_all = TRUE) %>% 
    select(-comp_het_id, -priority)
    all <- rbind(all, comphets) } 
}

if (file.size(xdenovo_file) == 0) {
  xdenovo <- data.frame("family_id" = family_name, "note" = "Empty xdenovo query")
} else {
  xdenovo <- sortFilterGemini(xdenovo_file)
  if (nrow(xdenovo) > 0) { all <- rbind(all, xdenovo) } 
}

if (file.size(xd_file) == 0) {
  xd <- data.frame("family_id" = family_name, "note" = "Empty xd query")
} else {
  xd <- sortFilterGemini(xd_file)
  if (nrow(xd) > 0) { all <- rbind(all, xd) } 
}

if (file.size(xr_file) == 0) {
  xr <- data.frame("family_id" = family_name, "note" = "Empty xr query")
} else {
  xr <- sortFilterGemini(xr_file)
  if (nrow(xr) > 0) { all <- rbind(all, xr) } 
}

if (file.size(mendel_errors_file) == 0) {
  mendel_errors <- data.frame("family_id" = family_name, "note" = "Empty mendel_errors query")
} else {
  mendel_errors <- sortFilterGemini(mendel_errors_file)
  if (nrow(mendel_errors) > 0) { 
    mendel_errors <- mendel_errors %>% select(-violation)
    all <- rbind(all, mendel_errors) 
    } 
}

#summaryInfo <- data.frame("family_id" = family_name, "DxOutcome"= NA, "variant" = NA, "reviewer" = NA, "date" = NA, "SecondReviewer" = NA, "SecondReviewDate" = NA)
summaryInfo <- data.frame("family_id" = family_name, "PatientDxPhenotype" = NA, "DxOutcome"= NA, "variant" = NA, "reviewer" = NA, "date" = NA) %>% 
  add_row("family_id" = family_name, "PatientDxPhenotype" = NA, "DxOutcome"= NA, "variant" = NA, "reviewer" = NA, "date" = NA)

acmg_genes <- read_xlsx(geneCategory_file, sheet = "ACMG", na = c("NA", "", "None", "NONE", ".")) %>% pull(Gene) %>% unique()
if (dim(all)[1] == 0) {
  acmg <- data.frame("family_id" = family_name, "note" = "Empty ACMG")
} else {
  acmg <- all %>% filter(ref_gene %in% acmg_genes, priority_score > 4) %>% distinct(chr_variant_id, .keep_all = TRUE)
}

#in the next version, consider TTN truncating, HFE C..Y hmz etc.
config <- read_tsv(config_file, col_names = FALSE, na = c("NA", ""), col_types = cols(.default = col_character())) %>% 
  separate("X1", c("tool", "version", "note"), sep = "\\:|\\#", remove = TRUE)

openxlsx::write.xlsx(list( "de_novo" = denovo, 
                           "AD" = ad, 
                           "AR" = ar, 
                           "comp_hets" = comphets, 
                           "X_denovo" = xdenovo, 
                           "XD" = xd, 
                           "XR" = xr, 
                           "mendel_errors" = mendel_errors, 
                           "ACMG" = acmg,
                           "config" = config,
                           "summary" = summaryInfo), 
                     file = output_xlsx, firstRow = TRUE, firstCol = TRUE)
