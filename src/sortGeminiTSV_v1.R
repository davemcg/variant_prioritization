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

# args <- c("W:/ddl_nisc_custom_capture/042020/freebayesPrioritization/gemini_tsv/108976P.20200420.freebayes.nisc100.gemini.tsv",
#           "Z:/resources/OGLpanelGeneDxORcandidate.xlsx", "rearranged.tsv", "filtered.tsv", "108976P", "filtered.xlsx", "0.5", "W:/ddl_nisc_custom_capture/042020/CoNVaDING/CNV_hiSens/108976P.b37.aligned.only.best.score.shortlist.txt")

gemini_file <- args[1]
geneCatogory_file <- args[2]
rearrangedGemini_file <- args[3]
filteredGemini_tsv_file <- args[4]
sampleName <- args[5]
gemini_xlsx_file <- args[6]
aafCutoff <- args[7]
manta_file <- args[8]
scramble_mei_file <- args[9]
scramble_del_file <- args[10]
convading_file <- args[11]
convading_LAF <- args[12]

library(tidyverse)
library(readxl)
library(RColorBrewer)

gemini_input <- read_tsv(gemini_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  type_convert() %>% mutate(exon = sub("^", " ", exon))
print("###gemini tsv loaded### 10%")
gemini <-  gemini_input %>% mutate( start_vcf = start + 1 ) %>% 
  unite("chr_variant_id", chrom, start_vcf, ref, alt, sep = "-", remove = FALSE ) %>% 
  mutate(gene = toupper(gene)) %>% 
  mutate(sample = sampleName) 
rm(gemini_input)
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
panelGene <- read_xlsx(geneCatogory_file, sheet = "analysis", na = c("NA", "", "None", ".")) %>% select(gene, panel_class) %>% rename(ref_gene = gene)

# get max_priority_score for each gene
max_priority_score <- select(gemini, c(ref_gene, priority_score)) %>% group_by(ref_gene) %>% summarize(maxpriorityscore = max(priority_score)) 
print("###max priority score### 20%")
#arrange by max_priority_score, then by gene, and priority score. None gene region?
gemini_max_priority_score <- left_join(gemini, max_priority_score, by=c("ref_gene"))
rm(gemini)
#VEP hg19 version's gene names are the same as in the IDT ordering design sheets. This is what used for left_join
gemini_rearrangeCol <- left_join(gemini_max_priority_score, panelGene, by = c("ref_gene")) %>% 
  mutate(note = "") %>% separate(vcf_id, c('caller', 'hg38_id'), sep = "_") %>% 
  mutate(gnomad_hom = gno_hom + gnog_hom, gnomad_ac = gno_ac_all + gnog_ac_all, gnomad_an = gno_an_all + gno_an_all ) %>% 
  mutate(gnomad_af = gnomad_ac/gnomad_an) %>% 
  unite("gnomad_acan", gnomad_ac, gnomad_an, sep = "/", remove = TRUE) %>% 
  mutate(eyeGene = case_when(panel_class == "Dx" ~ 2,
                             panel_class == "Candidate" ~ 1,
                             TRUE ~ 0)) %>% 
  select('chr_variant_id', 'sample', 'chrom', 'start_vcf', 'qual', 'filter', starts_with('gts'), starts_with('gt_'), 'aaf', 'caller', 'hg38_id',
         'panel_class', 'priority_score', 'priority_score_intervar', 'clinvar_hgmd_score', 'splice_score', 'other_predic_score', 'gnomad_af', 'gnomad_acan', 'pmaxaf', 'max_af', 'max_af_pops', 'gnomad_hom', 'ref_gene', 'note', 
         'exonicfunc_ensgene', 'refgenewithver', 'hgvsc', 'hgvsp', 'gene', 'exon', 'aa_length', 'omim_gene', 'omim_inheritance', 'omim_phenotype', 'pvs1', 'truncating_vep', 'hgmd_overlap', 'existing_variation', 'clinvar_intervar', 'intervar_and_evidence', 
         'clinvar_id', 'clinvar_pathogenic', 'clinvar_sig', 'clin_sig', 'existing_inframe_oorfs','existing_outofframe_oorfs','existing_uorfs','five_prime_utr_variant_annotation','five_prime_utr_variant_consequence',
         'spliceai', 'spliceai_maxscore', 'spliceaimasked50', 'spliceaimasked50max', 'dbscsnv_ada_score_intervar', 'dbscsnv_rf_score_intervar', 'dpsi_max_tissue_annovar', 'dpsi_zscore_annovar', 'genesplicer', 'maxentscan_diff', 'branchpoint_u2_binding_energy', 'branchpoint_prob', 'branchpoint_to_3prime', 
         'sift_pred', 'polyphen_pred', 'mutationassessor_pred', 'mutationtaster_pred', 'metasvm_pred','metasvm_score_intervar', 'clinpred_score', 'primatedl', 'revel', 'ccr_pct','mpc', 'cadd_raw', 'cadd_phred','remm', 'fathmm_xf_coding','fathmm_xf_noncoding','eigen_pc_raw', 'eigen_raw', 'gerp_rs_intervar', 'phylop46way_placental_intervar', 'phylop_100way', 
         'atac_rpe_score','atac_rpe_itemrgb', 'ft_ret_rpe_score', 'genedetail_ensgene', 'aachange_ensgene', 'gene_refgenewithver', 'func_refgene', 'func_refgenewithver', 'exonicfunc_refgenewithver', 'exonicfunc_refgene', 'avsnp150_annovar', 'interpro_domain_intervar', 
         'pfam_domain', 'tfbs', 'pli', 'lof_z', 'mis_z', 'sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01', 'atac_rpe_itemrgb', 'atac_rpe_score', 
         'eyeintegration_rnaseq_tpm_rpe_adulttissue', 'eyeintegration_rnaseq_tpm_rpe_cellline', 'eyeintegration_rnaseq_tpm_rpe_fetaltissue', 'eyeintegration_rnaseq_tpm_rpe_stemcellline', 'eyeintegration_rnaseq_tpm_retina_adulttissue', 'eyeintegration_rnaseq_tpm_retina_stemcellline', 'eyeintegration_rnaseq_tpm_wholeblood', 
         'start', 'end', 'ref', 'alt', 'popfreqmax_annovar', 'gnomad_exome_all_annovar', 'gnomad_genome_all_annovar', 'freq_esp6500siv2_all_annovar', 'freq_1000g2015aug_all_annovar', 'aaf_esp_all', 'aaf_1kg_all', 'pubmed', everything() )
#4/12/20: removed 'chr_annovar', 'start_annovar', 'ref_annovar', 'alt_annovar',
write_tsv(gemini_rearrangeCol, file = rearrangedGemini_file)
print("###rearranged file written### 30%")
gemini_filtered <- gemini_rearrangeCol %>% mutate(temp_group = ifelse(priority_score >= 3, 3, ifelse(priority_score >= -2, -2, -3))) %>% 
  filter(temp_group >= -2, pmaxaf < 0.2, aaf < aafCutoff) %>% arrange(desc(eyeGene), desc(temp_group), desc(maxpriorityscore), ref_gene, desc(priority_score)) %>% 
  select(-temp_group)
gemini_filtered0 <- gemini_filtered %>% select(-maxpriorityscore)
# consider change to filter(priority_score > 10 | (temp_group >= -2, pmaxaf < 0.05, aaf < aafCutoff))

#AnnotSV does not have GD_POPMAX_AF column as of 3/2/2021, if GD_POPMAX_AF in columnames, then use it in the next version.
write_tsv(gemini_filtered0, file = filteredGemini_tsv_file)
#found out whether blank is 0 after importing, use 1197 for filtering
#consider adding manta to the main gemini df for sorting/filtering after knowing the specificity of the manta calls. To better sort AR, AD, and ACMG 2nd.

if (file.size(manta_file) == 0) {
  manta_sort <- data.frame("sample" = sampleName, "note" = "Empty manta")
} else {
  manta <- read_tsv(manta_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
    filter(FILTER == 'PASS') %>%
    filter( is.na(Gene_name) | Annotation_mode == 'split' & !is.na(Gene_name) ) %>% 
    mutate(ACMG_class = sub("full=", "", ACMG_class)) %>% 
    mutate(ACMG_class = as.integer(ACMG_class), SV_length = as.numeric(SV_length)) %>% 
    filter(ACMG_class > 1, is.na(SV_length) | abs(SV_length) < 1000000) %>% 
    separate(Location, c('temp_location1', 'temp_location2'), sep = "-", remove = FALSE, convert = FALSE) %>% 
    filter(!(grepl("intron", Location) & temp_location1 == temp_location2)) %>% 
    select(-starts_with('temp_')) %>% 
    rename(ref_gene = "Gene_name")
  manta_sort <- left_join(manta, panelGene, by = c("ref_gene")) %>% 
    mutate(note = "") %>% 
    mutate(eyeGene = case_when(panel_class == "Dx" ~ 2,
                               panel_class == "Candidate" ~ 1,
                               TRUE ~ 0)) %>% 
    arrange(desc(eyeGene)) %>% 
    select(AnnotSV_ID:Annotation_mode,panel_class,note,everything() )
}

if (file.size(scramble_mei_file) == 0) {
  scramble_mei <- data.frame("sample" = sampleName, "note" = "Empty scramble mei")
} else {
  scramble_mei <- read_xlsx(scramble_mei_file, na = c("NA", "", "None", "."))
}

if (!is.na(scramble_del_file)) {
  if (file.size(scramble_del_file) == 0) {
    scramble_del_sort <- data.frame("sample" = sampleName, "note" = "Empty scramble del") } 
  else {
    scramble_del <- read_tsv(scramble_del_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
      filter( is.na(Gene_name) | Annotation_mode == 'split' & !is.na(Gene_name) ) %>% 
      mutate(ACMG_class = sub("full=", "", ACMG_class)) %>% 
      mutate(ACMG_class = as.integer(ACMG_class), SV_length = as.numeric(SV_length)) %>% 
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
      arrange(desc(eyeGene)) %>% 
      select(AnnotSV_ID:Annotation_mode,panel_class,note,everything() ) 
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

xR <- gemini_filtered3 %>% filter(chrom == "X", recessive_cnt >= 2) %>% select(-maxpriorityscore, -recessive_cnt)
xD <- gemini_filtered3 %>% filter(chrom == "X", recessive_cnt == 1, pmaxaf < 0.002) %>% select(-maxpriorityscore, -recessive_cnt)

#ar are those genes with homozygous or compound hets variants of ps >= 5. However, ps = 4 variants were also listed if there are 2 ps>=5.
#ar gene with 1 hit will not be here.
ar <- gemini_filtered3 %>% filter(!chrom %in% c("X", "Y"), !omim_inheritance %in% c("AD"), recessive_cnt >= 2) %>% 
  mutate(knownAR = ifelse(grepl("AR", omim_inheritance), 1, 0)) %>% 
  arrange(desc(eyeGene), desc(knownAR), desc(maxpriorityscore), ref_gene, desc(priority_score)) %>% 
  mutate(note = ifelse(ref_gene %in% c("PRPH2", "ROM1", "PCDH15", "CDH23", "CNGA1", "CNGB1", "CNGA3", "CNGB3"), "Digenic?", note) ) %>%
  select(-maxpriorityscore, -knownAR, -recessive_cnt) # digenic recessive

#ad are those omim unknown or AD inheritance. 
ad <- gemini_filtered3 %>% filter(!chrom %in% c("X", "Y"), 
                                  recessive_cnt == 1 & !omim_inheritance %in% c("AR") | recessive_cnt >=2 & grepl("AD", omim_inheritance),
                                  pmaxaf < 0.002, priority_score >= 5) %>% 
  arrange(desc(eyeGene), desc(maxpriorityscore), ref_gene, desc(priority_score)) %>% 
  select(-maxpriorityscore, -recessive_cnt)
print("###inheritance search done### 50%")
#grepl("AD", omim_inheritance) | is.na(omim_inheritance)

#AD: score > 4 , AR: score > 4, all: score >= 3

acmg_genes = c('ACTA2','ACTC1','APC','APOB','ATP7B','BMPR1A','BRCA1','BRCA2',
               'CACNA1S','COL3A1','DSC2','DSG2','DSP','FBN1','GLA','KCNH2','KCNQ1',
               'LDLR','LMNA','MEN1','MLH1','MSH2','MSH6','MUTYH','MYBPC3','MYH11',
               'MYH7','MYL2','MYL3','NF2','OTC','PCSK9','PKP2','PMS2','PRKAG2',
               'PTEN','RB1','RET','RYR1','RYR2','SCN5A','SDHAF2','SDHB','SDHC',
               'SDHD','SMAD3','SMAD4','STK11','TGFBR1','TGFBR2','TMEM43','TNNI3',
               'TNNT2','TP53','TPM1','TSC1','TSC2','VHL','WT1')
acmg <- gemini_filtered3 %>% filter(ref_gene %in% acmg_genes, priority_score > 4) %>% select(-maxpriorityscore, -recessive_cnt)
print("###acmg done### 70%")
summaryInfo <- data.frame("sample" = sampleName, "DxOutcome"= NA, "variant" = NA, "reviewer" = NA, "date" = NA, "2ndReviewer" = NA, "2ndReviewDate" = NA)
if (is.na(scramble_del_file)) {
  openxlsx::write.xlsx(list("AR" = ar, "AD" = ad, "XR" = xR, "XD" = xD, "ACMG59" = acmg, "all" = gemini_filtered1, "manta" = manta_sort, "scramble_mei" = scramble_mei, "summary" = summaryInfo), file = gemini_xlsx_file, firstRow = TRUE, firstCol = TRUE)
} else if (is.na(convading_file)) {
  openxlsx::write.xlsx(list("AR" = ar, "AD" = ad, "XR" = xR, "XD" = xD, "ACMG59" = acmg, "all" = gemini_filtered1, "manta" = manta_sort, "scramble_mei" = scramble_mei, "scramble_del" = scramble_del_sort, "summary" = summaryInfo), file = gemini_xlsx_file, firstRow = TRUE, firstCol = TRUE)
} else {
    cnv <- read_tsv(convading_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
      type_convert() 
    if (dim(cnv)[1] == 0) {
      openxlsx::write.xlsx(list("AR" = ar, "AD" = ad, "XR" = xR, "XD" = xD, "ACMG59" = acmg, "all" = gemini_filtered1, "manta" = manta_sort, "scramble_mei" = scramble_mei, "scramble_del" = scramble_del_sort, "summary" = summaryInfo), file = gemini_xlsx_file, firstRow = TRUE, firstCol = TRUE)
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
      openxlsx::write.xlsx(list("AR" = ar, "AD" = ad, "XR" = xR, "XD" = xD, "ACMG59" = acmg, "all" = gemini_filtered1, "CoNVaDING" = cnv, "CNV-variant" = cnv_variant, "manta" = manta_sort, "scramble_mei" = scramble_mei, "scramble_del" = scramble_del_sort, "summary" = summaryInfo), file = gemini_xlsx_file, firstRow = TRUE, firstCol = TRUE)
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

