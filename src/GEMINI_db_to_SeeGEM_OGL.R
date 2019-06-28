# Generate a full series of GEMINI queries for a sample/trio
args = commandArgs(trailingOnly=TRUE)

# argument handling
gemini_db <- args[1]
family_name <- args[2]
output_html <- args[3]
peddy_path <- args[4]
# cutoff cohort AF
if (is.null(args[5])) {
	aaf_freq <- 0.1 } else {
	aaf_freq <- args[5]
}
# run gemini with lenient handing?
if (is.null(args[6])) {
	lenient <- ''
} else {
	lenient <- args[6]
}
if (toupper(lenient) != 'YES') {
	lenient <- ''
} else {
	lenient <- '--lenient'
}
# output GEMINI as data frame also?
if (is.na(args[7])) {
    output_df <- '' } else {
    output_df <- args[7]
}

cur_dir <- getwd()

library(SeeGEM)
library(tidyr)

writeLines('\n\n\n\n\n##########################################################')
writeLines('Starting GEMINI queries')

GEMINI_list <- list()
GEMINI_list$ar <- gemini_test_wrapper(gemini_db, 
                                      test = 'autosomal_recessive', 
                                      min_gq = 0,
                                      "--allow-unaffected",
                                      filter = paste("aaf < ", aaf_freq, " AND priority_score_intervar > -4"),
                                      families = family_name, ... = lenient)
writeLines('Autosomal Recessive test done')
GEMINI_list$ad <- gemini_test_wrapper(gemini_db, 
                                      test = 'autosomal_dominant', 
                                      min_gq = 0,
                                      "--allow-unaffected",
                                      filter = paste("aaf < ", aaf_freq, " AND priority_score_intervar > -4"),
                                      families = family_name, ... = lenient)
writeLines('Autosomal Dominant test done')
GEMINI_list$dn <- gemini_test_wrapper(gemini_db, 
                                      test = 'de_novo',
                                      min_gq = 0,
                                      "--allow-unaffected",
                                      filter = paste("aaf < ", aaf_freq, " AND priority_score_intervar > -4"),
                                      families = family_name, ... = lenient)
writeLines('De novo test done')
GEMINI_list$xlr <- gemini_test_wrapper(gemini_db, 
                                       test = 'x_linked_recessive',
                                       min_gq = 0,
                                       "--allow-unaffected",
                                       filter = paste("aaf < ", aaf_freq, " AND priority_score_intervar > -4"),
                                       families = family_name)
writeLines('XL Recessive test done')
GEMINI_list$xld <- gemini_test_wrapper(gemini_db, 
                                       test = 'x_linked_dominant',
                                       min_gq = 0,
                                       "--allow-unaffected", 
                                       filter = paste("aaf < ", aaf_freq, " AND priority_score_intervar > -4"),
                                       families = family_name)
writeLines('XL Dominant test done')
GEMINI_list$xldn <- gemini_test_wrapper(gemini_db, 
                                        test = 'x_linked_de_novo',
                                        min_gq = 0,
                                        "--allow-unaffected", 
                                        filter = paste("aaf < ", aaf_freq, " AND priority_score_intervar > -4"),
                                        families = family_name)
writeLines('XL De Novo test done')
GEMINI_list$me <- gemini_test_wrapper(gemini_db, 
                                      test = 'mendel_errors', 
                                      min_gq = 0,
#                                      "--allow-unaffected",
                                      filter = paste("aaf < ", aaf_freq, " AND priority_score_intervar > -4"),
                                      families = family_name, ... = lenient)
writeLines('Mendelian Errors test done')
GEMINI_list$ch <- gemini_test_wrapper(gemini_db, 
                                      test = 'comp_hets', 
                                      min_gq = 0,
                                      "--allow-unaffected", 
                                      filter = paste("aaf < ", aaf_freq, " AND priority_score_intervar > -4"),
                                      families = family_name)
writeLines('Compound Hets test done')

acmg_genes = c('ACTA2','ACTC1','APC','APOB','ATP7B','BMPR1A','BRCA1','BRCA2',
               'CACNA1S','COL3A1','DSC2','DSG2','DSP','FBN1','GLA','KCNH2','KCNQ1',
               'LDLR','LMNA','MEN1','MLH1','MSH2','MSH6','MUTYH','MYBPC3','MYH11',
               'MYH7','MYL2','MYL3','NF2','OTC','PCSK9','PKP2','PMS2','PRKAG2',
               'PTEN','RB1','RET','RYR1','RYR2','SCN5A','SDHAF2','SDHB','SDHC',
               'SDHD','SMAD3','SMAD4','STK11','TGFBR1','TGFBR2','TMEM43','TNNI3',
               'TNNT2','TP53','TPM1','TSC1','TSC2','VHL','WT1')
# gemini select * is annoying in that it won't grab the sample genotypes
# they have to be explicitly given
# so we'll pull the samples names for the family
sample_ped <- gemini_query_wrapper(gemini_db,
                                   ... = paste0("\"SELECT * FROM samples WHERE family_id == '",
                                                family_name, "' \""))
gts = paste0(rep('gts.', length(sample_ped$name)), sample_ped$name, collapse = ',')
gts_var = paste0(rep('gts.', length(sample_ped$name)), sample_ped$name)
gt_types_het = paste0(rep('gt_types.', length(sample_ped$name)), sample_ped$name, rep(' == HET', length(sample_ped$name)), collapse = ' or ')
gt_types_hom_alt = paste0(rep('gt_types.', length(sample_ped$name)), sample_ped$name, rep(' == HOM_ALT', length(sample_ped$name)), collapse = ' or ')
GEMINI_list$acmg <- gemini_query_wrapper(gemini_db,
                             ... = paste0("\"SELECT *,", gts, " FROM variants WHERE (gene IN (\'",
                                          paste(acmg_genes, collapse="\',\'"),
                                          "\')) AND ((clinvar_sig LIKE '%pathogenic%' OR impact_severity='HIGH')
                                          AND (aaf <= ", aaf_freq, " AND aaf_esp_all < 0.01 AND
                                          aaf_1kg_all < 0.01 AND af_exac_all < 0.01))
                                          AND filter IS NULL \" --gt-filter \"", gt_types_het, " or ", gt_types_hom_alt, "\""),
                                          test_name = 'ACMG59')
                                          
# make the family genotypes column
GEMINI_list$acmg$family_members <- paste(sample_ped$name, collapse = ",")
GEMINI_list$acmg <- unite(GEMINI_list$acmg, family_genotypes, gts_var, sep = ",")
writeLines('ACMG test done')

# data.table rbindlist will collapse each element of the list into one data frame
# gemini_query_wrapper() and gemini_test_wrapper() will add the test name
# to each query, so you can distinguish them later (via the 'test' column)
my_GEMINI_data <- data.table::rbindlist(GEMINI_list, fill = TRUE)
if (output_df != '' & toupper(stringr::str_sub(output_df, -5,-1)) == 'RDATA' ){
	save(my_GEMINI_data, file = output_df)
} else if (output_df != '' & toupper(stringr::str_sub(output_df, -3,-1)) == 'TSV' ){
	readr::write_tsv(my_GEMINI_data, path = output_df)
} else {}

# now that you've created the core data, you can create the reactive document
# I'm assuming you've already run peddy on the same vcf you used to make the GEMINI
# db. 

# one wrinkle is that peddy doesn't give the family labels throughout the output,
# rather it uses the sample ids. So we need to get then from the `sample_ped` query above
writeLines('Create reactive document!')
# first decorate
decorated <- See_GEM_formatter(my_GEMINI_data,
								extra_columns_to_retain = "^gno|rankscore$|*num*|^clin|*domain*|*codon*|*annov*|*interv*|spliceai*|ClinPred_Score|branchpoint*|omim*|atac*|PrimateDL|GeneSplicer|MaxEntScan*|MPC|Existing_variation")
if (toupper(peddy_path) != 'NO'){
	knit_see_gem(GEMINI_data = decorated, 
    	         output_file = paste0(cur_dir, '/', output_html), 
        	     peddy_path_prefix = paste0(cur_dir,'/',peddy_path), 
            	 peddy_id = sample_ped$name, 
            	 sample_name = family_name,
				 decorate = FALSE)
} else {
	knit_see_gem(GEMINI_data = decorated, 
    	         output_file = paste0(cur_dir, '/', output_html), 
            	 skip_stats = 'yes', 
				 decorate = FALSE,
				 sample_name = family_name)
}

