### sort intervar output according to PVS/PS/PM/PP ###
### after InterVar and Annovar to prepare for vcfanno
## PVS = 8
## PS  = 6
## PM  = 3
## PP  = 1
## BA = -5
## BS = -3
## BP = -1
## if PVS = 0 & dbscSNV's ada > 0.97 and rf>0.63, Score=+3 (equals weight of PM)
## if PVS = 0 & frameshit or stop_gain or Start_loss (?intervar), Score=+3 (equals weight of PM)
## Intervar - select one gene for each variant

args <- commandArgs(trailingOnly=TRUE)

library(tidyverse)

# change the file name below
intervar <- read.delim(args[1], sep = "\t", header = TRUE, na.strings = c("."),
                       colClasses = c("factor","integer","integer","character","character","character","character","character","character",
										"character","character","character","character","character","numeric","numeric","numeric","numeric",
										"numeric","numeric","numeric","numeric","numeric","numeric","character","character","character","numeric",
										"character","integer","character","character","character","character") )

#does the read.delim evaluate all rows for the types? seems to be yes.
annovar <- read.delim(args[2], sep = "\t", header = TRUE, na.strings = c("."),
                      colClasses = c("factor","integer","integer","character","character","character","character","character","character","character",
                                     "numeric","numeric","character","numeric","numeric","character","numeric","numeric","character","numeric",
                                     "numeric","character","numeric","numeric","character","numeric","numeric","character","numeric","numeric",
                                     "character","numeric","numeric","character","numeric","numeric","character","numeric","numeric","numeric",
                                     "numeric","character","numeric","numeric","character","numeric","numeric","character","numeric","numeric",
                                     "numeric","numeric","numeric","numeric","numeric","character","character","numeric","numeric","numeric",
                                     "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                     "numeric","numeric","numeric","numeric","numeric","numeric","character","character","character","character",
                                     "character","character","character","character","numeric","numeric","numeric","numeric","numeric","numeric",
                                     "numeric","numeric","numeric","numeric","character","character","character","character","character","character",
                                     "character","character","character","character","character","character","character","character","character","character",
                                     "character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                     "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                     "numeric","character","factor","integer","character","character","character","numeric","factor") )


#pick one annotation for each variant that is of highest Priority score
intervar_for_sorting <- intervar %>% 
  separate('InterVar..InterVar.and.Evidence', 
           c('InterVarInterpretation','PVS1','PS','PM','PP','BA1','BS','BP'), 
           sep = '\\]{0,1} [A-Z]{2,3}\\d{0,1}\\=\\[{0,1}', remove = FALSE, convert = TRUE) %>% 
  separate(PS, c('PS1','PS2','PS3','PS4','PS5'), sep = ',', convert = TRUE) %>% 
  separate(PM, c('PM1','PM2','PM3','PM4','PM5','PM6','PM7'), sep =',', convert = TRUE) %>% 
  separate(PP, c('PP1','PP2','PP3','PP4','PP5','PP6'), sep = ',', convert = TRUE) %>% 
  separate(BS, c('BS1','BS2','BS3','BS4','BS5'), sep = ',', convert = TRUE) %>% 
  separate(BP, c('BP1','BP2','BP3','BP4','BP5','BP6','BP7',"BP8"), sep = ',', convert = TRUE) %>% 
  mutate(BP8 = str_sub(BP8, 2, 2)) %>% 
  mutate(BP8 = as.integer(BP8)) %>% 
  replace_na(list(dbscSNV_ADA_SCORE = 0, dbscSNV_RF_SCORE = 0)) %>%
  mutate(dbscSNV_value = ifelse((PVS1 == 0 & dbscSNV_ADA_SCORE>0.97 & dbscSNV_RF_SCORE>0.63), 1, 0)) %>% 
  mutate(truncating = ifelse((PVS1 == 0 & grepl("^frameshift|stop", ExonicFunc.refGene, ignore.case = TRUE)), 1, 0)) %>% 
  mutate(Priority.Score = (PVS1*8+(PS1+PS2+PS3+PS4+PS5)*6+(PM1+PM2+PM3+PM4+PM5+PM6+PM7+dbscSNV_value+truncating)*3+(PP1+PP2+PP3+PP4+PP5+PP6)-BA1*5-(BS1+BS2+BS3+BS4+BS5)*3-(BP1+BP2+BP3+BP4+BP5+BP6+BP7+BP8))) %>% 
  unite("variantkey", X.Chr:Alt, sep = "_", remove = FALSE ) %>%
  group_by(variantkey) %>%
  slice(which.max(Priority.Score))
  
annovar_inter <- merge(x = annovar, y = intervar_for_sorting, 
                       by.x = c("Chr", "Start", "End", "Ref", "Alt"), by.y = c("X.Chr", "Start", "End", "Ref", "Alt"), all.x = TRUE, 
                       sort = FALSE, suffixes = c(".annovar", ".intervar"), no.dups = TRUE,
                       incomparables = NULL) %>% 
  select(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "Chr", "Start", "Ref", "Alt", "End", "Priority.Score", "Ref.Gene",  
		"Func.refGene.intervar", "Gene.refGeneWithVer", "GeneDetail.refGeneWithVer", "ExonicFunc.refGeneWithVer", "AAChange.refGeneWithVer",
		"clinvar..Clinvar", "InterVar..InterVar.and.Evidence", "PopFreqMax", "gnomAD_exome_ALL", "gnomAD_genome_ALL", "Freq_ExAC_ALL", 
		"Freq_esp6500siv2_all", "Freq_1000g2015aug_all","dbscSNV_ADA_SCORE.intervar", "dbscSNV_RF_SCORE.intervar", "dpsi_max_tissue", "dpsi_zscore",
        "SIFT_score.intervar", "MetaSVM_score.intervar", "CADD_raw.intervar", "CADD_phred.intervar", "GERP.._RS.intervar", "phyloP46way_placental", 
		"Func.refGeneWithVer", "ExonicFunc.refGene.intervar", "avsnp150")) %>% 
  replace_na(list(ID = ".")) %>% 
  rename( Chr_annovar = Chr, Start_annovar = Start, Ref_annovar = Ref, Alt_annovar = Alt, End_annovar = End)

write_tsv(annovar_inter, file.path('.', args[3]))
