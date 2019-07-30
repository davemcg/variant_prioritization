#updated 7/28/19
args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
#args <- c("freebayes.filtered.20190722.miseq0626_0719.G05139.gemini.test.tsv", "G05139.sorted.tsv", "OGLv1_panel_DxORcandidate.txt")

library(tidyverse)

#gemini <- read_tsv(args[1], col_names = TRUE, na = c("NA", "", "None"), col_types = cols (.default = "c", time = "i") )
# Read in gemini_tsv file.
gemini <- read.delim(args[1], sep = "\t", header = T,  na.strings=c("", "NA", "None"), check.names=FALSE, colClasses = c("factor","integer","integer","character","character","numeric","factor","integer","character","character",
                                    "character","character","factor","character","character","factor","character","character","character","character",
                                    "character","character","character","character","factor","numeric","numeric","numeric","numeric","factor",
                                    "numeric","numeric","numeric","numeric","numeric","numeric","character","numeric","character","character",
                                    "numeric","character","character","character","numeric","numeric","character","numeric","numeric","numeric",
                                    "numeric","numeric","numeric","factor","factor","character","character","character","character","character",
                                    "character","character","numeric","numeric","numeric","numeric","character","character","character","character",
                                    "character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","factor","integer",
                                    "integer","character","character","character","character","numeric","numeric","numeric","numeric","numeric",
                                    "character","factor","character","factor","integer","numeric","numeric")) %>%
   mutate(exon = sub("^", " ", exon))



#Read in OGLv1_panel gene class of either Dx or Candidate

#sampleData <- read.tsv(args[1], header = TRUE, check.names=FALSE ) #nrows = 5
# classes <- sapply(gemini, class)
# largeData <- read.csv("huge-file.csv", header = TRUE, colClasses = classes)


OGLv1_gene_class <- read.delim(args[3], sep = "\t", header = T, colClasses = c("character","character") )

# get max_priority_score for each gene
max_priority_score <- select(gemini, c(ref_gene_annovar, priority_score_intervar)) %>% group_by(ref_gene_annovar) %>% summarize(maxpriorityscore = max(priority_score_intervar)) 

#arrange by max_priority_score, then by gene, and priority score. None gene region?
gemini_max_priority_score <- left_join(gemini, max_priority_score, by=c("ref_gene_annovar"))

#VEP hg19 version's gene names are the same as in the IDT ordering design sheets. This is what used for left_join
gemini_sorted <- left_join(gemini_max_priority_score, OGLv1_gene_class, by=c("gene")) %>% 
  filter(priority_score_intervar > 2) %>% 
  arrange(desc(maxpriorityscore), ref_gene_annovar, desc(priority_score_intervar)) %>% 
  select(-maxpriorityscore) %>% 
  select('chr_annovar', 'start_annovar', 'ref_annovar', 'alt_annovar', 'end_annovar', 'qual', 'filter', starts_with('gts'), starts_with('gt_'), 'panel_class', 'priority_score_intervar', everything() )

write_tsv(gemini_sorted, file.path('.', args[2]))
