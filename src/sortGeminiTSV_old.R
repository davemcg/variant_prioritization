
args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
#args <- c("1197.tsv", "1197.sorted.tsv", "OGLv1_panel_DxORcandidate.txt")

library(tidyverse)

# Read in gemini_tsv file.
gemini <- read.delim(args[1], sep = "\t", header = T,  na.strings=c("NA", "None"), check.names=FALSE, 
                     colClasses = c("factor","integer","integer","character","character","numeric","factor","factor","integer","integer",
                                    "character","character","integer","character","character","character","character","character","character","character",
                                    "character","character","character","character","character","character","character","character","character","character",
                                    "character","character","character","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                    "numeric","numeric","character","numeric","numeric","numeric","numeric","numeric","numeric","character",
                                    "numeric","character","character","numeric","character","character","character","numeric","numeric","character",
                                    "numeric","numeric","numeric","numeric","numeric","numeric","character","character","character","character",
                                    "character","character","character","character","character","numeric","numeric","numeric","numeric","character",
                                    "character","character","character","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                    "numeric","integer","character","integer","integer","numeric","numeric") ) %>% 
  mutate(exon = sub("^", " ", exon))

# Read in OGLv1_panel gene class of either Dx or Candidate

OGLv1_gene_class <- read.delim(args[3], sep = "\t", header = T, colClasses = c("character","character") )

# get max_priority_score for each gene
max_priority_score <- select(gemini, c(ref_gene_annovar, priority_score_intervar)) %>% group_by(ref_gene_annovar) %>% summarize(maxpriorityscore = max(priority_score_intervar)) 

#arrange by max_priority_score, then by gene, and priority score. None gene region?
gemini_max_priority_score <- left_join(gemini, max_priority_score, by=c("ref_gene_annovar"))

#VEP hg19 version's gene names are the same as in the IDT ordering design sheets. This is what used for left_join
gemini_sorted <- left_join(gemini_max_priority_score, OGLv1_gene_class, by=c("gene")) %>% 
  arrange(desc(maxpriorityscore), ref_gene_annovar, desc(priority_score_intervar)) %>% 
  select(-maxpriorityscore) %>% 
  select('chrom', 'start', 'end', 'ref', 'alt', 'qual', 'filter', 'chr_annovar', 'start_annovar', 'ref_annovar', 'alt_annovar', 'end_annovar', starts_with('gts'), starts_with('gt_'), 'panel_class', everything() )

write_tsv(gemini_sorted, file.path('.', args[2]))
