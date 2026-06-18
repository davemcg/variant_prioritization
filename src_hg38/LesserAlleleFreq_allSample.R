
## #!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("Please input the sample name and directory of output plot", call.=FALSE)
# } else if (length(args)==1) {
#   # default output file
#   args[2] = "out.txt"
# }

#I actually use ‘fread’ from the data.tables package for import, as it will switch the column type to chr if it detects a problem like this.
#But you can fix this issue with readr by manually setting the type. Problem is (I believe) you either set nothing and let it choose or you have to tell it what everything is. But there is a way to set a ‘default’ import and then change a few things like so:
#read_tsv(df, col_types = cols(.default = "i", chr = "c"))
#c = character, i = integer, n = number, d = double, l = logical, D = date, T = date time, t = time, ? = guess, or _/- to skip the column.

#library(bedr)
library(tidyverse)
library(ggsci)
rm(list=ls(all=TRUE))
# it seems that R sort the chr correctly, no need to sort in Terminal, but the following input file was sorted already.
vcf100 <- read_tsv('/BG/Panel/VCF/DDL_NISC_targeted_panel.hardFilterSNP-INDEL.vcf.PASS.noOPN.txt',
                   col_types = cols(.default = "c", POS = "i", QUAL = "d")) %>% mutate(variant_no = row_number())

long_vcf<- gather(vcf100, 
                         SAMPLE_NAME, GT_RD,
                         '107768M':'109862J')

long_vcf_sep <- long_vcf %>% separate(GT_RD, c('GT','AD','DP','GQ','PL'), sep = ':') %>%
  separate(AD, c('AD_REF', 'AD_ALT'), sep = ',') #%>%
#  filter(str_detect(GT, '0/1')) 
#yyz$b <- as.numeric(as.character(yyz$b))
long_vcf_sep$AD_REF <- as.integer(as.character(long_vcf_sep$AD_REF))
long_vcf_sep$AD_ALT <- as.integer(as.character(long_vcf_sep$AD_ALT))
long_vcf_sep$DP <- as.integer(as.character(long_vcf_sep$DP))


LessAlleleFreq <- long_vcf_sep %>% mutate(AD_LessAllele = ifelse(AD_REF <= AD_ALT, AD_REF, AD_ALT)) %>%
  mutate(LessAlleleFreq = AD_LessAllele/(AD_REF + AD_ALT))

# write_tsv(LessAlleleFreq_S1, file.path('/BG/MyR_BG/NISC_Panel', "LessAlleleFreq_108976P.txt"))

LessAlleleFreq_DP200 <- LessAlleleFreq %>% filter(DP>199, (AD_REF + AD_ALT)>199)

# The following is not necessary
# LessAlleleFreq_S1_DP200_1 <- LessAlleleFreq_S1_DP200 %>% mutate(VariantNo = 1:n()) %>% select(VariantNo, everything())
# facet_wrap(~CHROM, scales = 'free_x') +
ggplot(LessAlleleFreq_DP200, aes(x= variant_no, y = LessAlleleFreq, color = CHROM)) + 
  scale_color_discrete (breaks=c("1","2","3","4","5","6","7","8","9","10", "11", "12","13","14","15","16","17","19","20", "21","22","X")) +
  coord_cartesian(ylim = c(0.05, 0.52)) +
  labs(x= 'Het variants', y = 'Less-allele frequency') +
  facet_wrap(~SAMPLE_NAME) +
  geom_point(alpha = 0.5, size = 1) +
  theme_bw() +
  theme(axis.text.x  = element_text(size=8), axis.text.y  = element_text(size=8)) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16)) # +
#  scale_color_aaas()

#ggplot(LessAlleleFreq_S1, aes(x= 1:nrow(LessAlleleFreq_S1), y = LessAlleleFreq), alpha = 0.2) + 
#  geom_point()

#vcf <- read.vcf('the_vcf')
#vcf$vcf %>% head()
# this will work for a single sample VCF
# you will have to tweak to match the 'FORMAT' fields in your vcf
# in this example there are GT, AF, DP, and AD fields
# 10 refers to the 10th column, which would be the first sample in your vcf
#vcf$vcf %>% separate(10, c('GT','AF','DP','AD'), sep = ':') %>% head()
# you can then grab the DP field, which is what you need
#DP_values <- vcf$vcf %>% separate(10, c('GT','AF','DP','AD'), sep = ':') %>% select(CHROM, POS, REF, ALT, DP)
