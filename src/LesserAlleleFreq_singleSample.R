
args <- commandArgs(trailingOnly=TRUE)
#When testing, comment out line above and use the line below.
#args <- c("W:/Brooks_Coloboma/rd2/prioritization_freebayes/gemini_tsv/freebayes.combined.filtered.rd2.COL124_1.gemini.sorted.tsv",
#          "COL124_1", "COL124_1.LAF.jpeg")

library(tidyverse)
library(RColorBrewer)

gemini_input <- read_tsv(args[1], col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  select('chr_annovar', 'start_annovar', starts_with('gts'), starts_with('gt_'), 
         'panel_class', 'gene_refgenewithver_annovar') %>% 
  type_convert() %>% 
  mutate(chr_annovar = as.factor(chr_annovar)) %>% 
  rename_all(funs(str_replace(., args[2], ""))) %>% 
  filter(gt_depths. >= 30) %>% 
  group_by(chr_annovar) %>%
  arrange(start_annovar, .by_group = TRUE) %>% 
  mutate(LesserAlleleFreq = ifelse(gt_alt_freqs. <= 0.5, gt_alt_freqs., 1 - gt_alt_freqs.)) %>% 
  mutate(variant_no = row_number()) %>% 
  mutate(DepthGroup = ifelse(gt_depths. >= 100, "DP>=100", "DP>=30" ) )

gemini_input$chr_annovar = factor(gemini_input$chr_annovar,
                                            levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")) 

gemini_input$DepthGroup = factor(gemini_input$DepthGroup,levels = c("DP>=100", "DP>=30"))

ggplot(gemini_input, aes(x= variant_no, y = LesserAlleleFreq, color = DepthGroup)) + 
  scale_color_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(0, 0.52)) +
  labs(title = args[2], x= 'Variants', y = 'Lesser allele frequency') +
  facet_wrap(~chr_annovar, scales = "free_x") +
  geom_point(alpha = 0.5, size = 1) +
  theme_bw() +
  theme(axis.text.x  = element_text(size=8), axis.text.y  = element_text(size=8)) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16)) +
  theme(legend.position = c(0.9, 0.1)) 

#theme(plot.title = element_text(hjust = 0.5))

ggsave(args[3], path = ".", width = 32, height = 18, units = "cm")

#write_tsv(gemini_filtered, file.path('.', args[4]))
