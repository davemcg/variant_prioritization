manta_file <- "Z:/genome/broad/manta/manta.D1440_01.annotated.tsv"
manta_freq_file <- "Z:/resources/manta/manta.OGL.freq.2022-09.tsv"
geneCategory_file <- "Z:/resources/OGLpanelGeneDxORcandidate.xlsx"

manta_freq <- read_tsv(manta_freq_file, col_names = TRUE, na = c("NA", "full=NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
  type_convert()
panelGene <- read_xlsx(geneCategory_file, sheet = "analysis", na = c("NA", "", "None", "NONE", ".")) %>%
  mutate(ref_gene = toupper(gene)) %>% 
  select(ref_gene, panel_class) %>% distinct()

manta_original <- read_tsv(manta_file, col_names = TRUE, na = c("NA", "full=NA", "", "None", "NONE", "."), col_types = cols(.default = col_character())) %>%
  filter(FILTER == 'PASS') %>% 
  mutate(ACMG_class = sub("full=", "", ACMG_class)) %>% 
  type_convert() %>%
  mutate(temp_SV_start = round(SV_start, -3), temp_SV_end = round(SV_end, -3)) %>%
  unite("variant", SV_chrom, temp_SV_start, temp_SV_end, SV_type, sep = "-", remove = FALSE) 
manta <- left_join(manta_original, manta_freq, by = c("variant") ) %>% 
  filter(is.na(CohortFreq) | CohortFreq < 0.025) %>% 
  filter( Annotation_mode == "split" | Gene_count == 0) %>%
  mutate(ACMG_class = case_when(!is.na(ACMG_class) & SV_type == "DEL" & !is.na(B_loss_source) ~ ACMG_class - 2,
                                !is.na(ACMG_class) & SV_type == "DUP" & !is.na(B_gain_source) ~ ACMG_class - 1,
                                SV_type == "BND" & ( ( SV_chrom %in% c("X", "chrX") & between(SV_start, 140420272, 140421278) | ( grepl("chrX|X", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 140420272, 140421278) ) )) ~ 5,
                                ( SV_chrom %in% c("17", "chr17") & between(SV_start, 59105674, 59683460) ) | ( grepl("chr17|17", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 59105674, 59683460) )  ~ 5,
                                TRUE ~ ACMG_class ) ) %>% 
  filter(ACMG_class > 1 | is.na(ACMG_class), is.na(SV_length) | abs(SV_length) < 1000000) %>% #added is.na(ACMG_class) 11/17/2021, may need to remove this part if too many lines in the results
  separate(Location, c('temp_location1', 'temp_location2'), sep = "-", remove = FALSE, convert = FALSE) %>% 
  filter(!(grepl("intron", temp_location1) & temp_location1 == temp_location2 & (!is.na(B_gain_source) | !is.na(B_loss_source) | !is.na(B_ins_source)) )) %>% 
  select(-starts_with('temp_'), -variant) %>% 
  mutate(ref_gene = toupper(Gene_name))
manta_sort <- left_join(manta, panelGene, by = c("ref_gene")) %>% 
  mutate(note = "") %>% 
  mutate(panel_class == ifelse(SV_type == "BND" & ( ( SV_chrom %in% c("X", "chrX") & between(SV_start, 140420272, 140421278) | ( grepl("chrX|X", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 140420272, 140421278) ) )), "Dx", panel_class)) %>% 
  mutate(eyeGene = case_when(panel_class == "Dx" ~ 2,
                             panel_class == "Candidate" ~ 1,
                             SV_type == "BND" & ( ( SV_chrom %in% c("X", "chrX") & between(SV_start, 140420272, 140421278) | ( grepl("chrX|X", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 140420272, 140421278) ) )) ~ 3,
                             ( SV_chrom %in% c("17", "chr17") & between(SV_start, 59105674, 59683460) ) | ( grepl("chr17|17", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 59105674, 59683460) )  ~ 2.5,
                             TRUE ~ 0)) %>% # GRCh38 coordinates
  arrange(desc(eyeGene), desc(ACMG_class)) %>% 
  mutate(note = ifelse(SV_type == "BND" & ( ( SV_chrom %in% c("X", "chrX") & between(SV_start, 140420272, 140421278) | ( grepl("chrX|X", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 140420272, 140421278) ) )), "XLFD", note)) %>% 
  mutate(note = ifelse(( SV_chrom %in% c("17", "chr17") & between(SV_start, 59105674, 59683460) ) | ( grepl("chr17|17", ALT) & between(as.numeric(gsub("\\D", "", ALT)), 59105674, 59683460) ), "RP17", note)) %>% 
  select(AnnotSV_ID:Gene_name,panel_class,ACMG_class,note,CohortFreq,`NaltP/NtotalP`,OMIM_phenotype:GnomAD_pLI,Frameshift,everything() )

openxlsx::write.xlsx(list("manta" = manta_sort), 
                     file = "D1440.manta.xlsx", firstRow = TRUE, firstCol = TRUE)
openxlsx::write.xlsx(list("manta" = manta_original), 
                     file = "D1440.manta.original.xlsx", firstRow = TRUE, firstCol = TRUE)
