library(tidyverse)
library(readxl)

args <- commandArgs(trailingOnly=TRUE)
#setwd("W:/abca4/clinvar.hgmd")
args <- c("Z:/resources/OGLpanelGeneDxORcandidate.xlsx", 
          "Z:/resources/omim/omim_API_eye_phenotype/omim.20250916.csv",
          "Z:/resources/omim/omim_API_eye_phenotype/omim.20250916.DxORcandidate.xlsx")

geneCategory_file <- args[1]
omim_file <- args[2]
output_excel <- args[3]


panelGene <- read_xlsx(geneCategory_file, sheet = "analysis", na = c("NA", "", "None", "NONE", ".")) %>% 
  select(-gene_mim_number, -gene_ids,-eye_clinical_terms,-eye_phenotype)
#  select(gene, panel_class, GenePhenotypeCategory)

##check the phenotype category:
#merged <- panelGene %>% summarise(category = str_c(unique(na.omit(str_trim(GenePhenotypeCategory))), collapse=",")) %>% pull(category) %>% str_split(",") %>%  unlist() %>% str_trim() %>% discard(~ .x == "") %>% {paste(unique(.), collapse = ",")}

#2576 omim genes with eye phenotype on 9/16/2025
inheritance_order <- c("AD","AR","XL","XLD","XR","MT")

omim <- read_csv(omim_file, col_names = TRUE, na = c("NA", "", "None", "NONE", ".", "FALSE", "False"), col_types = cols(.default = col_character())) %>%
  mutate(temp_eye_clinical_terms = eye_clinical_terms) %>% 
  mutate(eye_clinical_terms = gsub("\\|", "|\r", eye_clinical_terms)) %>% 
  separate_rows(temp_eye_clinical_terms, sep = "\\|") %>% 
  separate(temp_eye_clinical_terms, c("Phenotype", "Phenotype_MIM_number", "Inheritance_omim", "Phenotype_mapping_key", "eye_phenotype"),
  sep = "::", remove = TRUE) %>%
  filter(!(gene_symbols == "ABCA4" & Inheritance_omim == "AD")) %>%
  select(-Phenotype, -Phenotype_MIM_number, -Phenotype_mapping_key) %>% 
  group_by(gene_mim_number, gene_ids, gene_symbols, eye_clinical_terms) %>% 
  summarise(eye_phenotype = str_c(unique(na.omit(str_trim(eye_phenotype))), collapse = " | "),
            Inheritance_omim = str_c(unique(na.omit(str_trim(Inheritance_omim))), collapse = "|"),
             .groups = "drop_last") %>% 
  mutate(
    Inheritance_omim = Inheritance_omim %>%
      str_replace_all(";", "|") %>%
      str_split("\\|") %>%
      map(~ .x %>%
            discard(is.na) %>%
            str_trim() %>%
            discard(~ .x == "") %>%
            unique() %>%
            { .[order(match(., inheritance_order, nomatch = length(inheritance_order) + 1),
                      toupper(.))] }   # preferred first, others alpha
      ) %>%
      map_chr(~ str_c(.x, collapse = "|"))
  ) %>% 
  rename(gene = gene_symbols) #%>% 
 # filter(!eye_phenotype %in% c("No ocular symptoms", "Normal eyes"
 #                              ))

##multiple MAC in the GenePhenotypeCategory for several genes, such GNB3, GRK1, why? 
newGeneTable <- full_join(panelGene, omim, by = "gene") %>%
  unite("Inheritance", c("Inheritance", "Inheritance_omim"), sep = "|", na.rm = TRUE, remove = TRUE) %>%
  #mutate(OGL_Phenotypes = ifelse(grepl("::", Phenotypes), NA, Phenotypes)) %>% #previous version included disease names by manual editing, #only use during the 1st trial to make v1 
  unite("temp_phenotype", c("OGL_Phenotypes", "eye_clinical_terms"), sep = "|", na.rm = TRUE, remove = FALSE) %>% #temp_phenotype includes the disease name and use for gene grouping below.
  unite("Phenotypes", c("OGL_Phenotypes", "eye_phenotype"), sep = "::", na.rm = TRUE, remove = FALSE) %>% #OGL annotation and omim separated by "::"
  mutate(
    Inheritance = Inheritance %>%
      str_replace_all(",", "|") %>%
      str_split("\\|") %>%
      map(~ .x %>%
            discard(is.na) %>%
            str_trim() %>%
            discard(~ .x == "") %>%
            unique() %>%
            { .[order(match(., inheritance_order, nomatch = length(inheritance_order) + 1),
                      toupper(.))] }   # preferred first, others alpha
      ) %>%
      map_chr(~ str_c(.x, collapse = "|"))
  ) %>% 
  mutate(GenePhenotypeCategory = sub("Retinopathy", "RD", GenePhenotypeCategory, ignore.case = TRUE),
         GenePhenotypeCategory = sub("Maculopathy", "MD", GenePhenotypeCategory, ignore.case = TRUE)) %>% 
  mutate(panel_class = case_when( !is.na(panel_class) ~ panel_class,
                                  is.na(panel_class) & grepl("Joubert|Usher", eye_clinical_terms) ~ "Dx",
                                  is.na(panel_class) & grepl("retinopathy|dystrophy|degeneration|pigmentosa|microphthalmia|anophthalmia|coloboma|Ectopia lentis|aniridia|night|cataract|myopia|atrophy|hypoplasia|glaucoma|nystagmus|strabismus", eye_clinical_terms, ignore.case = TRUE) ~ "Candidate-High",
                                  TRUE ~ NA ))  %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("coloboma|anophthalmia|microphthalmia", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)MAC(,|$)", GenePhenotypeCategory), 
                                           ifelse(is.na(GenePhenotypeCategory), "MAC",
                                                  paste0(GenePhenotypeCategory, ",MAC")),
                                           GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("retinopathy|retinal dystrophy|rod dystrophy|cone dystrophy|degeneration|pigmentosa|Night blindness|Vision loss, progressive", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)RD(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "RD",
                                               paste0(GenePhenotypeCategory, ",RD")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("macular dystrophy|maculopathy", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)MD(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "MD",
                                               paste0(GenePhenotypeCategory, ",MD")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Macular degeneration, age-related|Age-related macular degeneration", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)AMD(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "AMD",
                                               paste0(GenePhenotypeCategory, ",AMD")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Foveal hypoplasia", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Foveal hypoplasia(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Foveal hypoplasia",
                                               paste0(GenePhenotypeCategory, ",Foveal hypoplasia")),
                                        GenePhenotypeCategory )) %>%
  mutate(GenePhenotypeCategory = ifelse(grepl("Glaucoma", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Glaucoma(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Glaucoma",
                                               paste0(GenePhenotypeCategory, ",Glaucoma")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Usher", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Usher(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Usher",
                                               paste0(GenePhenotypeCategory, ",Usher")),
                                        GenePhenotypeCategory )) %>%
  mutate(GenePhenotypeCategory = ifelse(grepl("Joubert", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Joubert(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Joubert",
                                               paste0(GenePhenotypeCategory, ",Joubert")),
                                        GenePhenotypeCategory )) %>%
  mutate(GenePhenotypeCategory = ifelse(grepl("Bardet-Biedl syndrome", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)BBS(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "BBS",
                                               paste0(GenePhenotypeCategory, ",BBS")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Hermansky-Pudlak syndrome", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)HPS(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "HPS",
                                               paste0(GenePhenotypeCategory, ",HPS")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("CHARGE", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)CHARGE(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "CHARGE",
                                               paste0(GenePhenotypeCategory, ",CHARGE")),
                                        GenePhenotypeCategory )) %>%
  mutate(GenePhenotypeCategory = ifelse(grepl("Albinism", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Albinism(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Albinism",
                                               paste0(GenePhenotypeCategory, ",Albinism")),
                                        GenePhenotypeCategory )) %>%
  mutate(GenePhenotypeCategory = ifelse(grepl("Optic Atrophy|Optic neuropathy|Optic nerve atrophy", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)OA/ON(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "OA/ON",
                                               paste0(GenePhenotypeCategory, ",OA/ON")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Achromatopsia", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Achromatopsia(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Achromatopsia",
                                               paste0(GenePhenotypeCategory, ",Achromatopsia")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Aniridia", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Aniridia(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Aniridia",
                                               paste0(GenePhenotypeCategory, ",Aniridia")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("FEVR", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)FEVR(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "FEVR",
                                               paste0(GenePhenotypeCategory, ",FEVR")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Nanophthalmos", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Nanophthalmos(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Nanophthalmos",
                                               paste0(GenePhenotypeCategory, ",Nanophthalmos")),
                                        GenePhenotypeCategory )) %>%
  mutate(GenePhenotypeCategory = ifelse(grepl("High hyperopia", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)High hyperopia(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "High hyperopia",
                                               paste0(GenePhenotypeCategory, ",High hyperopia")),
                                        GenePhenotypeCategory )) %>%
  mutate(GenePhenotypeCategory = ifelse(grepl("Cataract", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Cataract(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Cataract",
                                               paste0(GenePhenotypeCategory, ",Cataract")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Nystagmus", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Nystagmus(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Nystagmus",
                                               paste0(GenePhenotypeCategory, ",Nystagmus")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Strabismus", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Strabismus(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Strabismus",
                                               paste0(GenePhenotypeCategory, ",Strabismus")),
                                        GenePhenotypeCategory )) %>%
  mutate(GenePhenotypeCategory = ifelse(grepl("Ectopia lentis", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Ectopia lentis(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Ectopia lentis",
                                               paste0(GenePhenotypeCategory, ",Ectopia lentis")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Stickler", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Stickler(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Stickler",
                                               paste0(GenePhenotypeCategory, ",Stickler")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("HPS", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)HPS(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "HPS",
                                               paste0(GenePhenotypeCategory, ",HPS")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Anterior segment dysgenesis", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)ASD(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "ASD",
                                               paste0(GenePhenotypeCategory, ",ASD")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Myopia", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Myopia(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Myopia",
                                               paste0(GenePhenotypeCategory, ",Myopia")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Corneal dystrophy", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Corneal dystrophy(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Corneal dystrophy",
                                               paste0(GenePhenotypeCategory, ",Corneal dystrophy")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Ciliopathy", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Ciliopathy(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Ciliopathy",
                                               paste0(GenePhenotypeCategory, ",Ciliopathy")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("congenital stationary|CSNB", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)CSNB(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "CSNB",
                                               paste0(GenePhenotypeCategory, ",MAC")),
                                        GenePhenotypeCategory )) %>% 
  mutate(GenePhenotypeCategory = ifelse(grepl("Mitochondrial", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Mito(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Mito",
                                               paste0(GenePhenotypeCategory, ",Mito")),
                                        GenePhenotypeCategory )) %>%
  mutate(GenePhenotypeCategory = ifelse(grepl("Cortical visual impairment", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)CVI(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "CVI",
                                               paste0(GenePhenotypeCategory, ",CVI")),
                                        GenePhenotypeCategory )) %>%
  mutate(GenePhenotypeCategory = ifelse(grepl("Ophthalmoplegia", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Ophthalmoplegia(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Ophthalmoplegia",
                                               paste0(GenePhenotypeCategory, ",Ophthalmoplegia")),
                                        GenePhenotypeCategory )) %>%
  mutate(GenePhenotypeCategory = ifelse(grepl("Syndrome", temp_phenotype, ignore.case = TRUE) & !grepl("(^|,\\s*)Syndrome(,|$)", GenePhenotypeCategory), 
                                        ifelse(is.na(GenePhenotypeCategory), "Syndrome",
                                               paste0(GenePhenotypeCategory, ",Syndrome")),
                                        GenePhenotypeCategory )) %>%
  mutate(GenePhenotypeCategory = ifelse(gene == "ADAM9", sub("MAC,", "", GenePhenotypeCategory), GenePhenotypeCategory)) %>% 
  select(-starts_with("temp_")) %>% 
  select(Gene_Queue_by_time:Inheritance, OGL_Phenotypes,Phenotypes, everything())

#intermediate step to add OGL_Phenotype column to possibly manually edited file
# Intermediate_panelGene <- read_xlsx(geneCategory_file, sheet = "analysis", na = c("NA", "", "None", "NONE", ".")) %>% 
#   select(-gene_mim_number, -gene_ids,-eye_clinical_terms,-eye_phenotype)
# 
# Intermediate_newGeneTable <- left_join(Intermediate_panelGene, select(newGeneTable, gene, OGL_Phenotypes))
# openxlsx::write.xlsx(list("updated" = Intermediate_newGeneTable), file = "test2.xlsx", firstRow = TRUE, firstCol = TRUE)

openxlsx::write.xlsx(list("updated" = newGeneTable), file = "test.xlsx", firstRow = TRUE, firstCol = TRUE)

