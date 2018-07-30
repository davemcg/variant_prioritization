# move GTEx / eyeIntegration data into matrix for eye_var_Pathogenicity
# files on mcgaughey 'eyeMac' tower
# use output in build_eyeIntegration_TPM02.py

library(dplyr)
library(tidyr)
library(RSQLite)
#eye_and_gtex_lsTPM$GeneName <- row.names(eye_and_gtex_lsTPM)
# shiny_data <- eye_and_gtex_lsTPM
sqlite_file <- "~/git/Human_eyeIntegration_App/www/eyeIntegration_human.sqlite"
sqldb <- dbConnect(SQLite(), dbname=sqlite_file)
load('~/git/Human_eyeIntegration_App/www/core_tight.Rdata')
data = dbGetQuery(sqldb, 'select * from lsTPM') %>% left_join(core_tight)

sum_data <- data %>% group_by(GeneName, Sub_Tissue) %>% summarise(mean_lsTPM = mean(value)) %>% rowwise %>% mutate(Sub_Tissue = trimws(Sub_Tissue))

table_eyeIntegration_TPM <- sum_data %>% spread(Sub_Tissue, mean_lsTPM)
colnames(table_eyeIntegration_TPM) <- gsub('-','_', colnames(table_eyeIntegration_TPM))
colnames(table_eyeIntegration_TPM) <- gsub(' ','', colnames(table_eyeIntegration_TPM))
colnames(table_eyeIntegration_TPM) <- gsub('\\(','_', colnames(table_eyeIntegration_TPM))
colnames(table_eyeIntegration_TPM) <- gsub('\\)','', colnames(table_eyeIntegration_TPM))

# remove NA col
table_eyeIntegration_TPM <- table_eyeIntegration_TPM[,-56]
write.table(table_eyeIntegration_TPM, '~/git/variant_prioritization/data/eyeIntegration_TPM.tsv', quote = F, sep ="\t", row.names = F)
system('bgzip ~/git/variant_prioritization/data/eyeIntegration_TPM.tsv')
