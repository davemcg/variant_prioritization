#!/bin/bash
#SBATCH --gres=lscratch:50
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g

geneList=$1 #i.e.: MAC
geneList_file=$2

mkdir -p geneSearch
head -n 1 $(ls gemini_tsv_filtered/*.filtered.tsv | head -n 1) > geneSearch/"$geneList".tsv 


while read -r geneName;
do
echo "$geneName";
for file in gemini_tsv_filtered/*.filtered.tsv; do grep -e $'^'$geneName$'\t' -e $'\t'$geneName$'\t' $file >> geneSearch/"$geneList".tsv; done;
done < $geneList_file

module load R/3.6.3
Rscript ~/git/variant_prioritization/dev/cohort_variant_enrichment.R geneSearch/"$geneList".tsv geneSearch/"$geneList".enriched.variant.xlsx

