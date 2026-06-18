#!/bin/bash
#SBATCH --gres=lscratch:50
#SBATCH --cpus-per-task=2
#SBATCH --mem=128g

#genome of 22 samples requies 67GB mem
head -n 1 gemini_tsv_filtered/"$(ls gemini_tsv_filtered | head -n 1)" > gemini_tsv_filtered/variant.enrich.input.txt

for file in gemini_tsv_filtered/*.filtered.tsv; do
	tail -n +2 $file >> gemini_tsv_filtered/variant.enrich.input.txt
	done

module load R/3.6.3
Rscript ~/git/variant_prioritization/dev/cohort_variant_enrichment.R gemini_tsv_filtered/variant.enrich.input.txt gemini_tsv_filtered/enriched.variant.xlsx

rm gemini_tsv_filtered/variant.enrich.input.txt
