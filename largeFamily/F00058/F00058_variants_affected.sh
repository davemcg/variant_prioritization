#!/bin/bash
#SBATCH -c8
#SBATCH --mem=64g
#SBATCH --gres=lscratch:50
#SBATCH --time=2:0:0

module load R/4.3.0
mkdir -p family
Rscript family_multiple_affected_filter_F00058.R gemini_tsv_filtered/G00582.prasov23-1.23sample.gemini.filtered.tsv gemini_tsv_filtered/G00594.prasov23-1.23sample.gemini.filtered.tsv gemini_tsv_filtered/G00595.prasov23-1.23sample.gemini.filtered.tsv gemini_tsv_filtered/G00614.prasov23-1.23sample.gemini.filtered.tsv  gemini_tsv_filtered/G00628.prasov23-1.23sample.gemini.filtered.tsv gemini_tsv_filtered/G00572.prasov23-1.23sample.gemini.filtered.tsv gemini_tsv_filtered/G00598.prasov23-1.23sample.gemini.filtered.tsv  gemini_tsv_filtered/G00602.prasov23-1.23sample.gemini.filtered.tsv gemini_tsv_filtered/G00609.prasov23-1.23sample.gemini.filtered.tsv gemini_tsv_filtered/G00626.prasov23-1.23sample.gemini.filtered.tsv
 