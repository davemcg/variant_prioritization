#!/bin/bash

# file downloaded from https://github.com/Illumina/PrimateAI

cd /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation

gzcat PrimateAI_scores_v0.2.tsv.gz | head -n 100 | grep "primateDL_score" > primate.header
gzcat PrimateAI_scores_v0.2.tsv.gz | grep . | grep -v "#" | grep -v "primateDL_score" > primate.temp
python3 ~/git/ChromosomeMappings/convert_notation.py -c ~/git/ChromosomeMappings/GRCh37_UCSC2ensembl.txt -f primate.temp | sort -k1,1 -k2,2n > primate.temp2

cat primate.header primate.temp2 | bgzip -f > PrimateAI_scores_v0.2.tsv.reformatted.gz
tabix -f -S 1 -s 1 -b 2 -e 2 PrimateAI_scores_v0.2.tsv.reformatted.gz 
