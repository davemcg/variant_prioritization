#!/bin/bash
geneName=$1
mkdir -p geneSearch
head -n 1 gemini_tsv_filtered/$(ls gemini_tsv_filtered | head -n 1) > geneSearch/"$geneName".tsv 
for file in gemini_tsv_filtered/*.tsv; do grep -P "\t$geneName\t" $file >> geneSearch/"$geneName".tsv; done