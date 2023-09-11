#!/bin/bash
geneName=$1
mkdir -p geneSearch
head -n 1 $(ls gemini_tsv_filtered/*.filtered.tsv | head -n 1) > geneSearch/"$geneName".tsv 
for file in gemini_tsv_filtered/*.tsv; do grep -e $'^'$geneName$'\t' -e $'\t'$geneName$'\t' $file >> geneSearch/"$geneName".tsv; done

#while read -r gene; do bash ~/git/variant_prioritization/geneSearch.sh $gene; done < genelist.tsv
#geneList.tsv file may require Unix file EOL symbols.
#for gene in X,Y,Z; do bash ~/git/variant_prioritization/geneSearch.sh $gene; done