#!/bin/bash
geneName=$1
mkdir -p geneSearch
head -n 1 $(ls exome_tsv/*.tsv | head -n 1) > geneSearch/"$geneName".tsv 
for file in exome_tsv/*.tsv; do grep $'\t'$geneName$'\t' $file >> geneSearch/"$geneName".tsv; done
for file in genome_tsv/*.tsv; do grep $'\t'$geneName$'\t' $file >> geneSearch/"$geneName".tsv; done

#while read -r gene; do bash ~/git/variant_prioritization/geneSearch.sh $gene; done < genelist.tsv
#geneList.tsv file may require Unix file EOL symbols.
#for gene in X,Y,Z; do bash ~/git/variant_prioritization/geneSearch.sh $gene; done