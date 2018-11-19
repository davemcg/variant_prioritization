#!/bin/bash

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
zcat refFlat.txt.gz | awk -v OFS='\t' '{print $3, $5, $6, $1}' | sort -k1,1 -k2,2n | uniq > temp_file
~/git/ChromosomeMappings/convert_notation.py -c ~/git/ChromosomeMappings/GRCh37_UCSC2ensembl.txt -f temp_file > ensembl_genes.bed
rm temp_file
rm refFlat.txt.gz
