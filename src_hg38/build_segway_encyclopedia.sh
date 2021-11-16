#!/bin/bash

# https://noble.gs.washington.edu/proj/encyclopedia/
# http://dx.doi.org/10.1101/086025

cd /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation

wget https://noble.gs.washington.edu/proj/encyclopedia/segway_encyclopedia.bed.gz .

zgrep -v '^chrom' segway_encyclopedia.bed.gz | \
	~/git/ChromosomeMappings/convert_notation.py -c ~/git/ChromosomeMappings/GRCh37_UCSC2ensembl.txt -f - | \
	sort -k1,1 -k2,2n > segway_encyclopedia.TEMP

zgrep '^chrom' segway_encyclopedia.bed.gz > header.temp
cat header.temp segway_encyclopedia.TEMP | bgzip > segway_encyclopedia.ready.bed.gz

tabix -S 1 -f -s 1 -b 2 -e 3 -0 segway_encyclopedia.ready.bed.gz

	
