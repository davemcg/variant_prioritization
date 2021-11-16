#!/bin/bash

~/git/ChromosomeMappings/convert_notation.py -c ~/git/ChromosomeMappings/GRCh37_UCSC2ensembl.txt -f /fdb/gemini/data/hg19.gerp.elements.bed.gz | bgzip -c > /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/hg19_ensembl.gerp.elements.bed.gz
tabix -p bed /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/hg19_ensembl.gerp.elements.bed.gz
