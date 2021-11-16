#!/bin/bash

cd /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/
wget https://epilogos.altiusinstitute.org/assets/epilogos/v06_16_2017/hg19/15/group/all.KL.bed.gz .

python3 ~/git/variant_prioritization/src/build_epilogos.py \
		/data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/all.KL.bed.gz | \
	~/git/ChromosomeMappings/convert_notation.py \
		-c ~/git/ChromosomeMappings/GRCh37_UCSC2ensembl.txt -f - | \
    bgzip > all.KL.reformatted.bed.gz

tabix -f -s 1 -b 2 -e 3 -0 all.KL.reformatted.bed.gz

