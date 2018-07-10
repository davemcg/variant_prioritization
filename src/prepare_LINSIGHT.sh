$!/bin/bash

# https://github.com/CshlSiepelLab/LINSIGHT
# http://biorxiv.org/content/early/2016/08/15/069682
wget http://compgen.cshl.edu/%7Eyihuang/tracks/LINSIGHT.bw . 
module load UCSC
bigWigToBedGraph LINSIGHT.bw LINSIGHT.bedGraph
# run with 32G of memory
~/git/ChromosomeMappings/convert_notation.py -c ~/git/ChromosomeMappings/GRCh37_UCSC2ensembl.txt -f LINSIGHT.bedGraph | awk -v OFS="\t" '{print $1, $2, $3, "LINSIGHT", $4}' | bgzip -f > LINSIGHT.GRCh37.bed.gz
tabix -p bed LINSIGHT.GRCh37.bed.gz
mv LINSIGHT.GRCh37.bed.gz* /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/
