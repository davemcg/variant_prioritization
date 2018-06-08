#!/bin/bash

revel_file="revel_all_chromosomes.csv"
cd /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/
if [ ! -f $revel_file ]; then
    wget https://rothsj06.u.hpc.mssm.edu/revel/revel_all_chromosomes.csv.zip
	unzip revel_all_chromosomes.csv.zip
fi


echo "##INFO=<ID=REVEL_SCORE, Number=1,Type=Float>" > revel_all_chromosomes.vcf
echo -e "#CHROM\tPOS\tREF\tALT\tINFO" >> revel_all_chromosomes.vcf
cat revel_all_chromosomes.csv | grep -v '^chr,' | awk -F"," -v OFS="" '{print $1,"\t", $2, "\t", $3, "\t",  $4, "\tREVEL_SCORE=", $7}' >> revel_all_chromosomes.vcf
bgzip revel_all_chromosomes.vcf
tabix -p vcf revel_all_chromosomes.vcf.gz

