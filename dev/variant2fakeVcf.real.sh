#!/bin/bash
#SBATCH --gres=lscratch:50
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g

#$1: tsv file, $2: vcf.gz file
#the first line was skipped
module load samtools/1.17
awk -F"\t" 'BEGIN{OFS="\t"} NR>1 {print $1,$2,$3,$4,$5,"30","PASS",".","GT:GQ:DP","0/1:50:100"}' $1 | cat ~/git/variant_prioritization/dev/vcf.header.hg38 - | bgzip -c > $2

tabix -p vcf -f $2
