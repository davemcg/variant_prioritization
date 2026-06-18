#!/bin/bash
#SBATCH --gres=lscratch:50
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g

#$1: tsv file, $2: vcf.gz file

module load samtools/1.17
awk -F"\t" 'BEGIN{OFS="\t"} {$2 = $2 +1; print "chr"$1,$2,$1"-"$2"-"$4"-"$5,$4,$5,"30","PASS",".","GT:GQ:DP","0/1:50:100"}' $1 | cat ~/git/variant_prioritization/dev/vcf.header.hg38 - | bgzip -c > $2

tabix -p vcf $2