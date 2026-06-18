#!/bin/bash
#SBATCH --gres=lscratch:50
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g

#$1: tsv file, $2: vcf.gz file
#the first line was skipped
mkdir -p temp
grep -v "^#" $1 > temp/validator_output.tsv
module load R/4.3.0
Rscript ~/git/variant_prioritization/dev/variant_validator2vcfS1.v1.R temp/validator_output.tsv temp/validator_vcf.tsv
module load samtools/1.17
awk -F"\t" 'BEGIN{OFS="\t"} NR>1 {print $1,$2,$3,$4,$5,"30","PASS",".","GT:GQ:DP","0/1:50:100"}' temp/validator_vcf.tsv | cat ~/git/variant_prioritization/dev/vcf.header.hg38 - | bgzip -c > $2

tabix -p vcf -f $2
rm -rf temp