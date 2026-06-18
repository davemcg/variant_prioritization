#!/bin/bash
#SBATCH --gres=lscratch:50
#SBATCH --cpus-per-task=8
#SBATCH --mem=8g
#SBATCH --partition=quick
#SBATCH --time=2:0:0

# Input parameters
INPUT_VCF=$1
SAMPLE_LIST=$2 #"comma separated list"
REGION=$3 #chr14:80761505-80914366
OUTPUT_tsv=$4
module load samtools
#--min-ac 1 was used to remove the variants absent in either.

bcftools view -s $SAMPLE_LIST -r $REGION --min-ac 1 -i 'ID~"fbDv"' -Ov $INPUT_VCF | bcftools view --genotype het | grep -v "^##" | awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print $0,"dist_POS"; next} {d = (NR==2 || $1!=pchr) ? 0 : ($2-ppos); print $0,d; pchr=$1; ppos=$2}' > $OUTPUT_tsv

#the max distance
echo "The two largest distances between het variants are"
tail -n +2 $OUTPUT_tsv | cut -f 11 | sort -n | tail -n 2