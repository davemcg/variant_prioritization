#!/bin/bash
#SBATCH --gres=lscratch:50
#SBATCH --cpus-per-task=8
#SBATCH --mem=8g
#SBATCH --partition=quick
#SBATCH --time=2:0:0

# Input parameters
INPUT_VCF=$1
SAMPLE_LIST=$2
OUTPUT_VCF=$3

module load bcftools

# Step 1: Extract variants for specific samples
# Read the sample list into a comma-separated string
SAMPLES=$(paste -sd, $SAMPLE_LIST)

# Combine all steps into a single bcftools command
bcftools view -s $SAMPLES $INPUT_VCF \
    | bcftools +fill-tags -Ou -- -t AC,AC_Hom,AC_Het,AN,AF \
    | bcftools view --threads 4 -e 'AC=0' -Ou \
    | bcftools view --threads 4 -e 'F_MISSING=1' -Ou \
    | bcftools view --threads 4 -i 'ID~"fbDv"' -Oz -o $OUTPUT_VCF

# Index the output VCF file
tabix -f -p vcf $OUTPUT_VCF

echo "Variants extraction and filtering completed. Output saved to $OUTPUT_VCF"
