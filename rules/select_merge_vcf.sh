#!/bin/bash
#SBATCH --gres=lscratch:100
#SBATCH --cpus-per-task=8
#SBATCH --mem=8g
#SBATCH --partition=norm,quick
#SBATCH --time=2:0:0

#took 20 min for three genome samples
key=$1 #sample_vcf.key.tsv, first column: sampleID; 2nd column: vcf file name. Separate multiple samples with ",". The vcf should be in /data/OGL/resources/OGLsample/annotatedVCF
vcf_output=$2 #output.vcf.gz

set -e

module load samtools/1.21 parallel

LWORK_DIR=/lscratch/$SLURM_JOB_ID
mkdir -p $LWORK_DIR/vcf

grep -v "^#" $key | sed 's/\r$//' | parallel -C "\t" -j $SLURM_CPUS_PER_TASK "bcftools view --samples {1} /data/OGL/resources/OGLsample/annotatedVCF/{2} --output-type z -o $LWORK_DIR/vcf/{1}.vcf.gz && tabix -p vcf $LWORK_DIR/vcf/{1}.vcf.gz"

#sed 's/\r$//' |

bcftools merge --threads 4 --merge none --missing-to-ref --output-type u $LWORK_DIR/vcf/*.vcf.gz \
 | bcftools +fill-tags - -Ou -- -t AC,AC_Hom,AC_Het,AN,AF \
 | bcftools filter -i 'INFO/AC > 0' -Ov \
 | sed 's#0/0:\.:\.:\.#0/0:10:10:10,0#g' - \
 | bgzip -@ 4 > $vcf_output

tabix -@ 4 -p vcf $vcf_output
