#!/bin/bash

module load picard/2.1.1
module load GATK/3.5-0

###############################################################
# Picard
# cleaning, verify mate-pair, marking dups, and creating index
#
# Input Bam  must be aligned against GRCh37 (e.g. chromsomes
# are 1, not chr1. 
###############################################################

# "Soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads"
java -Xmx20g -jar $PICARDJARPATH/picard.jar CleanSam \
    INPUT=$1 OUTPUT=${1%.bam}.CleanSam.bam
rm $1

# "Verify mate-pair information between mates and fix if needed."
# also coord sorts
java -Xmx20g -jar $PICARDJARPATH/picard.jar FixMateInformation \
    INPUT=${1%.bam}.CleanSam.bam OUTPUT=${1%.bam}.sorted.bam SORT_ORDER=coordinate
rm ${1%.bam}.CleanSam.bam

# name for easier downstream use
sorted_bam=${1%.bam}.sorted.bam

# Mark dups
java -Xmx20g -jar $PICARDJARPATH/picard.jar MarkDuplicates \
    INPUT=$sorted_bam OUTPUT=${sorted_bam%.bam}.markDup.bam METRICS_FILE=${sorted_bam%.bam}.markDup.metrics
sorted_markDup_bam=${sorted_bam%.bam}.markDup.bam

# Build bam index
java -Xmx20g -jar $PICARDJARPATH/picard.jar BuildBamIndex \
    INPUT=$sorted_markDup_bam OUTPUT=$sorted_markDup_bam.bai
rm ${1%.bam}.sorted.bam


###############################################################
# GATK
# realign, recalibrate, call GVCF
###############################################################
input_bam=$sorted_markDup_bam
exome_bait_bed=$2

# Takes ~ 90 minutes
GATK -m 8g RealignerTargetCreator \
    -R /fdb/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta  \
    -I $input_bam \
    --known /fdb/GATK_resource_bundle/b37-2.8/1000G_phase1.indels.b37.vcf \
    --known /fdb/GATK_resource_bundle/b37-2.8/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -o ${input_bam%.bam}.forIndexRealigner.intervals \
    -L $2 \
    --interval_padding 100

# Takes ~ 100 minutes
GATK -m 8g IndelRealigner \
    -R /fdb/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
    -I $input_bam \
    --knownAlleles /fdb/GATK_resource_bundle/b37-2.8/1000G_phase1.indels.b37.vcf \
    --knownAlleles /fdb/GATK_resource_bundle/b37-2.8/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -targetIntervals ${input_bam%.bam}.forIndexRealigner.intervals \
    -o ${input_bam%.bam}.realigned.bam
rm $input_bam

GATK -m 8g BaseRecalibrator \
    -R /fdb/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
    -I ${input_bam%.bam}.realigned.bam \
    --knownSites /fdb/GATK_resource_bundle/b37-2.8/dbsnp_138.b37.excluding_sites_after_129.vcf \
    --knownSites /fdb/GATK_resource_bundle/b37-2.8/1000G_phase1.indels.b37.vcf \
    --knownSites /fdb/GATK_resource_bundle/b37-2.8/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -L $2 \
    --interval_padding 100 \
    -o ${input_bam%.bam}.recal_data.table1

GATK -m 8g PrintReads \
    -R /fdb/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
    -I ${input_bam%.bam}.realigned.bam \
    -BQSR ${input_bam%.bam}.recal_data.table1 \
    -o ${input_bam%.bam}.realigned.recalibrated.bam

GATK -m 8g BaseRecalibrator \
    -R /fdb/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
    -I ${input_bam%.bam}.realigned.bam \
    --knownSites /fdb/GATK_resource_bundle/b37-2.8/dbsnp_138.b37.excluding_sites_after_129.vcf \
    --knownSites /fdb/GATK_resource_bundle/b37-2.8/1000G_phase1.indels.b37.vcf \
    --knownSites /fdb/GATK_resource_bundle/b37-2.8/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -BQSR ${input_bam%.bam}.recal_data.table1 \
    -o ${input_bam%.bam}.recal_data.table2
rm ${input_bam%.bam}.realigned.bam

GATK -m 8g AnalyzeCovariates \
    -R /fdb/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
    -before ${input_bam%.bam}.recal_data.table1 \
    -after  ${input_bam%.bam}.recal_data.table2 \
    -plots ${input_bam%.bam}.BQSRplots.pdf

# Takes ~ 180 minutes
GATK -m 8g HaplotypeCaller \
    -R /fdb/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
    -I ${input_bam%.bam}.realigned.recalibrated.bam \
    -L $2 \
    --interval_padding 100 \
    --emitRefConfidence GVCF \
    -BQSR ${input_bam%.bam}.recal_data.table1 \
    -o ${input_bam%.bam}.realigned.raw.g.vcf.gz
