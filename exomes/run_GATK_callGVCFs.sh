#!/bin/bash

module load GATK/3.5-0

input_bam=$1
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
