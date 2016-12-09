#!/bin/bash

module load picard

# Sort bam, then marks duplicates, then creates index

echo $1

# "Soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads"
java -Xmx20g -jar $PICARDJARPATH/picard.jar CleanSam \
	INPUT=$1 OUTPUT=${1%.bam}.CleanSam.bam


# "Verify mate-pair information between mates and fix if needed."
# also coord sorts
java -Xmx20g -jar $PICARDJARPATH/picard.jar FixMateInformation \
	INPUT=${1%.bam}.CleanSam.bam OUTPUT=${1%.bam}.sorted.bam SORT_ORDER=coordinate

# name for easier downstream use
sorted_bam=${1%.bam}.sorted.bam

# Mark dups
java -Xmx20g -jar $PICARDJARPATH/picard.jar MarkDuplicates \
	INPUT=$sorted_bam OUTPUT=${sorted_bam%.bam}.markDup.bam METRICS_FILE=${sorted_bam%.bam}.markDup.metrics

sorted_markDup_bam=${sorted_bam%.bam}.markDup.bam

# Build bam index
java -Xmx20g -jar $PICARDJARPATH/picard.jar BuildBamIndex \
	INPUT=$sorted_markDup_bam OUTPUT=$sorted_markDup_bam.bai

# Delete intermediate files
rm $1
rm ${1%.bam}.CleanSam.bam
rm ${1%.bam}.sorted.bam
#rm ${sorted_bam%.bam}.markDup.bam
