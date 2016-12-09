#!/bin/bash

##################
# Exome Workflow v02
# David McGaughey, NEI, NIH

# Transfers lane bam from NISC to subfolder per sample
# Takes the lane bam files, extracts fastq, aligns with bwa-mem, then creates merged bam
# Processes the merged bam with various Picard tools to sort, mark dups, and index
# Calls GVCF
# The resulting GVCF can be combined with other GVCFs with run_GATK_GVCFtoFilteredVCF.sh 

# IT IS HIGHLY RECOMMENDED THAT YOU ONLY GIVE ONE SAMPLE AT A TIME
# This workflow does the bwa mem alignment for the lane bams in serial

##################

# modules needed for python j1 job
module load samtools/1.3
module load bwa/0.7.12
module load picard/2.1.1

NISC_laneBam_matrix=$1 
sbatch_job_name=$2
exome_bait_bed=$3 #Give full path

############
# PREP
# Create scp job and sbatch job (j2)
~/bin/exome_workflow_v02/create_scp_and_sbatch_jobs_for_NISC_laneBams.py -f $1 -j $2 -b $3
chmod +x $2.scp.sh	 # make executable
./$2.scp.sh &		 # execute in background
wait				 # don't proceed until all transfers are complete
############

############
#RE-ALIGNMENT
# pulls bwa-formatted info (read group info, bam file, etc) then hands over to sbatch 
j1=$(sbatch --job-name bwa.$2 --time=12:00:00 --mem=30g --cpus-per-task 10 ~/bin/exome_workflow_v02/realign_NISC_laneBams_with_bwa.py $1)
############

############
# BAM processing and GVCF calling
# sort, mark duplicates, and index bam
sbatch --dependency=afterok:$j1 $2.sh
############

