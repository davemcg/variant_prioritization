#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import subprocess
import collections
from collections import defaultdict
import sys
import glob
import os

parser = argparse.ArgumentParser()
parser.add_argument('filelist', help= 
	"Takes in list of NISC lane bams, transfer the bams over \
	from Trek, extracts fastq, realigns fastq with bwa-mem \
	to 1000g phaseII hs37d5 with read group info extracted \
	from the lane bam. \
	\
	I HIGHLY SUGGEST YOU ONLY DO ONE SAMPLE AT A TIME!!!! \
	This script does NOT parallelize any of the tasks. \
	\
	Invoke script with sbatch --mem=30G --cpus-per-task 10 \
	\
	Example (of the input list): \
		CCGO_800014	150223_OPTIMUS_C6HGHANXX.8.9477663 \
		CCGO_800014	150306_YOSHI_C6HDMANXX.5.9477663 \
		CCGO_800014	150223_OPTIMUS_C67T9ANXX.7.9477663 \
		CCGO_800014	150223_OPTIMUS_C67T9ANXX.8.9477663 \
		CCGO_800014	150223_OPTIMUS_C67T9ANXX.6.9477663 \
		CCGO_800015	150223_OPTIMUS_C67T9ANXX.7.9477645 \
		CCGO_800015	150223_OPTIMUS_C67T9ANXX.8.9477645 \
		CCGO_800015	150223_OPTIMUS_C67T9ANXX.6.9477645 \
		CCGO_800015	150306_YOSHI_C6HDMANXX.5.9477645 \
		CCGO_800015	150223_OPTIMUS_C6HGHANXX.8.9477645 \
	\
	Example (exome): \
   		sbatch --mem=20G --cpus-per-task 10 this_script.py list_of_laneBams.txt")

args = parser.parse_args()
bamlist = args.filelist

# loop through the bamlist, making a dict
sample_laneBam = collections.defaultdict(list)
for line in open(bamlist):
	line = line.split()
	if len(line) != 2:
		print("Input list formatting issue\n\n", line)
		sys.exit(0)
	sample_laneBam[line[0]].append(line[1])

# As I can't figure out a way to 
# create subfolders for each sample
# Use scp to move lane Bams from Trek to biowulf2
samples = []
[samples.append(k) for k,v in sample_laneBam.items()]
samples.sort()
for one_sample in samples:
	if not os.path.isdir(one_sample):
		print("Subfolder for ", one_sample, " missing!!!!")
		sys.exit(0)
	# create Trek directory structure
	base_dir = "/cluster/ifs/projects/solexa/reads/"
	# loop through each laneBam and do the scp
	for laneBam in sample_laneBam[one_sample]:
		# extract read group info from bam
		samtools_header_call = 'samtools view -H ' + one_sample + '/' + laneBam + '.bam ' +  '| grep ^@RG'
		read_group = (subprocess.check_output(samtools_header_call, shell=True)).decode('utf-8')
		# example read_group:
		# '@RG\tID:3\tSM:AMS0013\tLB:AMS0013_1\tPU:160129_OPTIMUS_C7NP2ANXX.3.11457313\tCN:NISC\tDT:2016-02-04T20:49:40-05:00\tPL:HiSeq2000\n'
		# replace SM with sample name
		read_group = read_group.split('\t')
		read_group[2] = 'SM:' + one_sample
		read_group[-1] = 'PL:ILLUMINA'
		# https://github.com/IARCbioinfo/BAM-tricks
		# https://github.com/samtools/samtools/issues/532#issuecomment-205877064
		# now shuffle, extract interleaved fastq, and run bwa-mem
		big_run = 'samtools collate -uOn 128 ' + one_sample + '/' + laneBam + '.bam' + ' /tmp/TMP-' + laneBam + \
				  '| samtools fastq - | bwa mem -t 10 -B 4 -O 6 -E 1 -M -p -R ' + \
				  "\"" + '\\t'.join(read_group) + "\"" + \
			 	  ' /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta -\
				  | samtools view -1 - > ' + one_sample + '/' + laneBam + '.bwa-mem.b37.bam' 
		bwa_commands = open(one_sample + '/' + one_sample + '.realignmentCommands.txt','a')
		bwa_commands.write(big_run + '\n')
		subprocess.check_call(big_run,shell=True)
		# rm NISC lane bams
		rm_call = 'rm ' + one_sample + '/' + laneBam + '.bam*'
		subprocess.check_call(rm_call,shell=True)
# merge laneBams into one bam for downstream GATK workflow
for one_sample in samples:
	bams = glob.glob(one_sample + '/*b37.bam')
	gatk_bam_input = []
	for i in bams:
		gatk_bam_input.append("I=" + i)
	MergeSamFiles_call = 'java -Xmx20g -jar $PICARDJARPATH/picard.jar MergeSamFiles ' + \
						 ' '.join(gatk_bam_input) + ' O=' + one_sample + '/' + one_sample + \
						 '.bwa-mem.b37.merged.bam'
	subprocess.check_call(MergeSamFiles_call, shell=True)
	bwa_commands = open(one_sample + '/' + one_sample + '.realignmentCommands.txt','a')
	bwa_commands.write('\n\n' + MergeSamFiles_call)
