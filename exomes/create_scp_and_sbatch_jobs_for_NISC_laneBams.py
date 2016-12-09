#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import subprocess
import collections
from collections import defaultdict
import sys
import glob
import os

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filelist', help=
    "Takes in list of NISC lane bams, then creates a \
	bash job transfer the bams over \
    from Trek via scp into subfolders. \
	\
	Can also create sbatch job file. \
	\
    Example (of the input list): \
        CCGO_800014 150223_OPTIMUS_C6HGHANXX.8.9477663 \
        CCGO_800015 150306_YOSHI_C6HDMANXX.5.9477645 \
        CCGO_800015 150223_OPTIMUS_C6HGHANXX.8.9477645") 
parser.add_argument('-j', '--job_name', help =
    "Give name for a GATK job sbatch job file. \
	(It will also be used for the sh script\
	for the scp transfer) \
	\
    The sbatch job can be invoked by a shell wrapper script to \
    continue with the realigned and merged BAM file \
    created by this script and further process and call \
    a GATK GVCF file.")
parser.add_argument('-b', '--exome_target_bed_file', help =
    "Give full path for exome target bed file used \
    for the capture process. They should be located \
    in biowulf2:/data/mcgaugheyd/genomes/1000G_phase2_GRCh37/converted_exome_bait_beds/")

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

# Grab all samples
samples = []
[samples.append(k) for k,v in sample_laneBam.items()]
samples = list(set(samples))
samples.sort()

# Pull names and create files
sbatch_file_name = args.job_name + '.sh'
scp_file_name = args.job_name + ".scp.sh"
sbatch_commands = open(sbatch_file_name, 'w')
scp_commands = open(scp_file_name, 'w')
bed_path = args.exome_target_bed_file


# loop to create the sh script for scp 
# and the sbatch script to run GATK after realignment
for one_sample in samples:
	# make sure the directory doesn't already exist
	if not os.path.isdir(one_sample):
		mkdir_call = "mkdir " + one_sample
		subprocess.check_call(mkdir_call, shell=True)
    # create Trek directory structure
	base_dir = "/cluster/ifs/projects/solexa/reads/"	
	# loop through each laneBam and create scp command
	for laneBam in sample_laneBam[one_sample]:
		folder = laneBam.split('.')[0]
		full_dir = base_dir + folder + '/' + laneBam + '*'
		# create scp command
		scp_call = 'scp trek.nhgri.nih.gov:' + full_dir + ' ' + one_sample + '/\n'
		scp_commands.write(scp_call)
	# create sbatch command(s)
	sbatch_commands.write("#!/bin/bash\n")
	sbatch_call = 'sbatch -J ' + one_sample + 'GVCFcall --mem=20G --time=36:00:00 \
				  /home/mcgaugheyd/bin/exome_workflow_v02/process_and_callGVCF.sh ' + \
			      one_sample + '/' + one_sample + '.bwa-mem.b37.merged.bam ' + bed_path + '\n'
	sbatch_commands.write(sbatch_call)
	sbatch_commands.write("sleep 1\n")



