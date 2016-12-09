#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('file', help= 
	"Takes in NISC bam file, extracts fastsq, \
	realigns fastq with bwa-mem to 1000g phaseII hs37d5  with \
	read group info extracted from NISC bam \
	Requires bam2fastq/1.1.0 and samtools/1.2 to be loaded \
	\
	Invoke script with sbatch --mem=30G --cpus-per-task 10 \
	\
	Example (exome): \
   		sbatch --mem=30G --cpus-per-task 10 this_script.py A_BAM_FILE_001.bam")
args = parser.parse_args()
bamfile = args.file
bam_name = bamfile.split('.bam')[0]

# Extract fastq
print("Beginning fastq extraction")
bam2fastq_call = "bam2fastq -o " + bam_name + "#.fastq " + bamfile
subprocess.check_call(bam2fastq_call, shell=True)
print("Done")

# gzip
gzip_call_1 = "gzip " + bam_name + "_1.fastq"
gzip_call_2 = "gzip " + bam_name + "_2.fastq"
subprocess.check_call(gzip_call_1, shell=True)
subprocess.check_call(gzip_call_2, shell=True)

# Runs samtools view -h and captures output
samtools_input = 'samtools view -h ' + bamfile + '| head -n 100 | grep ^@PG'
samtools_view = (subprocess.check_output(samtools_input, shell=True)).decode('utf-8')

info = samtools_view.split('\t')

# Builds the new RG from file name and NISC provided info from their bam
# CL field has the path to the fastq file
cl = [x for x in info if "CL:" in x]
# example cl field: ['CL:novoalign -F STDFQ -t 400 -o SAM @RG\\tID:1\\tSM:unknown\\tLB:unknown\\tPU:160128_YOSHI_C7NNRANXX.7.11457312.000.1\\tPL:Illumina -c 7 --tags PU- LB- PG- -d /cluster/ifs/projects/Exomes/cliagt/resources/reference/fasta/hg19/hg19.ndx -f /cluster/ifs/projects/Exomes/WE_Small_Projects/hg19.new/Brooks_RetDeg/1045/160128_YOSHI_C7NNRANXX.7.11457312/160128_YOSHI_C7NNRANXX.7.11457312.000.1.fq /cluster/ifs/projects/Exomes/WE_Small_Projects/hg19.new/Brooks_RetDeg/1045/160128_YOSHI_C7NNRANXX.7.11457312/160128_YOSHI_C7NNRANXX.7.11457312.000.2.fq\n']
cl = cl[0].split(' ')
# example: 160128_YOSHI_C7NNRANXX.7.11457312.000.1.fq
fastq_file = cl[cl.index('-f') + 1].split('/')[-1]

# ID is basically the NISC file name for the fastq (minus the .fq at the end)
ID = 'ID:' + fastq_file[0:-3]
# SM is the bam file name or sample name. Needs to the same for each sample!
SM = 'SM:' + bam_name
# LB is the library
LB = 'LB:' + fastq_file.split('.')[2]
PL = 'PL:Illumina\\" \\'

Output = SM + '.bwa-mem.b37.bam'
# Joins all together
RG_core = '\\\\t'.join(['\\"\@RG',ID, SM, LB, PL])


# bwa alignment
print("BWA run beginning")
run_bwa =   ('/home/mcgaugheyd/bin/exome_workflow_v02/run_bwa-mem_hg37d5.sh ' +
            bam_name + '_1.fastq.gz ' + bam_name + '_2.fastq.gz ' +
            '\\@RG\\\\t' + ID + '\\\\t' + SM + '\\\\t' + LB + '\\\\t' + 'PL:Illumina ' +
            bam_name + '.bwa-mem.b37.bam')
print(run_bwa)
subprocess.check_call(run_bwa, shell=True)
print("All done!")
