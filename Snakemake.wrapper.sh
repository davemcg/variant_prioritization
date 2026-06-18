#!/bin/bash
#SBATCH --gres=lscratch:200
#SBATCH --cpus-per-task=28
#SBATCH --mem=32g
#SBATCH --time=24:0:0

# to run snakemake as batch job
# run in the data folder for this project
# $1 - configfile
# $2 - --notemp --dryrun --unlock --rerun-triggers mtime
# $3 non-default json file # currently not used

set -e
module load $(grep "^snakemake_version:" $1 | head -n 1 | cut -d"'" -f 2) || exit 1
#snakemake/7.19.1 1/3/2025 updated this to config_generic.yaml
#7.19.1 works with InterVar, but does not work with crossmap/0.6.5 if region has ","
#7.7.0 does not have --rerun-triggers mtime option.
#previous version 5.24.1, intervar/2.1.3 does not work with snakemake/6.0.5 version.

mkdir -p 00log

sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"

if [[ $(grep "^genomeBuild" $1 | grep -i "GRCh38" | wc -l) < 1 ]]; then
	json="/home/$USER/git/variant_prioritization/src/cluster.json"
	snakefile="/home/$USER/git/variant_prioritization/src/Snakefile"
else
	json="/home/$USER/git/variant_prioritization/src_hg38/cluster.json"
	snakefile="/home/$USER/git/variant_prioritization/src_hg38/Snakefile"
fi

sed -i 's/\r$//' $(grep "^ped:" $1 | head -n 1 | cut -d"'" -f 2)
# add json file to config file if needed.
#if [ ! -z "$3" ]; then
#	json="$3"
# otherwise use the default
#else
#	json="/home/$USER/git/variant_prioritization/src_hg38/cluster.json"
#fi

WORK_DIR=$PWD
if (( $(echo $WORK_DIR | grep "/data/OGL" | wc -l) > 0 )); then
 find $WORK_DIR -type f ! -group OGL -print -exec chgrp OGL -- {} +
fi

check=$(echo $@ | grep "dryrun\|dry-run\|unlock\|touch" | wc -l)
if (( $check > 0 )); then
	echo "Argument contains unlock or dry-run"
else
	if [ -e ../*.yaml ]; then
		tail -n 1 ../*.yaml >> $WORK_DIR/$1
	fi
	cd ~/git/variant_prioritization
	echo "variant_prioritization.git version: '$(git describe --tags --abbrev=0)'" >> $WORK_DIR/$1
	cd $WORK_DIR
fi

snakemake -s $snakefile \
-pr --local-cores 4 --jobs 1999 --max-jobs-per-second 1 \
--cluster-config $json \
--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
-k --restart-times 1 --resources res=1 \
--configfile $@

if (( $(echo $WORK_DIR | grep "/data/OGL" | wc -l) > 0 )); then
 find $WORK_DIR -type f ! -group OGL -print -exec chgrp OGL -- {} +
fi
#$SLURM_JOB_CPUS_PER_NODE When local-cores set at 8 for genome 99 coordinates -- host machine used 24 cpus at the vcf_split step;
#Thus the local-cores flag did not  limit the number of local jobs.
#does res=1 limit the number of crossmap rule?
# --notemp Ignore temp() declaration;
# --dryrun
# --unlock
