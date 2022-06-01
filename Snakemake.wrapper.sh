#!/bin/bash
#SBATCH --gres=lscratch:100
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g

# to run snakemake as batch job
# run in the data folder for this project
# $1 - configfile
# $2 - --notemp --dryrun --unlock
# $3 non-default json file

module load snakemake/5.24.1 || exit 1
#previous version 5.7.4/5.24.1, intervar/2.1.3 does not work with snakemake/6.0.5 version.

cp /data/OGL/resources/variant_prioritization.git.log .
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

# add json file to config file if needed.
#if [ ! -z "$3" ]; then
#	json="$3"
# otherwise use the default
#else
#	json="/home/$USER/git/variant_prioritization/src_hg38/cluster.json"
#fi

snakemake -s $snakefile \
-pr --local-cores 2 --jobs 1999 \
--configfile $1 \
--cluster-config $json \
--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
-k --restart-times 0 --resources res=1 $2

# --notemp Ignore temp() declaration;
# --dryrun
# --unlock
