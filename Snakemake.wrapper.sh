#!/bin/bash

# to run snakemake as batch job
# run in the data folder for this project
# $2 - --notemp --dryrun --unlock
# $3 non-default json file

module load snakemake/5.7.4 || exit 1

cp /data/OGL/resources/variant_prioritization.git.log .
mkdir -p 00log

sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"


# if json given, then use it
if [ ! -z "$3" ]; then
	json="$3"
# otherwise use the default
else
	json="/home/$USER/git/variant_prioritization/src/cluster.json"
fi

snakemake -s /home/$USER/git/variant_prioritization/src/Snakefile \
-pr --local-cores 2 --jobs 1999 \
--configfile $1 \
--cluster-config $json \
--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
-k --restart-times 1 --resources res=1 $2

# --notemp Ignore temp() declaration;
# --dryrun 
# --unlock
