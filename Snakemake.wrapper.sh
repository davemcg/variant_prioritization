#!/bin/bash

# to run snakemake as batch job
# run in the data folder for this project

module load snakemake || exit 1
mkdir -p 00log

sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"


# if json given, then use it
if [ ! -z "$2" ]; then
	json="$2"
# otherwise use the default
else
	json="/home/$USER/git/variant_prioritization/src/cluster.json"
fi

snakemake -s /home/$USER/git/variant_prioritization/src/Snakefile \
-pr --local-cores 2 --jobs 1999 \
--configfile $1 \
--cluster-config $json \
--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
-k --restart-times 0 --resources res=1 --notemp

# --notemp Ignore temp() declaration;
# --dryrun 
# --unlock