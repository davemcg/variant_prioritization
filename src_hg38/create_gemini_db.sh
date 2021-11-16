#!/bin/bash

module load gemini/0.19.0

VCF=$1
PED=$2
DBNAME=$3

gemini load --cores $SLURM_CPUS_PER_TASK -t VEP -v $VCF -p $PED $DBNAME
