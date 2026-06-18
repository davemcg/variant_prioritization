#!/bin/bash
#SBATCH --gres=lscratch:10
#SBATCH --time=48:0:0

counter=0
while [ $counter -lt $1 ]
do
  echo "Counter is $counter"
  squeue -u $USER | grep "gpu" | grep " R " | sed -e 's/^[[:space:]]*//g' | cut -d " " -f 1 | awk '{print "newwall -j " $0 " --time=48:0:0"}' > gpu.running.jobs.newwall.sh
  bash gpu.running.jobs.newwall.sh
  sleep 7h
  sleep 30m
  ((counter++)) # Increments counter by 1
done

