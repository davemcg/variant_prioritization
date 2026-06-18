#!/bin/bash
#SBATCH --gres=lscratch:200
#SBATCH --cpus-per-task=56
#SBATCH --mem=128g
#SBATCH --time=2:0:0

#while read -r gene; do sbatch ~/git/variant_prioritization/geneSearch_v1.sh $gene; done < genelist.tsv
#geneList.tsv file may require Unix file EOL symbols.
#for gene in X Y Z; do bash ~/git/variant_prioritization/geneSearch_v1.sh $gene; done
#for gene in X Y Z; do sbatch ~/git/variant_prioritization/geneSearch_v1.sh $gene; done
#finished in ~40 min.

set -e

geneName="novelGene"
TIMESTAMP=$(date "+%Y%m%d-%H%M%S")
YEARSTAMP=$(date "+%Y%m%d")

mkdir /lscratch/$SLURM_JOB_ID/temp-$TIMESTAMP

module load R/4.3.0 parallel
export geneName TIMESTAMP YEARSTAMP SLURM_JOB_ID SLURM_CPUS_PER_TASK


##gzip the tsv files already run.
# find gemini_tsv_filtered/ -type f -name "*.tsv" | \
# parallel -j 8 --dry-run '
    # base=$(basename {});
    # prefix=${base%%.*};
    # gzip -c "{}" > "/data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/v2/${prefix}.tsv.gz"
# '

##gzip and change the file names of the existing files in the GeneSearch folder 
# find . -maxdepth 1 -type f -name "*.tsv"  | \
# parallel -j 8 '
    # base=$(basename {});
    # prefix=${base%%.*};
    # gzip "{}" && mv "{}".gz ${prefix}.tsv.gz && chgrp OGL ${prefix}.tsv.gz
# '

#exome
#declare -i probandNo
probandNo=$(find /data/OGL/resources/GeneSearch/exome/gemini_tsv_filtered/ -type f -name "*.tsv.gz" -printf '%f\n' | awk -F '[_x]' '{print $1}' | sort -u | wc -l) 
find /data/OGL/resources/GeneSearch/exome/gemini_tsv_filtered/ -type f -name "*.tsv*" \
 | parallel -j $SLURM_CPUS_PER_TASK -k '
 filename=$(basename {});
 sample=${filename%%.*};
 Rscript ~/git/variant_prioritization/dev/novel_geneSearch_exome.R {} $geneName $sample /lscratch/$SLURM_JOB_ID/temp-$TIMESTAMP/$sample.tsv
 '

# for file in /data/OGL/resources/GeneSearch/exome/gemini_tsv_filtered/*.tsv; do
# 	Rscript ~/git/variant_prioritization/dev/geneSearch_exome.R $file $geneName $(basename $file | cut -d. -f 1) /lscratch/$SLURM_JOB_ID/temp-$TIMESTAMP/$(basename $file);
# 	done
head -n 1 $(ls /lscratch/$SLURM_JOB_ID/temp-$TIMESTAMP/*.tsv | head -n 1) > /lscratch/$SLURM_JOB_ID/"$geneName".exome.tsv
for file in /lscratch/$SLURM_JOB_ID/temp-$TIMESTAMP/*.tsv; do tail -n +2 $file >> /lscratch/$SLURM_JOB_ID/"$geneName".exome.tsv; done

if [[ $(cat /lscratch/$SLURM_JOB_ID/"$geneName".exome.tsv | wc -l) = 1 ]]; 
then
 echo "empty exome search"
 touch $YEARSTAMP."$geneName".exome.xlsx
else 
 Rscript ~/git/variant_prioritization/dev/geneSearch_combine_sample.R /lscratch/$SLURM_JOB_ID/"$geneName".exome.tsv $probandNo $YEARSTAMP."$geneName".exome.xlsx
fi
rm /lscratch/$SLURM_JOB_ID/temp-$TIMESTAMP/*

#genome
probandNo=$(find /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/ -type f -name "*.tsv.gz" -printf '%f\n' | awk -F '[_x]' '{print $1}' | sort -u | wc -l)
find /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/ -name "*.tsv*" \
 | parallel -j $SLURM_CPUS_PER_TASK -k '
 filename=$(basename {});
 sample=${filename%%.*};
 Rscript ~/git/variant_prioritization/dev/novel_geneSearch.R {} $geneName /lscratch/$SLURM_JOB_ID/temp-$TIMESTAMP/$sample.tsv'

# for file in /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/*.tsv; do
# 	Rscript ~/git/variant_prioritization/dev/geneSearch.R $file $geneName /lscratch/$SLURM_JOB_ID/temp-$TIMESTAMP/$(basename $file);
# 	done
head -n 1 $(ls /lscratch/$SLURM_JOB_ID/temp-$TIMESTAMP/*.tsv | head -n 1) > /lscratch/$SLURM_JOB_ID/"$geneName".genome.tsv
for file in /lscratch/$SLURM_JOB_ID/temp-$TIMESTAMP/*.tsv; do tail -n +2 $file >> /lscratch/$SLURM_JOB_ID/"$geneName".genome.tsv; done

if [[ $(cat /lscratch/$SLURM_JOB_ID/"$geneName".genome.tsv | wc -l) = 1 ]]; 
then 
 echo "empty genome search"
 touch $YEARSTAMP."$geneName".genome.xlsx
else
 Rscript ~/git/variant_prioritization/dev/geneSearch_combine_sample.R /lscratch/$SLURM_JOB_ID/"$geneName".genome.tsv $probandNo $YEARSTAMP."$geneName".genome.xlsx
fi
#rm /lscratch/$SLURM_JOB_ID/temp-$TIMESTAMP/*
chgrp OGL $YEARSTAMP."$geneName".exome.xlsx $YEARSTAMP."$geneName".genome.xlsx

# ## if column names are identical
# #!/bin/bash
# geneName=$1
# mkdir -p geneSearch
# head -n 1 $(ls gemini_tsv_filtered/*.filtered.tsv | head -n 1) > geneSearch/"$geneName".tsv 
# for file in gemini_tsv_filtered/*.tsv; do grep -e $'^'$geneName$'\t' -e $'\t'$geneName$'\t' $file >> geneSearch/"$geneName".tsv; done
