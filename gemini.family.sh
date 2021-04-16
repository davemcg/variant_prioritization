#!/bin/bash
#SBATCH --gres=lscratch:150
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --time=6:0:0

#cut -f 1 *.ped | sort | uniq -d > family.txt
#while read -r family; do sbatch $family; done < family.txt
#$1 family name
module load gemini/0.20.1 R/3.6.3

cp novogene2021.GRCh37.PED_novogene.gemini.db /lscratch/$SLURM_JOB_ID/.

WKDIR=/lscratch/$SLURM_JOB_ID/
LENIENT=Yes
time gemini de_novo -d 9 --min-gq 5 --lenient --filter "priority_score > 3" --families $1 /lscratch/$SLURM_JOB_ID/novogene2021.GRCh37.PED_novogene.gemini.db > $WKDIR/denovo.tsv
time gemini autosomal_dominant -d 9 --min-gq 5 --filter "priority_score > 3" --families $1 /lscratch/$SLURM_JOB_ID/novogene2021.GRCh37.PED_novogene.gemini.db > $WKDIR/ad.tsv
time gemini autosomal_recessive -d 9 --min-gq 5 --lenient --filter "priority_score > 3" --families $1 /lscratch/$SLURM_JOB_ID/novogene2021.GRCh37.PED_novogene.gemini.db > $WKDIR/ar.tsv
time gemini comp_hets -d 9 --min-gq 5 --gene-where "priority_score >= 5" --families $1 /lscratch/$SLURM_JOB_ID/novogene2021.GRCh37.PED_novogene.gemini.db > $WKDIR/comphets.tsv
time gemini x_linked_de_novo -d 9 --min-gq 5 --filter "priority_score > 3" --families $1 /lscratch/$SLURM_JOB_ID/novogene2021.GRCh37.PED_novogene.gemini.db > $WKDIR/xdenovo.tsv
time gemini x_linked_dominant -d 9 --min-gq 5 --filter "priority_score > 3" --families $1 /lscratch/$SLURM_JOB_ID/novogene2021.GRCh37.PED_novogene.gemini.db > $WKDIR/xd.tsv
time gemini x_linked_recessive -d 9 --min-gq 5 --filter "priority_score > 3" --families $1 /lscratch/$SLURM_JOB_ID/novogene2021.GRCh37.PED_novogene.gemini.db > $WKDIR/xr.tsv
time gemini mendel_errors -d 9 --min-gq 5 --lenient --filter "priority_score > 3" --only-affected --families $1 /lscratch/$SLURM_JOB_ID/novogene2021.GRCh37.PED_novogene.gemini.db > $WKDIR/mendel_errors.tsv

Rscript /home/$USER/git/variant_prioritization/src/sortGeminiFamily.R /data/OGL/resources/OGLpanelGeneDxORcandidate.xlsx 0.5 gemini_xlsx/$1.novogene2021.GRCh37.novogene.lenientYes.xlsx $1 $WKDIR/denovo.tsv $WKDIR/ad.tsv $WKDIR/ar.tsv $WKDIR/comphets.tsv $WKDIR/xdenovo.tsv $WKDIR/xd.tsv $WKDIR/xr.tsv $WKDIR/mendel_errors.tsv


