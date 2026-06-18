#!/bin/bash
#SBATCH --gres=lscratch:50
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --time=1:0:0

#cut -f 2 *.ped > sample.txt
#while read -r sample; do sbatch $sample; done < sample.txt

module load R/3.6.3

sample=$1

Rscript /home/$USER/git/variant_prioritization/src/sortGeminiTSV_v1.R gemini_tsv/$sample.novogene2021.GRCh37.novogene.gemini.tsv /data/OGL/resources/OGLpanelGeneDxORcandidate.xlsx gemini_tsv/$sample.novogene2021.GRCh37.novogene.gemini.rearranged.tsv gemini_tsv_filtered/$sample.novogene2021.GRCh37.novogene.gemini.filtered.tsv $sample gemini_xlsx/$sample.novogene2021.GRCh37.novogene.gemini.filtered.xlsx 0.5 ../manta/manta.$sample.annotated.tsv ../scramble_anno/$sample.scramble.xlsx
