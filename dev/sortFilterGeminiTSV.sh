#!/bin/bash
#SBATCH --gres=lscratch:10
#SBATCH --partition=quick
#SBATCH --mem=8g
#SBATCH --time=0:30:00

#$1=sample
#$2=vcf.ped file stem
mkdir -p gemini_tsv_filtered gemini_xlsx
module load R/3.6.3
Rscript /home/$USER/git/variant_prioritization/src_hg38/sortGeminiTSV_v1.R gemini_tsv/$1.$2.gemini.tsv gemini_tsv/$1.$2.gemini.ref.tsv /data/OGL/resources/OGLpanelGeneDxORcandidate.xlsx gemini_tsv/$1.$2.gemini.rearranged.tsv gemini_tsv_filtered/$1.$2.gemini.filtered.tsv $1 gemini_xlsx/$1.$2.gemini.filtered.xlsx 1.1 ../manta/manta.$1.annotated.tsv /data/OGL/resources/manta/manta.OGL.freq.2022-09.tsv ../AutoMap/$1/$1.HomRegions.annot.tsv ../scramble_anno/$1.scramble.xlsx filePlaceholder config_variant_prioritization.yaml

#while read -r sample; do sbatch  ~/git/variant_prioritization/dev/sortFilterGeminiTSV.sh $sample lp22-09.nano; done < sampleList.txt
#??parallel -j 8 "bash ~/git/variant_prioritization/dev/sortFilterGeminiTSV.sh {} lp22-09.nano" < sampleList.txt 
