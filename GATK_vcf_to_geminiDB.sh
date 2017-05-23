#!/bin/bash

# Assumes a GATK processed VCF after GATK-recommended filtering
# Hard coded against grch37
module load vcftools
module load vcfanno/0.1.1
module load vcf2db/7dfc48a
VCF=$1
ped=$2
gemini_db_name=$3


#vt to "left align and trim alternative variants"
cat $VCF \
	| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
	| ~/git/vt/./vt decompose -s - \
	| ~/git/vt/./vt normalize -r /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta - \
	> /scratch/mcgaugheyd/${VCF%.gz}

# annotate with VEP
/home/mcgaugheyd/git/variant_prioritization/src/run_VEP.sh /scratch/mcgaugheyd/${VCF%.gz} GRCh37 $SLURM_CPUS_PER_TASK
echo 'VEP done'

# compress and index
bgzip /scratch/mcgaugheyd/${VCF%.vcf.gz}.VEP.GRCh37.vcf
tabix -p vcf /scratch/mcgaugheyd/${VCF%.vcf.gz}.VEP.GRCh37.vcf.gz
echo 'bgzip and tabix done'

# annotate with custom annotations
vcfanno -p $SLURM_CPUS_PER_TASK -lua /home/mcgaugheyd/git/variant_prioritization/src/vcfanno_custom.lua \
	/home/mcgaugheyd/git/variant_prioritization/src/vcfanno_exomes.conf \
	/scratch/mcgaugheyd/${VCF%.vcf.gz}.VEP.GRCh37.vcf.gz | \
	bgzip > /scratch/mcgaugheyd/${VCF%.vcf.gz}.VEP.GRCh37.anno.vcf.gz
tabix -p vcf /scratch/mcgaugheyd/${VCF%.vcf.gz}.VEP.GRCh37.anno.vcf.gz
echo 'vcfanno done'

# move out of tmp folder
mv /scratch/mcgaugheyd/${VCF%.vcf.gz}.VEP.GRCh37.anno.vcf.gz* . 

# create gemini db
vcf2db.py ${VCF%.vcf.gz}.VEP.GRCh37.anno.vcf.gz $ped $gemini_db_name
