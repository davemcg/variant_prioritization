#!/bin/bash

# Assumes a GATK processed VCF after GATK-recommended filtering
# Hard coded against grch37
module load vcftools
module load gemini/0.19.0

VCF=$1
PED=$2
DBNAME=$3

mkdir tmp

#only keep AF > 0.25
#vcffilter -f "AF > 0.25" $VCF > tmp/$VCF.AFfiltered

#vt to "left align and trim alternative variants"
cat $VCF \
	| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
	| ~/git/vt/./vt decompose -s - \
	| ~/git/vt/./vt normalize -r /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta - \
	> tmp/${VCF%.gz}

# annotate with VEP
/home/mcgaugheyd/bin/run_VEP.sh tmp/${VCF%.gz} GRCh37 $SLURM_CPUS_PER_TASK

# compress and index
bgzip tmp/${VCF%.vcf.gz}.VEP.GRCh37.vcf
tabix -p vcf tmp/${VCF%.vcf.gz}.VEP.GRCh37.vcf.gz

# annotate with Rob's eye gene list overlap and exac gene scores
vcfanno -p $SLURM_CPUS_PER_TASK ~/bin/gemini_exomes.toml tmp/${VCF%.vcf.gz}.VEP.GRCh37.vcf.gz | bgzip > tmp/${VCF%.vcf.gz}.VEP.GRCh37.anno.vcf.gz
tabix -p vcf tmp/${VCF%.vcf.gz}.VEP.GRCh37.anno.vcf.gz

# move out of tmp folder
mv tmp/${VCF%.vcf.gz}.VEP.GRCh37.anno.vcf.* . 

# delete temp files
rm -rf tmp

######## Now create Gemini DB #############
gemini load --cores $SLURM_CPUS_PER_TASK -t VEP -v ${VCF%.vcf.gz}.VEP.GRCh37.anno.vcf.gz -p $PED $DBNAME

# add annotations to Gemini DB 
gemini annotate -f ${VCF%.vcf.gz}.VEP.GRCh37.anno.vcf.gz \
	-o uniq_list,uniq_list,uniq_list,uniq_list,uniq_list,uniq_list,uniq_list,uniq_list,uniq_list,uniq_list,uniq_list,uniq_list,uniq_list \
	-t text,text,text,text,text,text,text,text,text,text,text,text,text \
	-e Gene_EyeDiseaseClass,n_syn,adj_exp_syn,n_mis,adj_exp_mis,n_lof,adj_exp_lof,syn_z,mis_z,lof_z,pRecessive,pNull,pLI \
	$DBNAME
