#!/bin/bash

#module load R/3.6.3
module load samtools/1.13

#get gene coordinate
rm -f /lscratch/$SLURM_JOBID/gene.nochr.bed
rm -f /lscratch/$SLURM_JOBID/gene.wchr.bed
for gene in $@; do
	grep "$gene" /data/OGL/resources/omim/genemap2.txt | grep -v "^#" | awk -F"\t" 'BEGIN{OFS="\t"} {print $1,$2-1000,$3+1000}' | sed 's/^chr//' >> /lscratch/$SLURM_JOBID/gene.nochr.bed
	grep "$gene" /data/OGL/resources/omim/genemap2.txt | grep -v "^#" | awk -F"\t" 'BEGIN{OFS="\t"} {print $1,$2-1000,$3+1000}' >> /lscratch/$SLURM_JOBID/gene.wchr.bed
done

#need have chr for hg38 pipeline
#QUAL score is needed for the intervar_edit step
#ID has to be present for the crossmap left_join in the priority score step

echo "HGMD_Pro_2023.3" > readme.txt
bcftools view --no-header -Ov /data/OGL/resources/HGMD/hgmd-download/2023.3/HGMD_Pro_2023.3_hg38.bgzf.vcf.gz -R /lscratch/$SLURM_JOBID/gene.nochr.bed | grep -v "<DEL>" | awk -F"\t" 'BEGIN{OFS="\t"} {print "chr"$1,$2,".",$4,$5,"30","PASS",".","GT:GQ:DP","0/1:50:100"}' - | cat ~/git/variant_prioritization/dev/vcf.header.hg38 - | bgzip -c > /lscratch/$SLURM_JOBID/gene.hgmd.vcf.gz
tabix -f -p vcf /lscratch/$SLURM_JOBID/gene.hgmd.vcf.gz
#hgmd ABCA4: 1732 variants including 6 <DEL> (big del SV) i
echo "The ClinVar version:" >> readme.txt
zcat /data/OGL/resources/clinvar/clinvar.vcf.gz | head -n 2 | tail -n 1 >> readme.txt
bcftools view --no-header -Ov /data/OGL/resources/clinvar/clinvar.vcf.gz -R /lscratch/$SLURM_JOBID/gene.nochr.bed | awk -F"\t" 'BEGIN{OFS="\t"} { if (length($4) > 999 || length($5) > 999 || $5 == ".") {next} else {print "chr"$1,$2,".",$4,$5,"30","PASS",".","GT:GQ:DP","0/1:50:100"}}' - | cat ~/git/variant_prioritization/dev/vcf.header.hg38 - | bgzip -c > /lscratch/$SLURM_JOBID/gene.clinvar.vcf.gz
tabix -f -p vcf /lscratch/$SLURM_JOBID/gene.clinvar.vcf.gz
#clinvar ABCA4: 2702 variants, clinvar could have variant ALT as ".", remove, also remove the REF or ALT fields that have more nt than 2x500 as set in spliceAI. NNNNNNN in ALT did not affect annotation so far.

# bcftools concat -a --output-type u /lscratch/$SLURM_JOBID/gene.clinvar.vcf.gz /lscratch/$SLURM_JOBID/gene.hgmd.vcf.gz | bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u - | bcftools norm -d exact --output-type u | bcftools annotate --set-id '%CHROM\:%POS%REF\>%ALT' -Ou | bcftools sort -T /lscratch/$SLURM_JOB_ID/ -m 28G -O z -o gene.hgmd.clinvar.vcf.gz
# tabix -f -p vcf gene.hgmd.clinvar.vcf.gz
#hgmd+clinvar ABCA4: 3425 variants

bcftools concat -a --output-type u /lscratch/$SLURM_JOBID/gene.clinvar.vcf.gz /lscratch/$SLURM_JOBID/gene.hgmd.vcf.gz | bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u - \
	| bcftools norm -d exact --output-type u | bcftools annotate --set-id '%CHROM\:%POS%REF\>%ALT' -Ou | bcftools sort -T /lscratch/$SLURM_JOB_ID/ -m 28G -O z -o gene.vcf.gz
#ABCA4: 1140289 variants

tabix -f -p vcf gene.vcf.gz

echo "The variant no. in the gene.vcf.gz file:" >> readme.txt
echo $(bcftools index -n gene.vcf.gz) >> readme.txt
