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
bcftools view --no-header -Ov /fdb/spliceai/spliceai_scores.masked.snv.hg38.vcf.gz -R /lscratch/$SLURM_JOBID/gene.nochr.bed | awk -F"\t" 'BEGIN{OFS="\t"} {print "chr"$1,$2,".",$4,$5,"30","PASS",".","GT:GQ:DP","0/1:50:100"}' - | cat ~/git/variant_prioritization/dev/vcf.header.hg38 - | bgzip -c > /lscratch/$SLURM_JOBID/gene.snv.vcf.gz 
tabix -p vcf /lscratch/$SLURM_JOBID/gene.snv.vcf.gz 

#insertion of 1 nt and del of 1-4 nt
bcftools view --no-header -Ov /fdb/spliceai/spliceai_scores.masked.indel.hg38.vcf.gz -R /lscratch/$SLURM_JOBID/gene.nochr.bed | awk -F"\t" 'BEGIN{OFS="\t"} {print "chr"$1,$2,".",$4,$5,"30","PASS",".","GT:GQ:DP","0/1:50:100"}' - | cat ~/git/variant_prioritization/dev/vcf.header.hg38 - | bgzip -c > /lscratch/$SLURM_JOBID/gene.indel.vcf.gz
tabix -p vcf /lscratch/$SLURM_JOBID/gene.indel.vcf.gz

bcftools view --no-header -Ov /data/OGL/resources/gnomad/release-2.1.1/gnomad.exomes.r2.1.1.noVEP.sites.liftover_grch38.vcf.gz -R /lscratch/$SLURM_JOBID/gene.wchr.bed | awk -F"\t" 'BEGIN{OFS="\t"} {print $1,$2,".",$4,$5,"30","PASS",".","GT:GQ:DP","0/1:50:100"}' - | cat ~/git/variant_prioritization/dev/vcf.header.hg38 - | bgzip -c > /lscratch/$SLURM_JOBID/gene.gnomad.e2.vcf.gz
tabix -p vcf /lscratch/$SLURM_JOBID/gene.gnomad.e2.vcf.gz

bcftools view --no-header -Ov /data/OGL/resources/gnomad/release-3.1.2/gnomad.genomes.v3.1.2.selectedINFO.sites.vcf.gz -R /lscratch/$SLURM_JOBID/gene.wchr.bed | awk -F"\t" 'BEGIN{OFS="\t"} {print $1,$2,".",$4,$5,"30","PASS",".","GT:GQ:DP","0/1:50:100"}' - | cat ~/git/variant_prioritization/dev/vcf.header.hg38 - | bgzip -c > /lscratch/$SLURM_JOBID/gene.gnomad.g3.vcf.gz
tabix -p vcf /lscratch/$SLURM_JOBID/gene.gnomad.g3.vcf.gz

bcftools view --no-header -Ov /data/OGL/resources/HGMD/hgmd-download/2022.3/HGMD_Pro_2022.3_hg38.bgzf.vcf.gz -R /lscratch/$SLURM_JOBID/gene.nochr.bed | grep -v "<DEL>" | awk -F"\t" 'BEGIN{OFS="\t"} {print "chr"$1,$2,".",$4,$5,"30","PASS",".","GT:GQ:DP","0/1:50:100"}' - | cat ~/git/variant_prioritization/dev/vcf.header.hg38 - | bgzip -c > /lscratch/$SLURM_JOBID/gene.hgmd.vcf.gz
tabix -f -p vcf /lscratch/$SLURM_JOBID/gene.hgmd.vcf.gz
#hgmd ABCA4: 1732 variants including 6 <DEL> (big del SV) i

bcftools view --no-header -Ov /data/OGL/resources/clinvar/clinvar.vcf.gz -R /lscratch/$SLURM_JOBID/gene.nochr.bed | awk -F"\t" 'BEGIN{OFS="\t"} { if (length($4) > 999 || length($5) > 999 || $5 == ".") {next} else {print "chr"$1,$2,".",$4,$5,"30","PASS",".","GT:GQ:DP","0/1:50:100"}}' - | cat ~/git/variant_prioritization/dev/vcf.header.hg38 - | bgzip -c > /lscratch/$SLURM_JOBID/gene.clinvar.vcf.gz
tabix -f -p vcf /lscratch/$SLURM_JOBID/gene.clinvar.vcf.gz
#clinvar ABCA4: 2702 variants, clinvar could have variant ALT as ".", remove, also remove the REF or ALT fields that have more nt than 2x500 as set in spliceAI. NNNNNNN in ALT did not affect annotation so far.

bcftools concat -a --output-type u /lscratch/$SLURM_JOBID/gene.clinvar.vcf.gz /lscratch/$SLURM_JOBID/gene.hgmd.vcf.gz | bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u - | bcftools norm -d exact --output-type u | bcftools annotate --set-id '%CHROM\:%POS%REF\>%ALT' -Ou | bcftools sort -T /lscratch/$SLURM_JOB_ID/ -m 28G -O z -o gene.hgmd.clinvar.vcf.gz
tabix -f -p vcf gene.hgmd.clinvar.vcf.gz
#hgmd+clinvar ABCA4: 3425 variants

bcftools concat -a --output-type u /lscratch/$SLURM_JOBID/gene.snv.vcf.gz /lscratch/$SLURM_JOBID/gene.indel.vcf.gz /lscratch/$SLURM_JOBID/gene.gnomad.e2.vcf.gz /lscratch/$SLURM_JOBID/gene.gnomad.g3.vcf.gz | bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u - \
	| bcftools norm -d exact --output-type u | bcftools annotate --set-id '%CHROM\:%POS%REF\>%ALT' -Ou | bcftools sort -T /lscratch/$SLURM_JOB_ID/ -m 28G -O z -o gene.vcf.gz
#ABCA4: 1140289 variants

tabix -f -p vcf gene.vcf.gz


region_size=10000 #~10 pieces
awk -v region_size="$region_size" -F"\t" 'BEGIN{OFS="\t"} {start=$2; chunk=int(($3-$2)/region_size+1); size=int(($3-$2)/chunk)+1;  for(m=0; m<chunk; m++){start=$2+size*m; end=start+size-1; if (end < $3) {print $1":"start"-"end} else {print $1":"start"-"$3} } }' /lscratch/$SLURM_JOBID/gene.wchr.bed > abca4.10.region

#chunk=10
#awk -v chunk="$chunk" -F"\t" 'BEGIN{OFS="\t"} {start=$2; size=int(($3-$2)/chunk)+1; for(m=0; m<10; m++){start=$2+size*m; end=start+size-1; if (end < $3) {print $1":"start"-"end} else {print $1":"start"-"$3} } }' gene.wchr.bed

#ABCA4 SNV+INDEL+gnomAD annotate with current version of ClinVar
bcftools annotate --threads 8 -x INFO/CLNALLELEID,INFO/CLNDN,INFO/CLNDISDB,INFO/CLNREVSTAT,INFO/CLNSIG -Oz -o abca4.SORTED.VT.VEP.VCFANNO.vcf.gz gene.SORTED.VT.VEP.VCFANNO.vcf.gz 

#ABCA4
#Lines   total/split/realigned/skipped:  1445347/0/273027/0
#REF/ALT total/modified/added:   1445347/0/0
#Lines   total/split/realigned/skipped:  1445347/0/0/0
#make two vcf files for annotation
#spliceAI runs ~25k variants per hour on a GPU, thus taking ~22 if spliting to 2 files. Need to split to 10 instead (4 hours)

##reference=GRCh38
##ID=<Description="ClinVar Variation ID">
##INFO=<ID=AF_ESP,Number=1,Type=Float,Description="allele frequencies from GO-ESP">
##INFO=<ID=AF_EXAC,Number=1,Type=Float,Description="allele frequencies from ExAC">
##INFO=<ID=AF_TGP,Number=1,Type=Float,Description="allele frequencies from TGP">
##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
##INFO=<ID=CLNDNINCL,Number=.,Type=String,Description="For included Variant : ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
##INFO=<ID=CLNDISDB,Number=.,Type=String,Description="Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN">
##INFO=<ID=CLNDISDBINCL,Number=.,Type=String,Description="For included Variant: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN">
##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status for the Variation ID">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance for this single variant; multiple values are separated by a vertical bar">
##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="Conflicting clinical significance for this single variant; multiple values are separated by a vertical bar">
##INFO=<ID=CLNSIGINCL,Number=.,Type=String,Description="Clinical significance for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:clinical significance; multiple values are separated by a vertical bar">
##INFO=<ID=CLNVC,Number=1,Type=String,Description="Variant type">
##INFO=<ID=CLNVCSO,Number=1,Type=String,Description="Sequence Ontology id for variant type">
##INFO=<ID=CLNVI,Number=.,Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier">
##INFO=<ID=DBVARID,Number=.,Type=String,Description="nsv accessions from dbVar for the variant">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)">
##INFO=<ID=MC,Number=.,Type=String,Description="comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence">
##INFO=<ID=ORIGIN,Number=.,Type=String,Description="Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other">
##INFO=<ID=RS,Number=.,Type=String,Description="dbSNP ID (i.e. rs number)">
##INFO=<ID=SSR,Number=1,Type=Integer,Description="Variant Suspect Reason Codes. One or more of the following values may be added: 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other">