module load samtools/1.13
#--min-ac 1 was used to obtain the variants absent in either.
bcftools view -s D847_001,D847_002 -r chrX:71109000-71149000 --min-ac 1 -Ov  broad.SORTED.VT.VEP.VCFANNO.vcf.gz | grep -v "^##" > MED12.20k.tsv