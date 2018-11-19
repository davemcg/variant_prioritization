#!/usr/bin/env python3
from cyvcf2 import VCF, Writer
import subprocess 

# http://www.internationalgenome.org/phase-3-structural-variant-dataset
# wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz" .

vcf = VCF('ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz')
# adjust the header to contain the new field
# the keys 'ID', 'Description', 'Type', and 'Number' are required.
vcf.add_info_to_header({'ID': 'num_hom_alt', 
	'Description': 'Number of samples that are homozygous alternative',
    'Type':'Integer', 'Number': '1'})
vcf.add_info_to_header({'ID': 'num_het', 
	'Description': 'Number of samples that are heterozygous',
    'Type':'Integer', 'Number': '1'})

# create a new vcf Writer using the input vcf as a template.
f = 'ALL.wgs.mergedSV.v8.20130502.svs.genotypes.counts.vcf'
w = Writer(f, vcf)

for v in vcf:
	num_hom_alt = (v.gt_types == 3).sum()
	num_het = (v.gt_types == 1).sum()
	v.INFO["num_hom_alt"] = str(num_hom_alt.item())
	v.INFO["num_het"] = str(num_het.item())
	w.write_record(v)

w.close(); vcf.close()
subprocess.call(['bgzip', f])
subprocess.call(['tabix', '-p', 'vcf', f + '.gz'])
