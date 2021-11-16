#!/bin/python3

import gzip
import subprocess

gtf = gzip.open('/data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gencode.v28lift37.annotation.gtf.gz')

# build gene <-> coord dictionary
gtf_dict = {}
for line in gtf:
    line = line.decode('utf-8')
    if line[0] == '#':
        continue
    if line.split()[2] == 'gene':
        gene_name = [x for x in line.split(';') if 'gene_name' in x][0].split('"')[1].upper()
        chrom = line.split()[0]
        start = line.split()[3]
        end = line.split()[4]
        gtf_dict[gene_name] = chrom + '_' + start + '_' + end


eyeIntegration = gzip.open('/home/mcgaugheyd/git/variant_prioritization/data/eyeIntegration_TPM.tsv.gz')

header_file = open('/data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/eyeIntegration_TPM.header', 'w')
out_file = open('/data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/eyeIntegration_TPM.tsv', 'w')

for line in eyeIntegration:
	try:
		line = line.decode('utf-8')
	except:
		print(line + ' decode fail!')
	if 'GeneName' in line:
		header = line
		out = "chr\tstart\tend\t" + line
		header_file.write(out)
		continue

	sline = line.split()
	try:
		coords = '\t'.join(gtf_dict[sline[0].upper()].split('_'))
	except:
		print("Can't match " + sline[0] + " to GTF")
		continue
	out = coords + '\t' + line
	out_file.write(out)

eyeIntegration.close()
out_file.close()
header_file.close()

subprocess.call('python3 ~/git/ChromosomeMappings/convert_notation.py -c ~/git/ChromosomeMappings/GRCh37_UCSC2ensembl.txt -f /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/eyeIntegration_TPM.tsv \
    | sort -k1,1 -k2,2n  > /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/eyeIntegration_TPM.sorted.tsv', shell = True)
subprocess.call('cat /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/eyeIntegration_TPM.header /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/eyeIntegration_TPM.sorted.tsv | bgzip -f > /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/eyeIntegration_TPM.sorted.tsv.gz', shell = True)
subprocess.call("tabix -f -S 1 -s 1 -b 2 -e 3 -0 /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/eyeIntegration_TPM.sorted.tsv.gz", shell = True)
#subprocess.call("rm /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/eyeIntegration_TPM.tsv", shell = True)
#subprocess.call("rm /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/eyeIntegration_TPM.sorted.tsv", shell = True)

