# python3

import os
import subprocess
import gzip

os.chdir('/data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/')
subprocess.call('rsync -av 165.112.66.100:/Volumes/Arges/2016-2018_John_Bryan_Archive/macbook_archive/Documents/git/rob_exac_gnomad/data/gnomad_genome_level_SumStatAnnotated.csv .', shell = True)
#subprocess.call('wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.annotation.gtf.gz ../', shell = True)

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

# read in sigmaAF csv (John Bryan)
LoF_0001 = {}
LoF_01 = {}
missense_0001 = {}
missense_01 = {}

gnomad_sigmaAF = open('/data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/gnomad_genome_level_SumStatAnnotated.csv')
for line in gnomad_sigmaAF:
	#line = line.decode('utf-8')
	line = line.split(',')
	gene = line[2].replace('"', '').upper()
	AF = line[1]
	sumAF = line[4]
	variation = line[3].replace('"', '')
	if AF == '1e-04' and variation == 'LoF':
		LoF_0001[gene] = sumAF
	if AF == '0.01' and variation == 'LoF':
		LoF_01[gene] = sumAF
	if AF == '1e-04' and variation == 'missense_variant':
		missense_0001[gene] = sumAF
	if AF == '0.01' and variation == 'missense_variant':
		missense_01[gene] = sumAF

output = []
for k, v in LoF_01.items():
	try:
		line = '\t'.join(gtf_dict[k].split('_')) + '\t' + \
			k + '\t' + \
			LoF_0001[k] + '\t' + \
			LoF_01[k] +  '\t' + \
			missense_0001[k] + '\t' + \
			missense_01[k] 
		output.append(line)
	except:
		print(k + ' is missing in gtf')

out_file = open('/data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/sigmaAF.bed', 'w')
for line in output:
	out_file.write(line)
	out_file.write('\n')
	
subprocess.call('python3 ~/git/ChromosomeMappings/convert_notation.py -c ~/git/ChromosomeMappings/GRCh37_UCSC2ensembl.txt -f /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/sigmaAF.bed \
	| sort -k1,1 -k2,2n | bgzip -f > /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/sigmaAF.ensembl.bed.gz', shell = True)
subprocess.call('tabix -f -s 1 -b 2 -e 3 -0 /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gemini_annotation/sigmaAF.ensembl.bed.gz', shell = True)
	

