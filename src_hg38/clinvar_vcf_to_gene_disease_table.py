#!/usr/local/Anaconda/envs/py3.5/bin/python

import sys
import gzip
import subprocess

vcf = gzip.open(sys.argv[1],'rt')

def gtf_reader(gtf):
	gtf_bed = {}
	for line in gtf:
		if line[0] == '#':
			continue
		line = line.split('\t')
		info = line[8]
		chr, start, end = line[0], int(line[3]), int(line[4])
		try:
			gene = [g.split('"')[1].upper() for g in info.split(';') if 'gene_name' in g][0]
		except:	
			pass
		if gene not in gtf_bed:
			gtf_bed[gene] = chr, start, end
		else:
			old_start, old_end = gtf_bed[gene][1:3]
			new_start = min(old_start, start)
			new_end = max(old_end, end)
			gtf_bed[gene] = chr, new_start, new_end
	return(gtf_bed)



# hand corrections
gene_corr = {'CYTB': 'MT-CYB','COX1': 'PTGS1', 'COX2': 'PTGS2', 'ZAK' : 'AC013461.1', 'ZFP112': 'ZNF112'}

# roll through vcf and grab gene and collapse diasease by variant to the the gene level
gene_disease = {}
for line in vcf:
	if line[0] == '#':
		continue
	info = line.split('\t')[7]
	try:
		gene = [gi.split('=')[1].split(':')[0] for gi in info.split(';') if 'GENEINFO' in gi][0]
		disease = [d.split('=')[1] for d in info.split(';') if 'CLNDBN' in d][0]
	except:
		sys.stderr.write('Gene or disease missing in line: ' + info)
	if gene in gene_corr:
		gene = gene_corr[gene]	
	if gene not in gene_disease: 
		gene_disease[gene] = disease
	else:
		existing_disease = gene_disease[gene]
		updated_disease = existing_disease + '|' + disease
		disease_set = list(set(updated_disease.split('|')))
		disease_set = [x for x in disease_set if x not in ['.','not_provided', 'not_specified']]
		disease_set = '|'.join(disease_set)
		gene_disease[gene] = disease_set

# pull in gtf to create bed
# one line for each gene, using min and max start/stop coordinates
gtf_bed = {}
gtf = gzip.open(sys.argv[2],'rt')
gtf_bed = gtf_reader(gtf)
gtf_alt = gzip.open(sys.argv[3],'rt')
gtf_bed_alt = gtf_reader(gtf_alt)

# match up gene_disease against gene_coordinates
for key in gene_disease:
	u_key = key.upper()
	if u_key in gtf_bed:
		start= str(gtf_bed[u_key][1])
		end = str(gtf_bed[u_key][2])
		print(gtf_bed[u_key][0] + '\t' + start + '\t' + end  + '\t' + gene_disease[key])
	elif u_key in gtf_bed_alt:
		start= str(gtf_bed_alt[u_key][1])
		end = str(gtf_bed_alt[u_key][2])
		print(gtf_bed_alt[u_key][0] + '\t' + start + '\t' + end  + '\t' + gene_disease[key])
	else:
		u_key = 'MT-' + u_key
		try:
			start= str(gtf_bed[u_key][1])
			end = str(gtf_bed[u_key][2])
			print(gtf_bed[u_key][0] + '\t' + start + '\t' + end  + '\t' + gene_disease[key])
		except:
			sys.stderr.write(key + ' not in GTFs\n')
