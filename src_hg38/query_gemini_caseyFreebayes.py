#!/usr/local/Anaconda/envs/py3.4.3/bin/python

# Yes, this uses python3. 3.4.3, to be specific. Should work on other
# versions of python3, but I haven't tested 



import argparse
from argparse import RawTextHelpFormatter
import subprocess
import xlsxwriter
from collections import Counter
import datetime
import sys
#########PARSER##############
parser = argparse.ArgumentParser(description=\
    """
	Queries gemini (v0.18) database to identify variants matching 
	models of automosomal recessive, de novo, mendelian error, 
	compound hets and autosomal dominant. 
	
	Returns an xlsx file of the results.
	
	Input: 
		database to query 
		family to analyze (optional) 
		name for output xlsx file.
	
	Examples (no need for sbatch): 
		query_gemini.py --database CCG0.gemini.db --family CCG0_800042 CCGO_800042.variants.xlsx 
		query_gemini.py --database CCGO.gemini.db --family all.variants.xlsx""", 

	formatter_class=RawTextHelpFormatter)

parser.add_argument('-d','--database',required=True)
parser.add_argument('-f','--family')
parser.add_argument('-o','--output_name', required=True)
parser.add_argument('-l','--lenient', default='Yes', help="Use '-l No' to not use lenient settings on Autosome Dominant.")
#########CODE#############
def autosomal_recessive(db, family):
	filter = " --filter \"aaf_esp_all < 0.01 AND aaf_1kg_all < 0.01 AND aaf_exac_all < 0.01 AND (is_coding=1 OR is_splicing=1) \
				AND filter IS NULL\" --gt-pl-max 10 "
	if family=='-':
		ar_query = "gemini autosomal_recessive" + columns + db + " " + filter
	else:
		new_columns = columns.replace('*','family_id=' + '\'' + family +'\'')
		ar_query = "gemini autosomal_recessive" + new_columns + \
					"--families " + family + " " + db + " " + filter
	ar = subprocess.check_output(ar_query,shell=True).decode('utf-8')
	ar = ar.split('\n')
	return(ar,ar_query)

def de_novo(db, family):
	filter = " --filter \"aaf_esp_all < 0.005 AND aaf_1kg_all < 0.005 AND aaf_exac_all < 0.005 AND (is_coding=1 OR is_splicing=1) \
				AND filter IS NULL\" --gt-pl-max 10 "
	if family=="-":
		dn_query = "gemini de_novo" + columns + db + " " + filter
	else:
		new_columns = columns.replace('*','family_id=' + '\'' + family +'\'')
		dn_query = "gemini de_novo" + new_columns + \
					"--families " + family + " " + db + " " + filter
	dn = subprocess.check_output(dn_query,shell=True).decode('utf-8')	
	dn = dn.split('\n')
	return(dn, dn_query)

def mendel_errors(db, family):
	# gemini v0.18 has a bug with this call:
		# Can't parse by family
		# Hence my workaround
	filter = " --filter \"aaf_esp_all < 0.005 AND aaf_1kg_all < 0.005 AND aaf_exac_all < 0.005 AND (is_coding=1 OR is_splicing=1) \
				AND filter IS NULL\" --gt-pl-max 1 "
	me_query = "gemini mendel_errors" + columns + db + " " + filter 	
	me = subprocess.check_output(me_query,shell=True).decode('utf-8')
	me = me.split('\n')
	if family == '-':
		me_out = me
	# i.e. there are no mendelian errors
	if me == ['']:
		me_out = me
	else:
		# Get header in, unformatted (will happen later)
		me_out = []
		me_out.append(me[0])
		# find family_id index 
		header = me[0].split('\t')	
		family_id_index = header.index('family_id') 
		# filter for only the family we want
		for line in me:
			s_line = line.split('\t')
			if line and s_line[family_id_index]==family:
				me_out.append(line)
	return(me_out, me_query)

def comp_hets(db, family):
	# can't call exac numbers in v0.18 (reported bug, fixed in next release)
	columns = 	" --columns \"chrom, start, end, codon_change, aa_change, type, impact, \
			impact_severity, gene, clinvar_gene_phenotype, pfam_domain, vep_hgvsp, \
			max_aaf_all, aaf_1kg_all, aaf_exac_all, \
			geno2mp_hpo_ct, gerp_bp_score, polyphen_score, cadd_scaled, sift_pred, \
			sift_score, vep_maxEntScan, vep_grantham \" "
	
	filter = " --filter \"aaf_esp_all < 0.01 AND aaf_1kg_all < 0.01 AND aaf_exac_all < 0.01 AND (is_coding=1 OR is_splicing=1) \
				AND filter IS NULL\" --gt-pl-max 10 --max-priority 2 "
	if family == "-":
		ch_query = "gemini comp_hets" + columns + db + " " + filter
	else:
		ch_query = "gemini comp_hets" + columns + \
					"--families " + family + " " + db + " " + filter
	ch = subprocess.check_output(ch_query,shell=True).decode('utf-8')
	ch = ch.split('\n')
	####
	# reorder to put common comp_het genes (more than 4 variants) at bottom
	####
	# find position of the gene column
	gene_index = ch[0].split('\t').index('gene')
	# get counts for genes (last item in ch is blank)
	gene_counts = Counter([x.split('\t')[gene_index] for x in ch[:-1]])
	# id genes that appear more than 4 times
	common_genes = [x[0] for x in gene_counts.items() if x[1]>4]
	unique_ch = [x for x in ch[:-1] if x.split('\t')[gene_index] not in common_genes]
	common_ch = [x for x in ch[:-1] if x.split('\t')[gene_index] in common_genes]
	new_ch = unique_ch
	new_ch.append('\n')
	new_ch.append('Below are likely false positives (more than four \
				variants in a gene are unlikely to be a deleterious comp het)')
	new_ch.append('\n')
	new_ch.extend(common_ch)
	return(new_ch, ch_query)

def autosomal_dominant(db, family, lenient):
	filter = " --filter \"aaf_esp_all < 0.0001 AND aaf_1kg_all < 0.0001 AND aaf_exac_all < 0.0001 AND (is_coding=1 OR is_splicing=1) \
				AND filter IS NULL\" --gt-pl-max 10 "
	if family == "-":
		ad_query = "gemini autosomal_dominant" + columns + db + " " + filter
	if lenient == 'yes':
		new_columns = columns.replace('*','family_id=' + '\'' + family +'\'')
		ad_query = "gemini autosomal_dominant --lenient" + new_columns + \
					"--families " + family + " " + db + " " + filter
	else:
		new_columns = columns.replace('*','family_id=' + '\'' + family +'\'')
		ad_query = "gemini autosomal_dominant" + new_columns + \
					"--families " + family + " " + db + " " + filter
	ad = subprocess.check_output(ad_query,shell=True).decode('utf-8')
	ad = ad.split('\n')
	return(ad, ad_query)

def overview(db, queries):
	# pull useful info and  parameters from vcf header
	vcf_header_query = "gemini query --header -q \"SELECT * FROM vcf_header\" " + db
	vcf_header = subprocess.check_output(vcf_header_query,shell=True).decode('utf-8')
	vcf_header = vcf_header.split('\n')
	vcf_header_bits = []
	vcf_header_bits.extend([x for x in vcf_header if x.startswith("##FILTER")])
	vcf_header_bits.extend([x for x in vcf_header if x.startswith("##GATKCommandLine")])
	vcf_header_bits.extend([x for x in vcf_header if x.startswith("##reference")])
	vcf_header_bits.extend([x for x in vcf_header if x.startswith("##VEP")])
	vcf_header_bits.extend([x for x in vcf_header if x.startswith("#CHROM")])
	# summary stats, queries used, ped file
	stats_query = "gemini stats --gts-by-sample " + db
	stats = subprocess.check_output(stats_query,shell=True).decode('utf-8')
	stats = stats.split('\n')
	ped_query = "gemini query --header -q \"SELECT * FROM samples\" " + db
	ped = subprocess.check_output(ped_query,shell=True).decode('utf-8')
	ped = ped.split('\n')
	output = []
	output.append("Date and Time this file was generated")
	output.append(str(datetime.datetime.now().date()) + '\t' + str(datetime.datetime.now().time()))
	output.append('Overall genotypes by sample')
	output.extend(stats)
	output.append('PED information used for calls')
	output.append('Gender: 1=male, 2=female, 0=unknown')
	output.append('Phenotype: 1=unaffected, 2=affected, 0=unknown')
	output.append("Any column after 'phenotype' is not used in these queries")
	output.extend(ped)
	output.append('Gemini Queries')
	output.extend(queries)
	output.append('')
	output.append("Info on GATK commands and filtering, reference used, VEP version, samples in VCF")
	output.extend(vcf_header_bits)
	return(output)

def output_to_xlsx(data,sheet_name):
	worksheet = workbook.add_worksheet(sheet_name)
	row = 0
	col = 0
	# Handling for nothing found. Don't want anyone thinking a messup happened
	# Will print first bit of info if this logic screws up
	if len(data) < 2:
		worksheet.write(0,0, "No variants found")	
		worksheet.write(1,0, data[0])
	else:		
		for line in data:
			line = line.split('\t')
			for unit in line: 
				worksheet.write(row, col, unit)
				col += 1
			col = 0
			row += 1

def main():
	db = args.database
	if args.family:
		family = args.family
	else:
		family = '-'
	lenient = args.lenient
	lenient = lenient.lower()
#	if lenient != 'no' or lenient != 'yes':
#		print("-l --lenient must be 'Yes' or 'No'")
#		sys.exit()
	# output time
	ar, ar_query = autosomal_recessive(db, family)
	output_to_xlsx(ar, "Autosomal Recessive")	

	ch, ch_query = comp_hets(db, family)
	output_to_xlsx(ch, "Compound Hets")

	ad, ad_query = autosomal_dominant(db, family, lenient)
	output_to_xlsx(ad, "Autosomal Dominant")

	# get all queries in one list
	queries = []
	queries.append(ar_query.replace('\t',' '))#, queries.append(dn_query.replace('\t',' '))
#	queries.append(me_query.replace('\t',' ')), 
	queries.append(ch_query.replace('\t',' '))
	queries.append(ad_query.replace('\t',' '))

	# Create the info worksheet
	overview_info = overview(db, queries)
	output_to_xlsx(overview_info, "Info")
	workbook.close()
	


# global stuff
args = parser.parse_args()
workbook = xlsxwriter.Workbook(args.output_name)
columns = 	" --columns \"chrom, start, end, codon_change, aa_change, type, impact, \
			impact_severity, gene, clinvar_gene_phenotype, pfam_domain, vep_hgvsp, \
			max_aaf_all, aaf_1kg_all, aaf_exac_all, exac_num_hom_alt, exac_num_het, \
			geno2mp_hpo_ct, gerp_bp_score, polyphen_score, cadd_scaled, sift_pred, \
			sift_score, vep_maxEntScan, vep_grantham, (gt_ref_depths).(*), (gt_alt_depths).(*) \" "

# run it!
main()
