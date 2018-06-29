#/usr/local/Anaconda/envs/py3.5/bin/python3

import ast
import sys
import gzip

# key from https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
key = {1: 'Tss',
	   2: 'TssAFlnk',
	   3: 'TxFlnk',
	   4: 'Tx',
	   5: 'TxWk',
	   6: 'EnhG',
	   7: 'Enh',
	   8: 'ZNF',
	   9: 'Het',
	   10: 'TssBiv',
	   11: 'BivFlnk',
	   12: 'EnhBiv',
	   13: 'ReprPC',
	   14: 'ReprPCWk',
	   15: 'Quies'}

input_file = gzip.open(sys.argv[1], 'rb')
for line in input_file:
	line = line.decode('utf-8')
	sline = line.split('\t')
	chrom = sline[0]
	start = sline[1]
	end = sline[2]
	qcat_field = sline[3].split('qcat:')[1]
	qcat = ast.literal_eval(qcat_field[:-1])
	qcat_dict = {}
	for field in qcat:
		qcat_dict[key[int(field[1])]] = str(field[0])
	output = chrom + '\t' + start + '\t' + end + '\t' + \
		qcat_dict['Tss'] + '\t' + \
		qcat_dict['TssAFlnk'] + '\t' + \
		qcat_dict['TxFlnk'] + '\t' + \
		qcat_dict['Tx'] + '\t' + \
		qcat_dict['TxWk'] + '\t' + \
		qcat_dict['EnhG'] + '\t' + \
		qcat_dict['Enh'] + '\t' + \
		qcat_dict['ZNF'] + '\t' + \
		qcat_dict['Het'] + '\t' + \
		qcat_dict['TssBiv'] + '\t' + \
		qcat_dict['BivFlnk'] + '\t' + \
		qcat_dict['EnhBiv'] + '\t' + \
		qcat_dict['ReprPC'] + '\t' + \
		qcat_dict['ReprPCWk'] + '\t' + \
		qcat_dict['Quies'] 
	print(output)	

input_file.close()



