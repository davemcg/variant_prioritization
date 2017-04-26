#!/bin/bash

module load gemini/0.19.0

gemini_db=$1
family_name=$2
output_name=$3
new_aaf=$4

if [ -z "$4" ]
	then
		echo 'Running with default cohort AF (<0.1)'
		~/git/variant_prioritization/src/query_gemini.py -d $gemini_db -f $family_name -o $output_name
fi
	echo 'Running with user given cohort AF filter value'
	~/git/variant_prioritization/src/query_gemini.py -d $gemini_db -f $family_name -o $output_name --af_change $new_aaf
