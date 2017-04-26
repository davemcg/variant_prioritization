#!/bin/bash

module load gemini/0.19.0

gemini_db=$1
family_name=$2
output_name=$3
new_aaf=$4

if [ -z "$4" ]
	then
		~/git/variant_prioritization/src/query_gemini.py -d $gemini_db -f $family_name -o $output_name
fi
	~/git/variant_prioritization/src/query_gemini.py -d $gemini_db -f $family_name -o $output_name --af_change $new_af
