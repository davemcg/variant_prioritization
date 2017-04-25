# Genomic Variant Prioritization
Worklow post-genotype calling to prioritize disease-causing variants.

# Input
VCF from GVCF_to_hardFilteredVCF.sh

# Process
Annotate with vcf_to_annotated_vcf.sh

Then create gemini db with vcf2db on cyclops (trying to get this installed on biowulf2)

Note: If you are only processing a trio, then you need to modify query_gemini.py to NOT filter on AF < 0.1
