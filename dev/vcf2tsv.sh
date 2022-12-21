module load samtools/1.13

vcf=$1
tsv=$2 #eys.OGLanno.tsv

bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQ\t%INFO/priority_score\t%INFO/PrScore_intervar\t%INFO/clinvar_hgmd_score\t%INFO/splice_score\t%INFO/insilico_score\t%INFO/pmaxaf\t%INFO/SpliceAI\t%INFO/spliceai_maxscore\t%INFO/SpliceAImasked50\t%INFO/SpliceAImasked50max\t%INFO/squirls_interpretation\t%INFO/squirls_maxscore\t%INFO/squirls_score\t%INFO/grch37variant_id\t%INFO/CLNID\t%INFO/CLNALLELEID\t%INFO/CLNDN\t%INFO/CLNDISDB\t%INFO/CLNREVSTAT\t%INFO/CLNSIG\t%INFO/CLNSIGCONF\t%INFO/InterVar_and_Evidence\t%INFO/Ref_Gene\t%INFO/Func_refGeneWithVer\t%INFO/ExonicFunc_refGeneWithVer\t%INFO/gno2x_ac_all\t%INFO/gno2x_af_all\t%INFO/gno2x_an_all\t%INFO/gno2x_popmax\t%INFO/gno2x_af_popmax\t%INFO/gno2x_hom\t%INFO/gno2x_filter\t%INFO/gno3_ac_all\t%INFO/gno3_af_all\t%INFO/gno3_an_all\t%INFO/gno3_popmax\t%INFO/gno3_af_popmax\t%INFO/gno3_nhomalt\t%INFO/gno3_filter\t%INFO/hgmd_id\t%INFO/hgmd_class\t%INFO/hgmd_phen\t%INFO/HGMD_Overlap4aa\t%INFO/omim_Gene\t%INFO/omim_Phen\t%INFO/omim_Inheritance\t%INFO/ClinPred_score\t%INFO/MetaSVM_pred\t%INFO/REVEL_score\t%INFO/MPC_score\t%INFO/PrimateAI_score\t%INFO/ccr_pct\t%INFO/remm\t%INFO/fathmm_XF_coding_score\t%INFO/fathmm_xf_noncoding\t%INFO/mutscore\t%INFO/hmc_score\n' $vcf \
| sed -e '1 s/^# //' -e '1 s/\[[0-9]*\]//g' > $tsv

##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Codons|Amino_acids|Gene|SYMBOL|MANE_SELECT|Feature|EXON|INTRON|HGVSc|HGVSp|MAX_AF|MAX_AF_POPS|Protein_position|BIOTYPE|CANONICAL|DOMAINS|Existing_variation|CLIN_SIG|PICK|PUBMED|Phenotypes|SIFT|PolyPhen|CADD_RAW|CADD_PHRED|GeneSplicer|SpliceRegion|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|existing_InFrame_oORFs|existing_OutOfFrame_oORFs|existing_uORFs|five_prime_UTR_variant_annotation|five_prime_UTR_variant_consequence|Mastermind_counts|Mastermind_MMID3|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE">

#'CLNALLELEID\t%INFO/CLNDN\t%INFO/CLNDISDB\t%INFO/CLNREVSTAT\t%INFO/CLNSIG\t%INFO/CLNSIGCONF'
