# Genomic Variant Prioritization
Snakemae worklow post-genotype calling to prioritize disease-causing variants.

# Input
- VCF from GVCF_to_hardFilteredVCF.sh
- PED with samples in VCF
- List of families to process

# Set up
Copy [src/config_variant_prioritization.yaml]() to your local folder and edit the `ped` field to give a path to your ped file. You will also need to edit the `family_name` to instruct Snakemake which families (must match ped family field, column 1) to create reports from. You can either give one family like so:

- family_name: 'gupta_fam'

Or a list of families to process like so:

- family_name: ['gupta_fam', 'smith_fam', 'chan_fam']

# Run (in biowulf2)
sbatch --time=12:00:00 ~/git/variant_prioritization/Snakemake.wrapper.sh COPIED_OVER_YAML_FILE.yaml
