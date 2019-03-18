


# Genomic Variant Prioritization
Snakemake workflow post-genotype calling to prioritize disease-causing variants on biowulf2.

# Quick-ish Test Start Using Demo vcf
- Log into your biowulf2 account.
- `sinteractive`
- `mkdir -p ~/R/3.5/library`
- `module load R/3.5.2`
- `R`
- `devtools::install_github('davemcg/see_gem', build_vignettes=T)`
- `# THE ABOVE MUST INSTALL WITHOUT ERROR. IF IT DOES FIGURE IT OUT / ASK ME TO HELP.`
- `q()`
- `cd ~/`
- `mkdir -p ~/git`
- `cd ~/git`
- `git clone https://github.com/davemcg/variant_prioritization.git`
- `sed -i 's/mcgaugheyd/YOUR_BIOWULF2_USERNAME/g' src/vcfanno_v3.conf`
- `cd variant_prioritization/tests`
- `sbatch ../Snakemake.wrapper.sh config_variant_prioritization.yaml`


# Input
- VCF from [NGS_genotype_calling](https://github.com/davemcg/NGS_genotype_calling/blob/master/GVCF_to_VCF_snakemake.wrapper.sh)
- PED with samples in VCF
- List of families to process

# Set up
Copy [src/config_variant_prioritization.yaml](https://github.com/davemcg/variant_prioritization/blob/master/src/config_variant_prioritization.yaml) to your local folder and edit the `ped` field to give a path to your ped file. You will also need to edit the `family_name` to instruct Snakemake which families (must match ped family field, column 1) to create reports from. You can either give one family like so:

- family_name: 'gupta_fam'
  - if you leave this blank (`family_name: ''`) then only the GEMINI database will be created (no family reports)

Or a list of families to process like so:

- family_name: ['gupta_fam', 'smith_fam', 'chan_fam']

Install [SeeGEM](https://github.com/davemcg/SeeGEM) in `R` on biowulf2 to produce the html report. 
  - `sinteractive`
  - `module load R`
  - `R`
  - `devtools::install_github('davemcg/see_gem', build_vignettes=T)`

Finally edit the first line of [src/config_variant_prioritization.yaml](https://github.com/davemcg/variant_prioritization/blob/master/src/config_variant_prioritization.yaml) to put your vcf (bgzip'ed and tabix'ed) in. 

# Run (in biowulf2)
sbatch --time=12:00:00 ~/git/variant_prioritization/Snakemake.wrapper.sh COPIED_OVER_YAML_FILE.yaml

# Visualization
![](variant_prioritization_dag.svg)
