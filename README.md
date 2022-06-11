


# Genomic Variant Prioritization
Snakemake workflow post-genotype calling to prioritize disease-causing variants on biowulf2.

# Quick-ish Test Start Using Demo vcf
- Log into your biowulf2 account.
- `sinteractive`
- `mkdir -p ~/git`
- `cd ~/git`
- `git clone https://github.com/NEI/OGL/variant_prioritization.git`(to be created)
- After NGS_genotype_calling, `cd prioritization`
- `cp ~/git/variant_prioritization/config_variant_prioritization.yaml .`
- `sbatch --time=4:0:0 ../Snakemake.wrapper.sh config_variant_prioritization.yaml`


# Input
- VCF from deepvariant/freebayes tested [NGS_genotype_calling](https://github.com/NEI/OGL/NGS_genotype_calling/)
  Has to be bgzipped.
- PED with the same set of samples in VCF. The samples in PED and VCF must match. PED file has to be "\t" delimited. If header in PED, it has to start with #.
- SampleID in fastq files and PED files CANNOT contain "-" or "_" if using default script creating metadata.csv file. Seems that sampleID with "-" willl be converted by Gemini to "_".
- "Default" Gemini quieries for indiviudal samples and families will be included.

# Set up
Copy [config_variant_prioritization.yaml](https://github.com/davemcg/variant_prioritization/blob/master/src/config_variant_prioritization.yaml) to your local folder and edit the `ped` field to give a path to your ped file. You will also need to edit the `family_name` to instruct Snakemake which families (must match ped family field, column 1) to create reports from. You can either give one family like so:

####- family_name: 'gupta_fam'  - if you leave this blank (`family_name: ''`) then only the GEMINI database will be created (no family reports) Or a list of families to process like so:- family_name: ['gupta_fam', 'smith_fam', 'chan_fam'] Alternatively, family_name will be generated from PED file by the pipeline

Finally edit [config_variant_prioritization.yaml](https://github.com/davemcg/variant_prioritization/blob/master/src/config_variant_prioritization.yaml) to put your vcf (bgzip'ed and tabix'ed) in. 
#log
After git commit, run git log | head -n 5 > /data/OGL/resources/variant_prioritization.git.log
This file will be copied to project folder in SnakeWrapper 
# Run (in biowulf2)

sbatch --time=12:00:00 ~/git/variant_prioritization/Snakemake.wrapper.sh COPIED_OVER_YAML_FILE.yaml

# Visualization
![](variant_prioritization_dag_ogl.svg)
