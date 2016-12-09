This workflow is based on:  

bwa/0.7.12 on 1000G phase II GRCh37 genome
GATK/3.5-0 with Haplotype Caller and hard filters (using GATK and bcbio best practices for exome)

NISC_laneBam_bwaRealign_callGVCF.sh is the wrapper script. It is fed a matrix file*, the sbatch job name, and the location of the exome bait bed file (specific to certain library prep kits). 

The script calls ~/bin/exome_workflow_v02/create_scp_and_sbatch_jobs_for_NISC_laneBams.py which creates scripts that are executed to scp the NISC lane bam files from Trek (1), process thems with BWA (2), and creates the script to call raw genotypes (1, 3). 

1. create_scp_and_sbatch_jobs_for_NISC_laneBams.py
2. realign_NISC_laneBams_with_bwa.py
3. process_and_callGVCF.sh 

Afterwards we have GVCF (raw genonotype) files for each exome sample. The next step is to run these in together, giving GVCF_to_hardFilteredVCF.sh the GVCFs in a \n separated list of the GVCFs you want to call together. 

Then you have a single VCF containing the genotypes for your cohort. The next steps are to run the Gemini workflow to prioritize variants for examination. 

Notes:

1. /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/converted_exome_bait_beds/Nimblegen_exome_v3_UTR_EZ.b37.bed is the current exome target file for NISC exomes  (2016-08-01)
2. /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta is used as the reference
- https://www.biostars.org/p/73100/
- http://www.1000genomes.org/faq/which-reference-assembly-do-you-use

\* The matrix file contains the 'common' and the lane bam file name, which can be used to find the exact file on Trek

```bash
sqlite3 ~/git/NGS_db/NISC_laneBam.sqlite3.db "SELECT NISC_LaneBams.Sample, LaneBam_File FROM NISC_LaneBams INNER JOIN Sample_Info ON NISC_LaneBams.Sample=Sample_Info.Sample WHERE Sample_Info.Project='CCGO' AND Sample_Info.DateAdded='2016-08-04'" | sort | cut -d"|" -f1,2,3 --output-delimiter=' ' > run_matrix.txt

# then split the run_matrix by sample with ~/bin/split_matrix_by_col.py
split_matrix_by_col.py run_matrix.txt 
cat CCGO_800365.laneBam.matrix
  CCGO_800365 160614_YOSHI_C90KTANXX.5.12551313
  CCGO_800365 160614_YOSHI_C90KTANXX.6.12551313
  CCGO_800365 160614_YOSHI_C90KTANXX.7.12551313
  CCGO_800365 160614_YOSHI_C90KTANXX.8.12551313
  CCGO_800365 160621_OPTIMUS_C90KRANXX.1.12551313
  CCGO_800365 160621_OPTIMUS_C90KRANXX.2.12551313
  CCGO_800365 160621_OPTIMUS_C90KRANXX.3.12551313
  CCGO_800365 160621_OPTIMUS_C90KRANXX.4.12551313
  CCGO_800365 160621_OPTIMUS_C90KRANXX.6.12551313
  CCGO_800365 160630_OPTIMUS_C90LGANXX.1.12551313
```
