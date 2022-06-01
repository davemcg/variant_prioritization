chunks=$1

python fasta_generate_regions.py /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa.fai --chunks $chunks > vcf_region_split_$chunks.txt
lineNo=$(awk -F: '$1 ~ "chrM" {print NR}' vcf_region_split_$chunks.txt)
head -n $(($lineNo - 1)) vcf_region_split_$chunks.txt > vcf_region_split_$(($chunks - 1))_coords.txt
echo "MT_contigs" | cat vcf_region_split_$(($chunks - 1))_coords.txt - > vcf_region_split_$((chunks))_coords.txt
rm vcf_region_split_$chunks.txt