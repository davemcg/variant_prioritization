#!/bin/bash
module load samtools
# dbNSFP after 3.0 (?) is on GRCh38 by default
# they give hg19 coords, though, so I can just build new tabix against positions 
#ln -s /fdb/dbNSFP2/3.5a/dbNSFP3.5a.txt.gz /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/vep_annotation/dbNSFP3.5a.txt.gz

#cd  /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/vep_annotation/
zcat dbNSFP3.5a.txt.gz | head -n1 > h
zcat dbNSFP3.5a.txt.gz | grep -v ^#chr | sort -S 64G -k8,8 -k9,9n - | cat h - | awk '$8!="." {print $0}' | gzip -@ 2 -c > dbNSFP3.5a.sort.txt.gz
rm h
tabix -s 8 -b 9 -e 9 dbNSFP3.5a.sort.txt.gz
