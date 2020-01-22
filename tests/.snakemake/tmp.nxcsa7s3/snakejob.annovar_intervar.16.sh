#!/bin/sh
# properties = {"type": "single", "rule": "annovar_intervar", "local": false, "input": ["temp/vt.test_trio__1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y.vcf.gz", "temp/vt.test_trio__1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y.vcf.gz.tbi"], "output": ["temp/test_trio__1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y.avinput", "temp/test_trio__1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y.avinput.hg19_multianno.txt", "temp/test_trio__1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y.avinput.hg19_multianno.txt.intervar"], "wildcards": {"sample": "test_trio", "region": "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"}, "params": {}, "log": [], "threads": 2, "resources": {}, "jobid": 16, "cluster": {"partition": "norm", "time": "24:00:00", "mem": "16g", "output": "00log/annovar_intervar_region=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,sample=test_trio.out", "error": "00log/annovar_intervar_region=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,sample=test_trio.err", "extra": ""}}
cd /spin1/home/linux/guanb/git/variant_prioritization/tests && \
/usr/local/Anaconda/envs_app/snakemake/5.4.4/bin/python3.6 \
-m snakemake temp/test_trio__1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y.avinput.hg19_multianno.txt --snakefile /home/guanb/git/variant_prioritization/src/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /spin1/home/linux/guanb/git/variant_prioritization/tests/.snakemake/tmp.nxcsa7s3 temp/vt.test_trio__1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y.vcf.gz temp/vt.test_trio__1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y.vcf.gz.tbi --latency-wait 120 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
 --configfile /spin1/home/linux/guanb/git/variant_prioritization/tests/config_variant_prioritization.yaml -p --allowed-rules annovar_intervar --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/spin1/home/linux/guanb/git/variant_prioritization/tests/.snakemake/tmp.nxcsa7s3/16.jobfinished" || (touch "/spin1/home/linux/guanb/git/variant_prioritization/tests/.snakemake/tmp.nxcsa7s3/16.jobfailed"; exit 1)

