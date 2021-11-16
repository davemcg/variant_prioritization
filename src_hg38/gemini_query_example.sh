

module load gemini/0.20.1
module load R/3.6.3

for sample in D1218_01 D1218_02 D1218_03 D1243_01 D1361_01 D1361_02 D1376_01 D1376_02 D1376_03; do time gemini query -q "SELECT *, gt_types.$sample, gts.$sample, gt_phases.$sample, gt_depths.$sample, gt_alt_freqs.$sample, gt_quals.$sample FROM variants" --header --gt-filter "gt_types.$sample != HOM_REF AND gts.$sample != './.' " --region 3:192310933-194415600 novogene133.glnexus.GRCh37.PED_novogene.gemini.db > OPA1/gemini/$sample.novogene133.glnexus.GRCh37.novogene.gemini.tsv; Rscript /home/$USER/git/variant_prioritization/src/sortGeminiTSV_v1.R OPA1/gemini/$sample.novogene133.glnexus.GRCh37.novogene.gemini.tsv /data/OGL/resources/OGLpanelGeneDxORcandidate.xlsx OPA1/gemini/$sample.novogene133.glnexus.novogene.gemini.rearranged.tsv OPA1/gemini_filtered_ps-3_aaf10pct/$sample.novogene133.glnexus.novogene.gemini.filtered.tsv $sample OPA1/gemini_xlsx/$sample.novogene133.glnexus.novogene.gemini.filtered.xlsx 0.1 ../manta/manta.$sample.annotated.tsv ../scramble_anno/$sample.scramble.xlsx; rm OPA1/gemini/$sample.novogene133.glnexus.GRCh37.novogene.gemini.tsv; done

WKDIR=/lscratch/$SLURM_JOB_ID/
for family in D1218 D1361 D1376; do \
time gemini de_novo -d 5 --min-gq 5 --lenient --filter "priority_score >= -3 and chrom == "3" and start < 194415600 and start > 192310933" --families $family  novogene133.glnexus.GRCh37.PED_novogene.gemini.db > $WKDIR/denovo.tsv; \
time gemini autosomal_dominant -d 5 --min-gq 5 --filter "priority_score >= 0 and chrom == "3" and start < 194415600 and start > 192310933" --families $family --allow-unaffected novogene133.glnexus.GRCh37.PED_novogene.gemini.db > $WKDIR/ad.tsv; \
time gemini autosomal_recessive -d 5 --min-gq 5 --lenient --filter "priority_score >= -3 and chrom == "3" and start < 194415600 and start > 192310933" --families $family novogene133.glnexus.GRCh37.PED_novogene.gemini.db > $WKDIR/ar.tsv; \
time gemini comp_hets -d 5 --min-gq 5 --gene-where "priority_score >= 0 and chrom == "3" and start < 194415600 and start > 192310933" --families $family novogene133.glnexus.GRCh37.PED_novogene.gemini.db > $WKDIR/comphets.tsv; \
touch $WKDIR/xdenovo.tsv; \
touch $WKDIR/xd.tsv; \
touch $WKDIR/xr.tsv; \
time gemini mendel_errors -d 5 --min-gq 5 --lenient --filter "priority_score >= -3 and chrom == "3" and start < 194415600 and start > 192310933" --only-affected --families $family novogene133.glnexus.GRCh37.PED_novogene.gemini.db > $WKDIR/mendel_errors.tsv; \
Rscript /home/$USER/git/variant_prioritization/src/sortGeminiFamily.R /data/OGL/resources/OGLpanelGeneDxORcandidate.xlsx 0.1 OPA1/gemini_xlsx/$family.novogene133.glnexus.GRCh37.novogene.lenientYes.xlsx $family $WKDIR/denovo.tsv $WKDIR/ad.tsv $WKDIR/ar.tsv $WKDIR/comphets.tsv $WKDIR/xdenovo.tsv $WKDIR/xd.tsv $WKDIR/xr.tsv $WKDIR/mendel_errors.tsv; done

for sample in D636_001 D636_002 D1218_0{1..3} D1243_01 D1361_01{1..2} D1376_0{1..3}; do Rscript /home/$USER/git/variant_prioritization/src/sortGeminiTSV_v1.R gemini_tsv/$sample.novogene2021.GRCh37.novogene.gemini.tsv /data/OGL/resources/OGLpanelGeneDxORcandidate.xlsx gemini_tsv/$sample.novogene2021.GRCh37.novogene.gemini.rearranged.tsv gemini_tsv_filtered/$sample.novogene2021.GRCh37.novogene.gemini.filtered.tsv $sample gemini_xlsx/$sample.novogene2021.GRCh37.novogene.gemini.filtered.xlsx 0.8 ../manta/manta.$sample.annotated.tsv ../scramble_anno/$sample.scramble.xlsx; done

for sample in D636_003; do Rscript /home/$USER/git/variant_prioritization/src/sortGeminiTSV_v1.R gemini_tsv/$sample.novogene2021.GRCh37.novogene.gemini.tsv /data/OGL/resources/OGLpanelGeneDxORcandidate.xlsx gemini_tsv/$sample.novogene2021.GRCh37.novogene.gemini.rearranged.tsv gemini_tsv_filtered/$sample.novogene2021.GRCh37.novogene.gemini.filtered.tsv $sample gemini_xlsx/$sample.novogene2021.GRCh37.novogene.gemini.filtered.xlsx 0.8 ../manta/manta.$sample.annotated.tsv ../scramble_anno/$sample.scramble.xlsx; done
