
# http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff
# USE VEP for annotation


APPTAINER=/opt/apptainer/bin/apptainer
VEP_CONTAINER=/opt/tools/vep.sif

# dataset=WES
# assembly=GRCh37

for dataset in WES WGS
do 
for assembly in GRCh37 GRCh38
do


INPUT=./$dataset.SNP.$assembly.vcf
OUTPUT=./$dataset.SNP.$assembly.annotated.vcf

echo $INPUT $OUTPUT

$APPTAINER exec \
  --pwd /data \
  --bind $PWD:/data \
  /opt/tools/vep.sif \
  vep -i $INPUT -o $OUTPUT --vcf \
    --species homo_sapiens --assembly $assembly \
    --cache --cache_version 113 --coding_only

done
done

# gatk VariantsToTable \
#     -V cohort.merged.flt.snpgap.dnm_gatk2.annotated.vcf \
#     -O cohort.merged.flt.snpgap.dnm_gatk2.annotated.table \
#     -F CHROM -F POS -F REF -F ALT -F TYPE -F CSQ