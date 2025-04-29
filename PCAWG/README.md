# PCAWG TCGA somatic mutations

Here we report the integrative analysis of 2,658 whole-cancer genomes and their matching normal tissues across 38 tumour types from the Pan-Cancer Analysis of Whole Genomes (PCAWG) Consortium of the International Cancer Genome Consortium (ICGC) and The Cancer Genome Atlas (TCGA)

## Download data

### Links

- https://www.nature.com/articles/s41586-020-1969-6#Sec3
- https://www.synapse.org/Synapse:syn11726601/files/
- https://docs.icgc-argo.org/docs/data-access/icgc-25k-data#open-release-data---object-bucket-details

### Commands

```bash
cd data

# download code for spectra calculation
aws s3 cp s3://icgc25k-open/PCAWG/mutational_signatures/Code/PCAWG7_Data_Preparation_Code/PCAWG7_data_preparation_version_1.5.zip . --endpoint-url https://object.genomeinformatics.org --no-sign-request

# download WES somatic mutations
aws s3 cp s3://icgc25k-open/PCAWG/mutational_signatures/Input_Data_PCAWG7_23K_Spectra_DB/vcf_like_simple_files/WES_Other.20180327.simple.gz . --endpoint-url https://object.genomeinformatics.org --no-sign-request

# download WGS somatic mutations
aws s3 cp s3://icgc25k-open/PCAWG/mutational_signatures/Input_Data_PCAWG7_23K_Spectra_DB/vcf_like_simple_files/WGS_Other.20180413.simple.gz . --endpoint-url https://object.genomeinformatics.org --no-sign-request
```

### VCF-like files headers

```
# simple file has 'disease id cohort assembly snp/indel chr start end ref var annotations'
disease_idx = 0 # cancer type
tumour_idx = 1
assembly_idx = 3
vartype_idx = 4
chrome_idx = 5
pos_idx = 6
refbase_idx = 8
varbase_idx = 9
annotations_idx = 10
strand_idx = -1
```

### vartype SBS 

- 'SNP' used by TCGA, 'single base substitution' used by ICGC
- 'SNV' used in CCA lit. paper

### Select only SNPs and SNVs in 3 most represented reference

```bash
cd data

grep -e "        SNP" -e "       SNV" WES_Other.20180327.simple > WES.SNP.simple
grep -e "        SNP" -e "       SNV" WGS_Other.20180413.simple > WGS.SNP.simple

grep GRCh37 WGS.SNP.simple > WGS.SNP.GRCh37.simple
grep GRCh38 WGS.SNP.simple > WGS.SNP.GRCh38.simple
grep hg19 WGS.SNP.simple > WGS.SNP.hg19.simple

grep GRCh37 WES.SNP.simple > WES.SNP.GRCh37.simple
grep GRCh38 WES.SNP.simple > WES.SNP.GRCh38.simple

rm WGS.SNP.simple WES.SNP.simple
```

## Annotate using VEP

```bash
cd data

# reformat to vcf
python3 simple2vcf.py WES.SNP.GRCh37.simple WES.SNP.GRCh37.vcf
python3 simple2vcf.py WES.SNP.GRCh38.simple WES.SNP.GRCh38.vcf

python3 simple2vcf.py WGS.SNP.GRCh37.simple WGS.SNP.GRCh37.vcf
python3 simple2vcf.py WGS.SNP.GRCh38.simple WGS.SNP.GRCh38.vcf
# python3 simple2vcf.py WGS.SNP.hg19.simple WGS.SNP.hg19.vcf

bash annotate_variants.sh
```

## Download reference and proteome

```bash
cd data

datasets download genome accession GCF_000001405.40 --include protein,genome,cds
```
