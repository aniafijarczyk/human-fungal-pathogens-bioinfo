### Steps for generating variant calls
1. Generating config
2. Variant calling for each chromosome (bcftools)
3. Merging vcf
4. Generating statistics
5. Basic variant filtering
6. Annotation with snpeff
7. Annotation with gff
8. Creating the nice table


### 1. Generating config file
```
BAMS=bam_files.txt
EX=bam_files_to_exclude.txt
python 10_00_get_bams.py $BAMS $EX
python 10_01_generate_config.py 10_00_get_bams.txt
```
Variables set in the config:
```
BAMS = PATH
PATH_REFERENCE = # path to reference fasta file
SNPEFF_CONFIG_PATH = # path to snpEff.config file
SNPEFF_GENOME_NAME = # name of the genome in snpEff confog
PATH_GFF = # path to gff file
# Bcftools
PLOIDY = 2 # diploid in case of albicans
# Filtering
N_SAMPLES = len(readTxtColumn(PATH, 0)) # number of samples - calculated from list of bams 
MIN_QUAL = 20 # min genotype quality
MIN_MQ = 40 # minimum mapping quality
AVG_SAMPLE_DP = 10 # mean average depth across samples
MAX_INFO_DP = N_SAMPLES * 200 # maximum sample depth - this filter is not used
MIN_AF = 0.0 # allele frequency of minor allele  - 0 means it is not set
BIALLELIC = False # biallelic only or multialleleic - biallelic used
THREADS=9
```

### 2. Variant calling for each chromosome (bcftools)

```
run_bcftools_chrom.sh
```
```
if [ $PLOIDY = "1" ] ; then
echo "Ploidy is 1"
bcftools mpileup -C50 -f ${REF} -r ${CHROM} -min-MQ 4 -min-BQ 13 --skip-any-set 1796 -a FORMAT/AD,FORMAT/ADR,FORMAT/ADF,FORMAT/DP --threads 14 -Ou -b ${BAMS} | \
bcftools call -mv -f gq -ploidy 1 -Oz -o ${DIR}/bcftools_${CHROM_NAME}.vcf.gz -
tabix -f -p vcf ${DIR}/bcftools_${CHROM_NAME}.vcf.gz

elif [ $PLOIDY = "2" ] ; then
echo "Ploidy is 2"
bcftools mpileup -C50 -f ${REF} -r ${CHROM} -q 4 -Q 13 --skip-any-set 1796 -a FORMAT/AD,FORMAT/ADR,FORMAT/ADF,FORMAT/DP --threads 14 -Ou -b ${BAMS} | \
bcftools call -mv -f gq -Oz -o ${DIR}/bcftools_${CHROM_NAME}.vcf.gz -
tabix -f -p vcf ${DIR}/bcftools_${CHROM_NAME}.vcf.gz
fi
```
