### Steps for generating variant calls
1. Generating config
2. Variant calling for each chromosome (bcftools)
3. Merging vcf
4. Generating statistics
5. Basic variant filtering
6. Annotation with snpeff
7. Annotation with gff & creating the nice table
8. Generating statistics
9. Final filters
10. Generating statistics
11. Creating table with annotations
12. Getting additional stats

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
MAX_INFO_DP = N_SAMPLES * 200 # maximum sample depth, set automatically
MIN_AF = 0.0 # allele frequency of minor allele  - 0 means it is not set
BIALLELIC = False # biallelic only or multialleleic - biallelic used
THREADS=9
```

### 2. Variant calling for each chromosome (bcftools)
bcftools v1.17
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
### 3. Merging vcf
bcftools v1.17
```
sh 10_02_02_merging.sh "bcftools"
```
```
bcftools concat --naive -Oz -o ${DIR}/${NAME}.vcf.gz ${DIR}/bcftools_Ca22chr1A.vcf.gz ${DIR}/bcftools_Ca22chr2A.vcf.gz ${DIR}/bcftools_Ca22chr3A.vcf.gz ${DIR}/bcftools_Ca22chr4A.vcf.gz ${DIR}/bcftools_Ca22chr5A.vcf.gz ${DIR}/bcftools_Ca22chr6A.vcf.gz ${DIR}/bcftools_Ca22chr7A.vcf.gz ${DIR}/bcftools_Ca22chrRA.vcf.gz ${DIR}/bcftools_Ca22chrM.vcf.gz
tabix -f -p vcf ${DIR}/${NAME}.vcf.gz
```
### 4. Generating statistics
bcftools v1.17
vcftools v0.1.16
```
sh 10_03_stats.sh "bcftools"
sh 10_04_missing.sh "bcftools"
```
```
bcftools stats -s - ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}.stats
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%QUAL\t%INFO/DP%INFO/AN\t%INFO/AC\n' ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}.INFO
gzip -f ${DIR}/${NAME}.INFO
vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --missing-indv --out ${DIR}/${NAME}
vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --missing-site --out ${DIR}/${NAME}
```
### 5. Basic variant filtering
bcftools v1.17
```
sh 10_05_filter.sh ${CONFIG}
```

```
if [ $BIALLELIC = "True" ]; then
echo "Selecting biallelic SNPs"
bcftools filter -e "AVG(FORMAT/DP)<"${AVG_SAMPLE_DP}" || INFO/DP>"${MAX_INFO_DP}" || QUAL<"${MIN_QUAL}" || MQ<"${MIN_MQ}" || (AF < "${MIN_AF}") || (AC == 0)" \
-Ou ${DIR}/bcftools.vcf.gz | bcftools view -m2 -M2 -c1 -Oz -o ${DIR}/bcftools_FLT.vcf.gz -

else
echo "Selecting multiallelic variants"
bcftools filter -e "AVG(FORMAT/DP)<"${AVG_SAMPLE_DP}" || INFO/DP>"${MAX_INFO_DP}" || QUAL<"${MIN_QUAL}" || MQ<"${MIN_MQ}" || (AF < "${MIN_AF}") || (AC == 0)" \
-Ou ${DIR}/bcftools.vcf.gz | bcftools view -m2 -c1 -Oz -o ${DIR}/bcftools_FLT.vcf.gz -
fi

bcftools stats ${DIR}/bcftools_FLT.vcf.gz > ${DIR}/bcftools_FLT.stats
tabix -f -p vcf ${DIR}/bcftools_FLT.vcf.gz
```
### 6. Annotation with snpeff
snpeff v5.0
```
sh 10_06_snpeff.sh ${CONFIG}
```
```
snpEff -c $SNPEFF_CONFIG $SNPEFF_GENOME ${DIR}/bcftools_FLT.vcf.gz > ${DIR}/bcftools_FLT_snpEff.vcf
bgzip -f ${DIR}/bcftools_FLT_snpEff.vcf
tabix -f -p vcf ${DIR}/bcftools_FLT_snpEff.vcf.gz
```
### 7. Annotation with gff & creating the nice table
bcftools v1.17
bedtools v2.31.0
```
sh 10_07_vcf2bed.sh ${CONFIG}
sh 10_08_annotate.sh ${CONFIG}
```
```
# Positions in bed
bcftools query -f '%CHROM\t%POS0\t%POS\n' ${VCF} > ${DIR}/bcftools_FLT_snpEff.bed

# Basic info (first 8 columns)
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' ${VCF} > ${DIR}/bcftools_FLT_snpEff_INFO.tab

# gff to bed
python 10_09_gff2bed.py ${GFF} ${DIR}

# Intersection of vcf with gff
VCF_core=$(basename $VCF | sed 's/.vcf.gz//g')
GFF_core=$(basename $GFF | sed 's/.gff//g')
sort -k1,1 -k2,2n ${DIR}/${VCF_core}.bed > ${DIR}/${VCF_core}_sorted.bed
sort -k1,1 -k2,2n ${DIR}/${GFF_core}.bed > ${DIR}/${GFF_core}_sorted.bed
bedtools intersect -wao -a ${DIR}/${VCF_core}_sorted.bed -b ${DIR}/${GFF_core}_sorted.bed > ${DIR}/${VCF_core}_gff.tab

# Get table with info
python 10_10_get_vcf_info.py ${DIR}/${VCF_core}_INFO.tab

# Format table with bedtools output and combine with info table
python 10_11_bedtools_nice_table.py ${DIR}/${VCF_core}_gff.tab ${DIR}/${VCF_core}_INFO_EXT.tab
```
### 8. Generating statistics & plots
vcftools v0.1.16
```
vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --missing-indv --out ${DIR}/${NAME}
vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --missing-site --out ${DIR}/${NAME}

python 10_13_overview.py 10_variants/bcftools_FLT_snpEff_nice_table.tab
```
### 9. Final filters
vcftools v0.1.16  
bcftools v1.17

```
sh 11_01_filter.sh
```

Masking genotypes with GQ < 20, removing samples of poor quality, and removing variants with less than 2 alleles
```
vcftools --gzvcf $DIR/bcftools_FLT_snpEff.vcf.gz --remove ${EXCL} --remove-filtered-all --minGQ ${min_GQ} --min-alleles ${MIN_ALL} --recode --recode-INFO-all --out ${DIR}/bcftools_FLT2_snpEff
bgzip ${DIR}/bcftools_FLT2_snpEff.recode.vcf
bcftools view -m2 -Oz -o ${DIR}/bcftools_FLT2_snpEff_with_dubliniensis.vcf.gz ${DIR}/bcftools_FLT2_snpEff.recode.vcf
tabix -f -p vcf ${DIR}/bcftools_FLT2_snpEff_with_dubliniensis.vcf.gz
```
Parameters set:
min_GQ=20  
MIN_ALL=2  
EXCL # list of samples to exclude  

Another variant set without the two outgroup samples
```
SAMP=samples_to_keep.txt
bcftools view -m2 -S $SAMP -Oz -o ${DIR}/bcftools_FLT2_snpEff.vcf.gz ${DIR}/bcftools_FLT2_snpEff.recode.vcf.gz
tabix -f -p vcf ${DIR}/bcftools_FLT2_snpEff.vcf.gz
```

### 10. Generating statistics
vcftools v0.1.16  
bcftools v1.17
```
sh 11_02_stats.sh
```

```
bcftools query -l ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}_samples.txt
bcftools stats -s - ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}.stats
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%QUAL\t%INFO\n' ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}.INFO
gzip -f ${DIR}/${NAME}.INFO
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/INDEL[\t%GT]\n' ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}_GT.tab
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/INDEL[\t%DP]\n' ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}_DP.tab
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/INDEL[\t%GQ]\n' ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}_GQ.tab
vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --missing-indv --out ${DIR}/${NAME}
vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --missing-site --out ${DIR}/${NAME}
vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --het --out ${DIR}/${NAME}
vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --depth --out ${DIR}/${NAME}
```

### 11. Creating table with annotations
bcftools v1.17
bedtools v2.31.0
```
sh 11_03_get_gff_annotation.sh
```

```
echo "1. Get variant positions in bed"
bcftools query -f '%CHROM\t%POS0\t%POS\n' ${VCF} > ${DIR}/${NAME}.bed

echo "2. gff to bed"
python 11_04_gff2bed.py ${GFF} ${DIR}

echo "3. Intersection of vcf with gff"
conda activate bedtools
VCF_core=$(basename $VCF | sed 's/.vcf.gz//g')
GFF_core=$(basename $GFF | sed 's/.gff//g')
sort -k1,1 -k2,2n ${DIR}/${VCF_core}.bed > ${DIR}/${VCF_core}_sorted.bed
sort -k1,1 -k2,2n ${DIR}/${GFF_core}.bed > ${DIR}/${GFF_core}_sorted.bed
bedtools intersect -wao -a ${DIR}/${VCF_core}_sorted.bed -b ${DIR}/${GFF_core}_sorted.bed > ${DIR}/${VCF_core}_gff.tab

echo "4. Getting info table with snpEff annotations"
python 11_05_get_ANN_info.py ${DIR}/${VCF_core}.INFO.gz

echo "5. Formatting table"
python 11_06_format_gff.py ${DIR}/${VCF_core}_gff.tab ${DIR}/${VCF_core}_INFO_ANN.tab
```

### 12. Getting additional stats
Such as the proportion of heterozygotes and mean DP and combining them to form the nice table
```
11_07_get_DP.ipynb
11_08_get_GT.ipynb
11_09_combine_tables.ipynb
```

