#!/bin/bash

. ${1}

# - good quality (mapping and base qual)
# - decent coverage across samples
# - biallelic only
# - removing low frequency alleles

DIR=10_variants


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



### Other filters
##bcftools filter -e "AVG(FORMAT/DP)<10 || INFO/DP>20000 || QUAL<20 || MQ<40 " -Ou snp_bcftools.vcf.gz | bcftools view -m2 -Oz -o snp_bcftools_FLT1.vcf.gz -
