#!/bin/bash

. ${1}

CHROM=${2}
REF=${REF}
BAMS=${BAMS}
DIR=10_variants
CHROM_NAME=$(echo $CHROM | cut -d"_" -f1)

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



# --skip-any-set: Skip reads with any of the FLAG bits set;
# 1796:
# read unmapped (0x4)
# not primary alignment (0x100)
# read fails platform/vendor quality checks (0x200)
# read is PCR or optical duplicate (0x400)
