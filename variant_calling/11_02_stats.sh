#!/bin/bash

NAME=bcftools_FLT2_snpEff
DIR=/mnt/c/Users/aniaf/Projects/PC/GAPP/DATA/AFR_prediction_albicans/10_variants_5sets

conda activate bcftools-1.17
#tabix -f -p vcf ${DIR}/${NAME}.vcf.gz
#bcftools query -l ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}_samples.txt
#bcftools stats -s - ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}.stats


bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%QUAL\t%INFO\n' ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}.INFO
gzip -f ${DIR}/${NAME}.INFO
echo "GT"
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/INDEL[\t%GT]\n' ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}_GT.tab
echo "DP"
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/INDEL[\t%DP]\n' ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}_DP.tab
echo "GQ"
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/INDEL[\t%GQ]\n' ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}_GQ.tab

conda deactivate

conda activate vcftools
echo "VCFTOOLS"
#vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --missing-indv --out ${DIR}/${NAME}
#vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --missing-site --out ${DIR}/${NAME}
#vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --het --out ${DIR}/${NAME}
#vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --depth --out ${DIR}/${NAME}
conda deactivate


