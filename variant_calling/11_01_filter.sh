#!/bin/bash

#. ${1}

# - good quality (mapping and base qual)

DIR=/mnt/c/Users/aniaf/Projects/PC/GAPP/DATA/AFR_prediction_albicans/10_variants_5sets/
EXCL=samples_to_remove.txt
MIN_ALL=2
min_GQ=20

#conda activate vcftools
#vcftools --gzvcf $DIR/bcftools_FLT_snpEff.vcf.gz --remove ${EXCL} --remove-filtered-all --minGQ ${min_GQ} --min-alleles ${MIN_ALL} --recode --recode-INFO-all --out ${DIR}/bcftools_FLT2_snpEff
#conda deactivate


#mv ${DIR}/bcftools_FLT2_snpEff.recode.vcf ${DIR}/bcftools_FLT2_snpEff.vcf
#vcftools --gzvcf ${DIR}/bcftools_FLT2_snpEff.vcf --missing-indv --out ${DIR}/bcftools_FLT2_snpEff
#vcftools --gzvcf ${DIR}/bcftools_FLT2_snpEff.vcf --missing-site --out ${DIR}/bcftools_FLT2_snpEff
#vcftools --gzvcf ${DIR}/${VCF_CORE}.vcf.gz --het --out ${DIR}/${VCF_CORE}_het.tab
#vcftools --gzvcf ${DIR}/${VCF_CORE}.vcf.gz --depth --out ${DIR}/${VCF_CORE}_mean_depth.tab


### Other filters
conda activate bcftools
echo "BGZIP"
#bgzip ${DIR}/bcftools_FLT2_snpEff.recode.vcf
echo "FILTERING"
#bcftools view -m2 -Oz -o ${DIR}/bcftools_FLT2_snpEff_with_dubliniensis.vcf.gz ${DIR}/bcftools_FLT2_snpEff.recode.vcf
#tabix -f -p vcf ${DIR}/bcftools_FLT2_snpEff_with_dubliniensis.vcf.gz

echo "FILTERING2"
# Keeping only samples from previous step - recode.vcf but without 2 dubliniensis samples BL015 and BL050
SAMP=samples_to_keep.txt
bcftools view -m2 -S $SAMP -Oz -o ${DIR}/bcftools_FLT2_snpEff.vcf.gz ${DIR}/bcftools_FLT2_snpEff.recode.vcf.gz
tabix -f -p vcf ${DIR}/bcftools_FLT2_snpEff.vcf.gz

conda deactivate





