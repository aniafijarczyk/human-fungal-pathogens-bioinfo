#!/bin/bash


. ${1}

DIR=10_variants

snpEff -c $SNPEFF_CONFIG $SNPEFF_GENOME ${DIR}/bcftools_FLT.vcf.gz > ${DIR}/bcftools_FLT_snpEff.vcf
bgzip -f ${DIR}/bcftools_FLT_snpEff.vcf
tabix -f -p vcf ${DIR}/bcftools_FLT_snpEff.vcf.gz

