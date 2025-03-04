#!/bin/bash

. ${1}

DIR=10_variants

VCF=${DIR}/bcftools_FLT_snpEff.vcf.gz


#bcftools query -r Ca22chr1A_C_albicans_SC5314:100000-500000 -f '%CHROM\t%POS0\t%POS\n' ${VCF} > ${DIR}/bcftools_FLT_snpEff_test.bed

# Positions in bed
bcftools query -f '%CHROM\t%POS0\t%POS\n' ${VCF} > ${DIR}/bcftools_FLT_snpEff.bed


# Basic info (forst 8 columns)
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' ${VCF} > ${DIR}/bcftools_FLT_snpEff_INFO.tab
