#!/bin/bash

NAME=${1}
DIR=10_variants

bcftools concat --naive -Oz -o ${DIR}/${NAME}.vcf.gz ${DIR}/bcftools_Ca22chr1A.vcf.gz ${DIR}/bcftools_Ca22chr2A.vcf.gz ${DIR}/bcftools_Ca22chr3A.vcf.gz ${DIR}/bcftools_Ca22chr4A.vcf.gz ${DIR}/bcftools_Ca22chr5A.vcf.gz ${DIR}/bcftools_Ca22chr6A.vcf.gz ${DIR}/bcftools_Ca22chr7A.vcf.gz ${DIR}/bcftools_Ca22chrRA.vcf.gz ${DIR}/bcftools_Ca22chrM.vcf.gz
tabix -f -p vcf ${DIR}/${NAME}.vcf.gz
