#!/bin/bash

NAME=${1}
DIR=10_variants

vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --missing-indv --out ${DIR}/${NAME}
vcftools --gzvcf ${DIR}/${NAME}.vcf.gz --missing-site --out ${DIR}/${NAME}

