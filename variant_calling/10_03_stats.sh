#!/bin/bash

NAME=${1}
DIR=10_variants

bcftools stats -s - ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}.stats
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%QUAL\t%INFO/DP%INFO/AN\t%INFO/AC\n' ${DIR}/${NAME}.vcf.gz > ${DIR}/${NAME}.INFO
gzip -f ${DIR}/${NAME}.INFO
