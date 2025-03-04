#!/bin/bash


. ${1}

DIR=10_variants

VCF=${DIR}/bcftools_FLT_snpEff.vcf.gz
GFF=${GFF}

### 1. vcf to bed
#python 10_08_vcf2bed.py ${VCF}

#bcftools query -f '%CHROM\t%POS\t%POS\n' ${VCF} > ${DIR}/bcftools_FLT_snpEff.bed


### 2. gff to bed
python 10_09_gff2bed.py ${GFF} ${DIR}

### 3. Intersection of vcf with gff
VCF_core=$(basename $VCF | sed 's/.vcf.gz//g')
GFF_core=$(basename $GFF | sed 's/.gff//g')
sort -k1,1 -k2,2n ${DIR}/${VCF_core}.bed > ${DIR}/${VCF_core}_sorted.bed
sort -k1,1 -k2,2n ${DIR}/${GFF_core}.bed > ${DIR}/${GFF_core}_sorted.bed
bedtools intersect -wao -a ${DIR}/${VCF_core}_sorted.bed -b ${DIR}/${GFF_core}_sorted.bed > ${DIR}/${VCF_core}_gff.tab

rm ${DIR}/${VCF_core}.bed
rm ${DIR}/${GFF_core}.bed
rm ${DIR}/${VCF_core}_sorted.bed
rm ${DIR}/${GFF_core}_sorted.bed

##bedtools intersect -wao -a ${VCF} -b ${GFF} > ${SNP_core}_gff_test.tab


### 4. Format table
# Get table with info
python 10_10_get_vcf_info.py ${DIR}/${VCF_core}_INFO.tab

# Format table with bedtools output and combine with info table
python 10_11_bedtools_nice_table.py ${DIR}/${VCF_core}_gff.tab ${DIR}/${VCF_core}_INFO_EXT.tab

rm ${DIR}/${VCF_core}_gff.tab
rm ${DIR}/${VCF_core}_INFO.tab
rm ${DIR}/${VCF_core}_INFO_EXT.tab
