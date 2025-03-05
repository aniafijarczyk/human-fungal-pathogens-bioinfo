#!/bin/bash

NAME=bcftools_FLT2_snpEff
DIR=/mnt/c/Users/aniaf/Projects/PC/GAPP/DATA/AFR_prediction_albicans/10_variants_5sets
VCF=${DIR}/${NAME}.vcf.gz
GFF=/mnt/c/Users/aniaf/Projects/DATA/Candida_albicans_reference/C_albicans_SC5314_version_A22-s07-m01-r183_features.gff

#bcftools query -r Ca22chr1A_C_albicans_SC5314:100000-500000 -f '%CHROM\t%POS0\t%POS\n' ${VCF} > ${DIR}/bcftools_FLT_snpEff_test.bed

### 1. Get variant positions in bed
echo "1. Get variant positions in bed"
conda activate bcftools-1.17
#bcftools query -f '%CHROM\t%POS0\t%POS\n' ${VCF} > ${DIR}/${NAME}.bed
conda deactivate

### 2. gff to bed
echo "2. gff to bed"
conda activate mbio
#python 11_04_gff2bed.py ${GFF} ${DIR}
conda deactivate

### 3. Intersection of vcf with gff
echo "3. Intersection of vcf with gff"
conda activate bedtools
VCF_core=$(basename $VCF | sed 's/.vcf.gz//g')
GFF_core=$(basename $GFF | sed 's/.gff//g')
#sort -k1,1 -k2,2n ${DIR}/${VCF_core}.bed > ${DIR}/${VCF_core}_sorted.bed
#sort -k1,1 -k2,2n ${DIR}/${GFF_core}.bed > ${DIR}/${GFF_core}_sorted.bed
#bedtools intersect -wao -a ${DIR}/${VCF_core}_sorted.bed -b ${DIR}/${GFF_core}_sorted.bed > ${DIR}/${VCF_core}_gff.tab
conda deactivate

rm ${DIR}/${VCF_core}.bed
rm ${DIR}/${GFF_core}.bed
rm ${DIR}/${VCF_core}_sorted.bed
rm ${DIR}/${GFF_core}_sorted.bed
##rm ${DIR}/${VCF_core}_gff.tab


### 4. Getting info table with snpEff annotations
conda activate mbio
#python 11_05_get_ANN_info.py ${DIR}/${VCF_core}.INFO.gz
conda deactivate

### 5. Formatting table
conda activate mbio
python 11_06_format_gff.py ${DIR}/${VCF_core}_gff.tab ${DIR}/${VCF_core}_INFO_ANN.tab
conda deactivate


# Basic info (forst 8 columns)
#bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' ${VCF} > ${DIR}/bcftools_FLT_snpEff_INFO.tab

