#!/bin/bash

BAMS=$1

### REMOVING SAMPLES
echo "REMOVING SAMPLES"
BAMS=bam_files.txt
EX=bam_files_to_exclude.txt
python 10_00_get_bams.py $BAMS $EX


### GENERATING CONFIG
echo "GENERATING CONFIG"
if [ -f "00_get_bams.txt" ]; then
python 10_01_generate_config.py 10_00_get_bams.txt
else
python 10_01_generate_config.py ${BAMS}
fi

CONFIG=`pwd`/config_variant_calling.txt
. $CONFIG

### RUNNING VARIANT CALLER
echo "RUNNING VARIANT CALLER"
mkdir -p 10_variants
source activate bcftools
##sh 10_02_bcftools.sh ${CONFIG}
# Here instead variant I do variant calling for each chromosome separately
# with run_bcftools_chrom.sh
sh 10_02_bcftools_chrom.sh bam_files.txt Ca22chr1A_C_albicans_SC5314
sh 10_02_bcftools_chrom.sh bam_files.txt Ca22chr2A_C_albicans_SC5314
sh 10_02_bcftools_chrom.sh bam_files.txt Ca22chr3A_C_albicans_SC5314
sh 10_02_bcftools_chrom.sh bam_files.txt Ca22chr4A_C_albicans_SC5314
sh 10_02_bcftools_chrom.sh bam_files.txt Ca22chr5A_C_albicans_SC5314
sh 10_02_bcftools_chrom.sh bam_files.txt Ca22chr6A_C_albicans_SC5314
sh 10_02_bcftools_chrom.sh bam_files.txt Ca22chrRA_C_albicans_SC5314
sh 10_02_bcftools_chrom.sh bam_files.txt Ca22chrM_C_albicans_SC5314

### VCF MERGING
echo "VCF MERGING"
sh 10_02_02_merging.sh "bcftools"

### CHECKING VCF
echo "CHECKING VCF"
sh 10_03_stats.sh "bcftools"
conda deactivate

source activate vcftools
sh 10_04_missing.sh "bcftools"
conda deactivate

### BASIC FILTER
echo "FILTERING"
source activate bcftools
sh 10_05_filter.sh ${CONFIG}
conda deactivate

### SNPEFF
echo "ANNOTATION WITH SNPEFF"
source activate snpeff-5.0
sh 10_06_snpeff.sh ${CONFIG}
conda deactivate

### ANNOTATION
echo "ANNOTATION WTH GFF"
source activate bcftools
sh 10_07_vcf2bed.sh ${CONFIG}
conda deactivate

source activate bedtools
sh 10_08_annotate.sh ${CONFIG}
conda deactivate

### STATS
echo "GETTING STATS"
source activate vcftools
sh 10_12_vcftools.sh "bcftools_FLT_snpEff"
conda deactivate

### PLOTS
echo "MAKING PLOTS"
mkdir -p ./10_overview_plots
source activate mbio
python 10_13_overview.py 10_variants/bcftools_FLT_snpEff_nice_table.tab
conda deactivate
