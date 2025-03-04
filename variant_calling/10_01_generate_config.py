import glob
import re
import sys


# txt file with bams
#PATH = 'bam_files.txt'
PATH = sys.argv[1]

    
def readTxtColumn(filename, n_column):
    fh = open(filename, 'r')
    linie = fh.readlines()
    k = [ele.split()[n_column] for ele in linie]
    return k
    
    
# Global variables
BAMS = PATH
PATH_REFERENCE = "/project/chlandry/projects/MutationAccum/Anna/2023_08_gapp_albicans/DATA/C_albicans_SC5314_version_A22-s07-m01-r183_chromosomes_hapA.fasta"
SNPEFF_CONFIG_PATH = "/home/anfij/anaconda3/envs/snpeff-5.0/share/snpeff-5.0-1/snpEff.config"
SNPEFF_GENOME_NAME = "SC5314"
PATH_GFF = "/project/chlandry/projects/MutationAccum/Anna/2023_08_gapp_albicans/DATA/C_albicans_SC5314_version_A22-s07-m01-r183_features_hapA.gff"
# Bcftools
PLOIDY = 2
# Filtering
N_SAMPLES = len(readTxtColumn(PATH, 0))
MIN_QUAL = 20
MIN_MQ = 40
AVG_SAMPLE_DP = 10
MAX_INFO_DP = N_SAMPLES * 200
MIN_AF = 0.0
BIALLELIC = False
THREADS=9


def writeConfig():
    
    config = """
BAMS=""" + PATH + """
# Global variables
REF=""" + PATH_REFERENCE + """
SNPEFF_CONFIG=""" + SNPEFF_CONFIG_PATH + """
SNPEFF_GENOME=""" + SNPEFF_GENOME_NAME + """
GFF=""" + PATH_GFF + """
# Bcftools
PLOIDY=""" + str(PLOIDY) + """
# Filtering
N_SAMPLES=""" + str(N_SAMPLES) + """
MIN_QUAL=""" + str(MIN_QUAL) + """
MIN_MQ=""" + str(MIN_MQ) + """
AVG_SAMPLE_DP=""" + str(AVG_SAMPLE_DP) + """
MAX_INFO_DP=""" + str(MAX_INFO_DP) + """
MIN_AF=""" + str(MIN_AF) + """
BIALLELIC=""" + str(BIALLELIC) + """
THREADS=""" + str(THREADS) + """
"""
    wh = open('config_variant_calling.txt', 'w')
    wh.write(config)
    wh.flush()
    wh.close()
    
if __name__ == '__main__':
    
    writeConfig()
