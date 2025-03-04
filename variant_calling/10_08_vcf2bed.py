import pandas as pd
import numpy as np
import gzip
import sys


#fname = "snp_bcftools_FLT_snpEff.vcf.gz"
fname = sys.argv[1]


def convertGzVcf2Bed(gzipvcf):
    fh = gzip.open(gzipvcf, 'rt')
    linie = fh.readlines()
    k = [ele.split() for ele in linie]
    n = [ele for ele in k if '#' not in ele[0]]
    bed = [[ele[0], str(int(ele[1])-1), str(int(ele[1])-1)] for ele in n]
    dBed = pd.DataFrame(bed)
    return(dBed)

def convertVcf2Bed(vcf):
    fh = open(vcf, 'r')
    linie = fh.readlines()
    k = [ele.split() for ele in linie]
    n = [ele for ele in k if '#' not in ele[0]]
    bed = [[ele[0], str(int(ele[1])-1), str(int(ele[1])-1)] for ele in n]
    dBed = pd.DataFrame(bed)
    return(dBed)

if __name__ == '__main__':
    
   # Converting vcf to bed dataframe
    if fname.endswith('vcf.gz'):
        bed = convertGzVcf2Bed(fname)
        fnameID = fname.replace('.vcf.gz','')
        
    elif fname.endswith('vcf'):
        bed = convertVcf2Bed(fname)
        fnameID = fname.replace('.vcf','')

    # Saving dataframe
    bed.to_csv(fnameID + '.bed', sep='\t', header=False, index=False)
