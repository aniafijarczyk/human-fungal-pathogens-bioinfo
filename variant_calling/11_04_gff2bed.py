import pandas as pd
import numpy as np
import gzip
import sys

#fname = "/home/DATA/Candida_albicans_reference/C_albicans_SC5314_version_A22-s07-m01-r183_features_hapA.gff"
fname = sys.argv[1]
dir = sys.argv[2]

def convertGff2Bed(filename):
    gff = pd.read_csv(filename, sep="\t", comment='#', header=None)
    gff['ID'] = gff[8].apply(lambda x: x.split(';')[0].split('=')[1])
    gff['name'] = gff.apply(lambda x: x[2] + "=" + x['ID'], axis=1)
    gff['start'] = gff[3].apply(lambda x: x-1)
    return(gff)


if __name__ == '__main__':
    
    # Converting gff to bed dataframe
    bed = convertGff2Bed(fname)
    
    # Saving dataframe
    fnameID = fname.split('/')[-1].replace('.gff','').replace('.gff3','')
    bed.to_csv(dir + "/" + fnameID + '.bed', sep='\t', columns = [0, 'start', 4, 'name'], header=False, index=False)
