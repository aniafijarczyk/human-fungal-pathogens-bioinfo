import pandas as pd
import numpy as np
import sys


#fname_bedtools = "snp_bcftools_FLT_snpEff_gff.tab"
#fname_info = "INFO_ANN.tab"

fname_bedtools = sys.argv[1]
#fname_info = sys.argv[2]


def reformatGffTable(filename):
    
    df = pd.read_csv(filename, sep = "\t", header=None, 
                     names = ['chrom','start','stop','Rchrom','Rstart','Rstop', 'annotation','ovrl'])
    df = df[['chrom', 'start', 'stop', 'annotation']].reset_index(drop=True)
    df['feature'] = df['annotation'].apply(lambda x: x.split('=')[0] if x != '.' else 'no_feat')
    df['featureID'] = df['annotation'].apply(lambda x: x.split('=')[1] if x != '.' else 'NA')
    df['pos'] = df['start'].apply(lambda x: x+1)
    
    df_features = df.groupby(['chrom','pos','feature']).agg(count = ('start','count')).astype(bool).unstack().fillna(False).reset_index()
    df_features.columns = ['_'.join(col).replace('count_','').replace(
        'chrom_','chrom').replace('pos_','pos') for col in df_features.columns.values]

    df_annotations = df.groupby(['chrom','pos']).agg(features = ('annotation','unique')).reset_index()
    df_annotations['feature_annots'] = df_annotations['features'].apply(lambda x: '|'.join(x))
    df_annotations = df_annotations[['chrom','pos','feature_annots']].reset_index(drop=True)

    dm = pd.merge(df_features, df_annotations, on = ['chrom','pos'], how = 'left')
    dm = dm.rename(columns = {'chrom':'#CHROM', 'pos':'POS'})
    
    return(dm)

    
if __name__ == '__main__':
    
    # Convertig bedtools output to nice dataframe
    df_bedtools = reformatGffTable(fname_bedtools)
    
    # Reading info table
    #df_info = pd.read_csv(fname_info, sep='\t', header=0)
    
    # Combining two
    #dm = pd.merge(df_bedtools, df_info, on = ['#CHROM', 'POS'], how = 'left')
    
    # Saving output
    fname_ID = fname_bedtools.replace("_gff.tab","")
    df_bedtools.to_csv(fname_ID + '_gff_table.tab', sep='\t', header=True, index=False)
