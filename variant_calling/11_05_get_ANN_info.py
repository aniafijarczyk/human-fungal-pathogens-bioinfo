import pandas as pd
import numpy as np
import gzip
import io
import sys


#fname = "bcftools_FLT2_snpEff.INFO.gz"
fname = sys.argv[1]

# Global
bool_ANN = True



def get_unique(lista):
    N = []
    for ele in lista:
        if ele not in N:
            N.append(ele)
    return(N)

def readINFO(filepath):
    df_tab = pd.read_csv(filepath, sep="\t", header=None, names=["#CHROM","POS","REF","ALT","QUAL","TYPE","INFO"], compression="gzip")
    df_tab['INFO'] = df_tab['INFO'].str.replace('INDEL','TYPE=INDEL')
    return(df_tab)



def getINFO(table):

    ### Spreading INFO filed
    dinfo = table['INFO'].apply(lambda x: {a:b for a,b in [ele.split('=') for ele in x.split(';')]}).apply(pd.Series)
    dinfo['INFO_AF'] = dinfo.apply(lambda x: int(x['AC'].split(',')[0])/int(x['AN']), axis=1)

    return(dinfo)

def getANN(table):
    ann_tab = table['ANN'].values.tolist()
    AT = []
    for snp in ann_tab:
        # Combination
        annotations = snp.split(',')        
        combo_info = '/'.join(['|'.join([ele.split('|')[0], ele.split('|')[1], ele.split('|')[2], ele.split('|')[3]]
                                       ) for ele in annotations])
        # Unique annotations
        combo_ann = '/'.join(['|'.join(get_unique([ele.split('|')[1] for ele in annotations]))])
        
        # First annotation
        ann1 = snp.split(',')[0]
        inf1 = ann1.split('|')
        AT.append([inf1[0], inf1[1], inf1[2], inf1[3], inf1[7], inf1[8], inf1[9], inf1[10], combo_ann, combo_info])

    dAT = pd.DataFrame(AT, columns = ['ANN_mut','ANN_type','ANN_effect','ANN_feature',
                                      'ANN_feature_type','ANN_gt','ANN_nc','ANN_aa',
                                      'ANN_types_all', 'ANN_info_all'])
    return(dAT)




if __name__ == '__main__':
    
    # Convertig tab to nice dataframe
    df_vcf = readINFO(fname)
    
    # Basic SNP info
    dcombo = df_vcf.iloc[:,:6].reset_index(drop=True)

    # Getting INFO fields
    df_info = getINFO(df_vcf)
    
    # Getting ANN fields
    df_ann = getANN(df_info)
    dcombo = pd.concat([dcombo, df_ann], axis=1)

    # Saving output
    #fnameID = fname.split('/')[-1].replace('.vcf.gz','')
    fnameID = fname.replace('.INFO.gz','')
    #dcombo.to_csv('getVcfInfo.tab', sep='\t', header=True, index=False)
    dcombo.to_csv(fnameID + '_INFO_ANN.tab', sep='\t', header=True, index=False)
    
    
    
