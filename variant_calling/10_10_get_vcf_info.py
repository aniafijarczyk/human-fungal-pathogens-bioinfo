import pandas as pd
import numpy as np
import gzip
import io
import sys


#fname = "snp_bcftools_FLT_snpEff_INFO.tab"
fname = sys.argv[1]

# Global
bool_DP4 = True
bool_ANN = True

# This takes too much time
bool_GT = False


missing = ['.', './.', '.|.']
heterozygotes = ['0/1', '0/2', '1/2', '0/3', '1/3', '2/3',
                '0|1', '0|2', '1|2', '0|3', '1|3', '2|3',
                '1|0', '2|0', '2|1', '3|0', '3|1', '3|2']


def get_unique(lista):
    N = []
    for ele in lista:
        if ele not in N:
            N.append(ele)
    return(N)

def readINFO(filepath):
    df_tab = pd.read_csv(filepath, sep="\t", header=None, names=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"])
    df_tab['INFO'] = df_tab['INFO'].str.replace('INDEL','TYPE=INDEL')
    return(df_tab)

def readVcfGzip(filepath):
    fh = gzip.open(filepath,'rt')
    lines = [l for l in fh if not l.startswith('##')]
    df_vcf = pd.read_csv(io.StringIO(''.join(lines)),sep='\t')
    df_vcf['INFO'] = df_vcf['INFO'].str.replace('INDEL','TYPE=INDEL')
    return(df_vcf)

def getINFO(table, DP4):
    
    ### Substituting INDEL into TYPE=INDEL
    #table['INFO'] = table['INFO'].str.replace('INDEL','TYPE=INDEL')

    ### Spreading INFO filed
    dinfo = table['INFO'].apply(lambda x: {a:b for a,b in [ele.split('=') for ele in x.split(';')]}).apply(pd.Series)
    dinfo['INFO_AF'] = dinfo.apply(lambda x: int(x['AC'].split(',')[0])/int(x['AN']), axis=1)
    
    if DP4:
        dinfo['INFO_sumDP4'] = dinfo['DP4'].apply(lambda x: sum([int(i) for i in x.split(',')]))

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


def getGT(table):
    gt_format = table.iloc[0, 8]
    gt_index = [x for x,n in enumerate(gt_format.split(':')) if n == 'GT'][0]
    dp_index = [x for x,n in enumerate(gt_format.split(':')) if n == 'DP'][0]

    #gt_tab = table.iloc[:4, 9:].values.tolist()
    gt_tab = table.iloc[:, 9:].values.tolist()
    GT = []
    for snp in gt_tab:
        gt = [ele.split(':')[gt_index] for ele in snp]
        gt_nm = [g for g in gt if g not in missing]
        gt_het = [h for h in gt_nm if h in heterozygotes]
        dp = [int(ele.split(':')[dp_index]) for ele in snp]
        dp_nm = [int(ele.split(':')[dp_index]) for ele in snp if ele.split(':')[gt_index] not in missing]

        # N non missing genotypes
        n_nm = len(gt_nm)
        f_nm = len(gt_nm)/float(len(gt))

        # N heteroygotes
        n_nm_hets = len(gt_het)
        f_nm_hets = len(gt_het)/float(len(gt_nm))

        # Mean dp
        mean_dp = np.mean(dp)
        mean_dp_nm = np.mean(dp_nm)

        GT.append([n_nm, f_nm, n_nm_hets, f_nm_hets, mean_dp, mean_dp_nm])

    dGT = pd.DataFrame(GT, columns = ['GT_N_nomissing', 'GT_F_nomissing', 'GT_N_hets_nomissing',
                                      'GT_F_hets_nomissing', 
                                      'GT_Mean_DP', 'GT_Mean_DP_nomissing'])
    return(dGT)


if __name__ == '__main__':
    
    # Convertig tab to nice dataframe
    df_vcf = readINFO(fname)
    
    # Basic SNP info
    dcombo = df_vcf.iloc[:,:6].reset_index(drop=True)
    
    # Getting INFO fields
    df_info = getINFO(df_vcf, DP4=bool_DP4)
    dcombo = pd.concat([dcombo, df_info], axis=1)
    
    # Getting ANN fields
    if bool_ANN:
        df_ann = getANN(df_info)
        dcombo = pd.concat([dcombo, df_ann], axis=1)
    
    # Getting genotype info (GT and DP tags)
    if bool_GT:
        df_GT = getGT(df_vcf)
        dcombo = pd.concat([dcombo, df_GT], axis=1)

    # Saving output
    #fnameID = fname.split('/')[-1].replace('.vcf.gz','')
    fnameID = fname.replace('_INFO.tab','')
    #dcombo.to_csv('getVcfInfo.tab', sep='\t', header=True, index=False)
    dcombo.to_csv(fnameID + '_INFO_EXT.tab', sep='\t', header=True, index=False)
    
    
    
