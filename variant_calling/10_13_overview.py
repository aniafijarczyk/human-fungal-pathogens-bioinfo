import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

fname = sys.argv[1]


# GLOBAL
stats_all = ['INFO_AF', 'INFO_sumDP4', 
         'GT_N_nomissing','GT_F_nomissing',
         'GT_N_hets_nomissing', 'GT_F_hets_nomissing',
         'GT_Mean_DP','GT_Mean_DP_nomissing']

features_to_exclude_all = [
 'centromere',
 'long_terminal_repeat',
 'ncRNA',
 'pseudogene',
 'rRNA',
 'repeat_region',
 'snRNA',
 'snoRNA']

def readTable(filename):
    df = pd.read_csv(filename, sep='\t', header=0, low_memory=False)
    return(df)

def getFeatureList(columns):
    F = []
    for i in df.columns.values[2:]:
        if i != 'feature_annots':
            F.append(i)
        else:
            break
    return F

def convert2Long(dataframe, stats_cols, features_cols):
    df_ = dataframe.set_index(['#CHROM', 'POS'] + stats_cols)[features_cols].stack().reset_index()
    new_col = [i for i in df_.columns.values if str(i).startswith('level_')][0]
    df_ = df_[df_[0]].reset_index(drop=True)
    df_ = df_.rename(columns = {new_col : 'feature'})
    return df_

def plotHists(dataframe, column_grid, x, ncols, color, title):
    g = sns.FacetGrid(dataframe, col=column_grid, col_wrap=ncols, sharey=False)
    g.map(sns.histplot, x, color = color)
    g.fig.suptitle(title)
    plt.subplots_adjust(top=0.9)
    plt.savefig('10_overview_plots/Histograms_' + x + '.png')

def plotScatters(dataframe, column_grid, x, y, ncols, color, title):
    g = sns.FacetGrid(dataframe, col=column_grid, col_wrap=ncols, sharey=False)
    g.map_dataframe(sns.scatterplot, x = x, y = y, color = color)
    g.fig.suptitle(title)
    plt.subplots_adjust(top=0.9)
    plt.savefig('10_overview_plots/Scatterplots_' + x + '_' + y + '.png')

def plotHist(dataframe, x, color, title):
    fig, ax = plt.subplots(figsize=(10,4))
    ax = sns.histplot(data = dataframe, x = x, color = color, ax=ax)
    ax.set_title(title)
    plt.savefig("10_overview_plots/Histogram_" + x + ".png")   

    

if __name__ == '__main__':
    
    # Reading stats
    df = readTable(fname)
    
    # Get lists
    features = getFeatureList(df.columns.values)
    stats = [i for i in df.columns.values if i in stats_all]
    features_to_exclude = [i for i in df.columns.values if i in features_to_exclude_all]
    
    # Converting into long format
    df_long = convert2Long(df, stats, features)
    
    # PLOTS
    
    # Allele frequency
    plotHists(df_long, "feature", "INFO_AF", 4, "tab:blue", title = 'INFO_AF (AF after removing missing data = AC/AN)')
    
    # Summed DP
    plotHists(df_long, "feature", "INFO_sumDP4", 4, "tab:blue", title = 'sum DP4 (Summed DP over all individuals after removing bad reads)')
    
    # Mean DP of non-missing genotypes
    plotHists(df_long, "feature", "GT_Mean_DP_nomissing", 4, "tab:blue", title = 'Mean DP of non-missing genotypes')
    
    # Non-missing genotypes
    plotHists(df_long, "feature", "GT_N_nomissing", 4, "tab:blue", title = 'Number of non-missing genotypes')
    
    # Fraction of heterozygotes
    plotHists(df_long, "feature", "GT_F_hets_nomissing", 4, "tab:blue", title = 'Fraction of heterozygotes')
    
    # AF vs. mean DP
    plotScatters(df_long, "feature", "GT_N_nomissing", "GT_Mean_DP_nomissing", 4, "tab:blue", title = "AF vs. mean DP")
    
    # AF vs. hets
    plotScatters(df_long, "feature", "GT_N_nomissing", "GT_F_hets_nomissing", 4, "tab:blue", title = "AF vs. hets")
    
    # Genic SNPs only
    df_gene = df[df['gene']].reset_index(drop=True)
    df_gene_long = convert2Long(df_gene, stats, features)
    pd.DataFrame(df_gene_long['feature'].value_counts()).to_csv('10_overview_plots/table_genic_features_variant_count.txt', sep='\t', header=True, index=True)

    # Genic SNPs, but excluding repeats, pseudogenes, retrotransposons, and RNAs (except mRNA and tRNAs)
    df_gene_flt = df[(df['gene']) & (~df[features_to_exclude].any(axis=1))].reset_index(drop=True)
    df_gene_flt_long = convert2Long(df_gene_flt, stats, features)
    pd.DataFrame(df_gene_flt_long['feature'].value_counts()).to_csv('10_overview_plots/table_genic_filtered_features_variant_count.txt', sep='\t', header=True, index=True)  
    
    # AF
    plotHist(df_gene_flt_long, 'INFO_AF', "tab:blue", "Genic filtered variants - allele frequency spectrum")
    
    # summed DP
    plotHist(df_gene_flt_long, 'INFO_sumDP4', "tab:blue", "Genic filtered variants - summed read depth")
    
    # mean DP
    plotHist(df_gene_flt_long, 'GT_Mean_DP_nomissing', "tab:blue", "Genic filtered variants - mean read depth")
    
    # Number of samples
    plotHist(df_gene_flt_long, 'GT_N_nomissing', "tab:blue", "Genic filtered variants - number of samples")
    
    # Heterozygotes
    plotHist(df_gene_flt_long, 'GT_F_hets_nomissing', "tab:blue", "Genic filtered variants - fraction of heterozygotes")
    
    
