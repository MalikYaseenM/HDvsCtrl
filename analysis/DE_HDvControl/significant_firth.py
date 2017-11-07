#!/usr/bin/env python
import pandas as pd
import sys
import re

def pd_tab(filename):
    """Read tab separated file"""
    df = pd.read_csv(filename, sep="\t")
    return df

def pd_csv(filename):
    """Read comma separated file"""
    df = pd.read_csv(filename, sep=",")
    return df

def main():
    if len(sys.argv) < 2:
        print("Usage: {0} <CSV File to be tab-delimited> <biomart text>".format(sys.arv[0]))
        sys.exit(1)
    # Read files
    firth = sys.argv[1]
    biomart = sys.argv[2]

    # Turn biomart file into DF and change column name
    ensemble = pd_tab(biomart)
    ensemble = ensemble.rename(columns = {'Gene stable ID':'gene_id'})

    ############################## FIRTH ####################################
    df = pd_csv(firth)
    # Rename first column to gene id
    df = df.rename(columns = {'Unnamed: 0':'gene_id'})
    # Take the gene id number before the period (get rid of the period)
    df['gene_id'] = df['gene_id'].str.split('.').str[0]
    # Drop other columns, only keep the counts
    col = list(df)
    cols = col[:1] + col[7:10]
    df.drop([col for col in df.columns if col not in cols],axis=1, inplace=True)

    # Merge with gene names from biomart
    df = pd.merge(df, ensemble, on='gene_id')

    # Re-index columns for easier read
    col = list(df)
    cols = col[:1] + col[-2:] + col[1:4]
    df = df[cols]

    # Find significant genes, sort in ascending padj order
    df = df.sort_values(['counts__padj'], ascending=[True])
    sig_fdr = df[df['counts__padj'] < 0.05]

    # Set new file names from firth file
    gene = re.sub('\_firth', '_firth_sig_gene', str(firth))
    counts = re.sub('\_firth', '_firth_sig_counts', str(firth))

    # Create new files with significant genes/counts
    sig_fdr.to_csv(counts, index=False)
    sig_fdr['Gene name'].to_csv(gene, index=False)

if __name__ == "__main__":
    main()
