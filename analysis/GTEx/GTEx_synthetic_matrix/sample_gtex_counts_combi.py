import pandas as pd, numpy as np, itertools as it, math, os, random, sys
from collections import OrderedDict, defaultdict

#Make new test counts with n1 samples and n2 samples from df counts

def readfile(dataframe):
    df = pd.read_csv(dataframe, sep='\t')
    return df

def random_sample(df, A, B):
    df_cols = list(df.columns)
    # Using only CAU counts
    CAU_id = [ x for x in df_cols if 'CAU' in x ]
    ############ Sampling ############
    AID = random.sample(CAU_id, A)
    BID = random.sample(CAU_id, B)
    #BID = random.sample(list(set(CAU_id) - set(AID)), B)

    ############ Filtering step ############
    df = df[['gene_id'] + AID + BID].copy()
    df = df[((df[AID] == 0).astype(int).sum(axis=1) <= len(AID)/2) & ((df[BID] == 0).astype(int).sum(axis=1) <= len(BID)/2)]
    df = df.reset_index(drop=True)
    
    ############ Dummy sample names ############
    colA = ['A_{}'.format(i) for i in range(1, len(AID)+1)]
    colB = ['B_{}'.format(i) for i in range(1, len(BID)+1)]
    gene = ['Gene_{}'.format(i) for i in range(1, len(df)+1)]

    ############ Replacing names into dummy names ############
    df.columns = ['gene_id'] + colA + colB
    CAU = pd.concat([pd.DataFrame({'Gene':gene}), df], axis=1)
    CAU = CAU[['Gene'] + colA + colB]
    
    output = '../../../samples/Analysis_Results/CAU_synth_matrix/counts_combi' +  str(A) + '_' + str(B) + '.csv'
    CAU.to_csv(os.path.abspath(output), index=False)

    ############ Calculating means and vars into dictionaries ############
    CAU['Mean_A'], CAU['Mean_B'] = CAU[colA].mean(axis=1), CAU[colB].mean(axis=1)
    CAU['Var_A'], CAU['Var_B'] = CAU[colA].var(axis=1), CAU[colB].var(axis=1)
    means_A, means_B = OrderedDict(zip(CAU['Gene'], CAU['Mean_A'])), OrderedDict(zip(CAU['Gene'], CAU['Mean_B']))
    var_A, var_B = OrderedDict(zip(CAU['Gene'], CAU['Var_A'])), OrderedDict(zip(CAU['Gene'], CAU['Var_B']))

    ############ Genes combination ############
    genes = CAU['Gene'].tolist()
    genes_combi = list(it.combinations(genes, 2))
    gA, gB = [x[0] for x in genes_combi], [x[1] for x in genes_combi]

    # New dataframe with all gene combinations
    pdf = pd.DataFrame({'geneA':gA, 'geneB':gB})

    # Compress dataset by changeing dtypes
    pdf.loc[:,'geneA'] = pdf['geneA'].astype('category')
    pdf.loc[:,'geneB'] = pdf['geneB'].astype('category')

    pdf['logfc_pair'] = np.log2(pdf['geneA'].map(means_A)) - np.log2(pdf['geneB'].map(means_B))
    pdf['logvar_pair'] = np.log2(pdf['geneA'].map(var_A)) - np.log2(pdf['geneB'].map(var_B))
    pdf['baseMean_pair'] = (pdf['geneA'].map(means_A)+pdf['geneB'].map(means_B))/2
    
    # Compress dataset by changeing dtypes
    pdf['logfc_pair'] = pdf['logfc_pair'].astype('float32')
    pdf['logvar_pair'] = pdf['logvar_pair'].astype('float32')
    pdf['baseMean_pair'] = pdf['baseMean_pair'].astype('float32')

    #output2 = output.replace('counts','df')
    #pdf.to_csv(os.path.abspath(output.replace('counts','df')), index=False)
    return(pdf)

def find_min(pdf, baseMeanExp, n1, n2):
    d = pd.DataFrame(columns=['geneA', 'geneB', 'logfc', 'log_var', 'logfc_pair', 'logvar_pair', 'baseMean_pair', 'score'])
    baseMean = 10**baseMeanExp
    for logfc in np.arange(-3,3.2,0.2):
        pdf['logfc'] = float(logfc)
        pdf['logfc'] = pdf['logfc'].astype('float32')
        for log_var in np.arange(-2,2.2,0.2):
            pdf['log_var'] = float(log_var)
            pdf['log_var'] = pdf['log_var'].astype('float32')
            pdf['score'] = ((pdf['baseMean_pair']-baseMean)**2 + (pdf['logfc_pair']-pdf['logfc'])**2 + (pdf['logvar_pair']-pdf['log_var'])**2)
            d = pd.concat([d, pd.DataFrame(pdf.loc[pdf['score'].idxmin()]).T], axis=0, sort=False)
    output = '../../../samples/Analysis_Results/CAU_synth_matrix/df_combi' + str(n1) + '_' + str(n2) +'_bm'+ str(baseMeanExp) + '.csv'
    d.to_csv(os.path.abspath(output), index=False)

def main():
    if len(sys.argv) < 5:
        print("you must call program as:  ")
        print("sample_gtex_counts_combi.py <File> <n1> <n2> <basemean>")
        print("")
        sys.exit(1)
    
    fl = sys.argv[1]
    n1 = int(sys.argv[2])
    n2 = int(sys.argv[3])
    bm = int(sys.argv[4])
    df = readfile(fl)
    rs = random_sample(df, n1, n2)
    fm = find_min(rs, bm, n1, n2)

if __name__ == "__main__":
    main()
