import pandas as pd, numpy as np, scipy, random, statistics, sys, os
import math, itertools as it
from collections import OrderedDict, defaultdict

def readfile(dataframe):
    # Optimizing dataset to use less memory
    col_types = {'baseMean_pair': 'float32', 'geneA': 'category', 'geneB': 'category', 'logfc_pair': 'float32', 'logvar_pair': 'float32'}
    df = pd.read_csv(dataframe, sep=",", dtype = col_types)
    return df

def max_pairs(name, pdf, bm):
    baseMean = 10**int(bm)
    pairs_df = pd.DataFrame(columns=['geneA', 'geneB','logfc', 'log_var', 'logfc_pair', 'logvar_pair','baseMean_pair','score'])
    for logfc in np.arange(-3,3.2,0.2):
        pdf['logfc'] = float(logfc)
        pdf['logfc'] = pdf['logfc'].astype('float32')
        for log_var in np.arange(-2,2.2,0.2):
            pdf['log_var'] = float(log_var)
            pdf['log_var'] = pdf['log_var'].astype('float32')
            pdf['score'] = ((pdf['baseMean_pair']-baseMean)**2 + (pdf['logfc_pair']-pdf['logfc'])**2 + (pdf['logvar_pair']-pdf['log_var'])**2)
            pairs_df = pd.concat([pairs_df, pd.DataFrame(pdf.iloc[pdf['score'].idxmin()]).T], axis=0, sort=False)
#    if 'df' in name:
        #d = re.findall(r'\d+', name)
        #output = '../../../samples/Analysis_Results/CAU_synth_matrix/df'+str(d[0])+'_'+str(d[1]) + '_pairs_baseExp' + str(bm) +'.csv'
#    output = '../../../samples/Analysis_Results/CAU_synth_matrix/df_combi2_2_bm' + str(bm) +'.csv'
    output = '../../../samples/Analysis_Results/CAU_synth_matrix/' + name + '_bm'+ str(bm) + '.csv'
    pairs_df.to_csv(os.path.abspath(output), index=False)
#    else:
#        output = '../../../samples/Analysis_Results/CAU_synth_matrix/pairs_baseExp' + str(bm) +'.csv'
#    pairs_df.to_csv(os.path.abspath(output), index=False)

def main():
    if len(sys.argv) < 3:
        print("you must call program as:  ")
        print("synth_mat.py <CSV file> <basemean power>")
        print("")
        sys.exit(1)

    fl = sys.argv[1]
    name = str(sys.argv[1])
#    name = output.split('/')[-1].split('.')[0]
    name = name.split('/')[-1].split('.')[0]
    bm = int(sys.argv[2])
    df = readfile(fl)
    pairs_df = max_pairs(name, df, bm)

if __name__ == "__main__":
    main()
