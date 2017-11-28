import pandas as pd
import os
RSE_path = os.path.abspath('../../samples/GTEx/results/RSE_deseq2_results.csv')
GTEX_path = os.path.abspath('../../samples/GTEx/results/GTEX_deseq2_results.csv')
cols = ['log2FoldChange', 'padj']
RSE = pd.read_csv(RSE_path, sep=',', index_col=0, usecols=cols)
GTEX = pd.read_csv(GTEX_path, sep=',',index_col=0, usecols=cols)
#new_RSE = RSE[RSE['log2FoldChange'] > 1]
print(RSE.head)
#new_RSE = RSE.sort_values(['log2FoldChange'])
#print(new_RSE.head)
#print(new_RSE.head)
#print('greater than 2 RSE', RSE[RSE['log2FoldChange'] > 2.0] | RSE[RSE['log2FoldChange']> -2.0].shape)
#print()
#print('greater than 2 GTEX', GTEX[GTEX['log2FoldChange'] > 2.0] | GTEX[GTEX['log2FoldChange']> -2.0].shape)

