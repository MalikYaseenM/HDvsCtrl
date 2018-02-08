import pandas as pd
import os
direct = os.path.abspath('../../samples/GTEx/results/')
os.chdir(direct)

RSE = 'filtered_RSE.csv'
GTEX = 'GTEX_salmon_filter.csv'
df_RSE = pd.read_csv(RSE, sep=',')
df_GTEX = pd.read_csv(GTEX, sep=',')

new_df_RSE = df_RSE.loc[df_GTEX.index]
new_df_RSE.to_csv('subset_RSE.csv', sep=',', index=False)


