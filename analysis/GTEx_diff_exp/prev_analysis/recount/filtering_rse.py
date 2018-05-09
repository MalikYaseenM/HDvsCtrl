import pandas as pd
import os
df = pd.read_csv('RSE_ID_Counts.tsv', sep='\t')

ba9 = df.filter(regex='^GTEX-.*BA9.*').columns.tolist()
cau = df.filter(regex='^GTEX-.*CAU.*').columns.tolist()

df['avg_ba9'] = df[ba9].mean(axis=1)
df['avg_cau'] = df[cau].mean(axis=1)

print('before filter', df.shape)
df = df[(df.avg_ba9 > 10) & (df.avg_cau > 10)]
df = df.drop('avg_ba9', axis=1)
df = df.drop('avg_cau', axis=1)
print(len(ba9))
print('after filter', df.shape)
loc = os.path.abspath("../../samples/GTEx/results/filtered_RSE.csv")
df.to_csv(loc, index=False)
print('done')
