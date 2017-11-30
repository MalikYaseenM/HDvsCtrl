import pandas as pd
import os


#hd =  os.path.abspath('../HD_mRNASeq_sample_info.csv')
#fn = '../../samples/GTEx/all_salmon_quant.tsv'
#fn = '../../samples/GTEx/CAU_BA9_Salmon.tsv'
# import tsv
fn = 'CAU_BA9_Salmon.tsv'
ALL_GTEX = pd.read_csv(fn, sep='\t')
# get ba9 and cau samples and put them in a list
ba9 = ALL_GTEX.filter(regex='^GTEX-.*BA9.*').columns.tolist()
cau = ALL_GTEX.filter(regex='^GTEX-.*CAU.*').columns.tolist()
# Drop the Unnamed column which is just an index
ALL_GTEX = ALL_GTEX.drop('Unnamed: 0', axis=1)

# Take the meas of ba9 and cau samples and create new columns
ALL_GTEX['avg_ba9'] = ALL_GTEX[ba9].mean(axis=1)
ALL_GTEX['avg_cau'] = ALL_GTEX[cau].mean(axis=1)

# number of rows before filtering 58037
print('before filter',ALL_GTEX.shape)
# Filter rows by means 10
ALL_GTEX = ALL_GTEX[(ALL_GTEX.avg_ba9 > 10) & (ALL_GTEX.avg_cau > 10)]
ALL_GTEX = ALL_GTEX.drop('avg_ba9', axis=1)
ALL_GTEX = ALL_GTEX.drop('avg_cau', axis=1)
print(len(ba9))
# number of rows after filtering 21654
print('after filter', ALL_GTEX.shape)
ALL_GTEX.to_csv("GTEX_salmon_filter.csv", index=False)

