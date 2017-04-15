import pandas as pd
import os

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')

col_names = ['Dataset.dataset_id','Dataset.protocol']

files = pd.read_csv(fn, sep=',', comment='#', usecols=col_names)

#files.to_csv('protocol.tsv', sep='\t')

files.reindex(columns=['Dataset.protocol','Dataset.dataset_id']).to_csv('sample_info.csv', sep=',', index=False,header=False)
