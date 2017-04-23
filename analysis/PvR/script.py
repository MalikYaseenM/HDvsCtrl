import pandas as pd
import os

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')

col_names = ['Dataset.dataset_id','Dataset.protocol']

readin = pd.read_csv(fn, sep=',', comment='#', usecols=col_names)
readin.reindex(columns=['Dataset.protocol','Dataset.dataset_id']).to_csv('sample_info.csv', sep=',', index=False,header=False)
