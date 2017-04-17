import pandas
import os

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
col_names = ['Dataset.dataset_id']

sample_info = pandas.read_csv('../HD_mRNASeq_sample_info.csv', comment='#')

files = pandas.read_csv(fn, sep=',', comment='#', usecols=col_names)

files.reindex(columns=['Dataset.dataset_id']).to_csv('sample_info.csv', index=False,header=False)
