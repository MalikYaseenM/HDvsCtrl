import pandas
import os

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
col_names = ['Dataset.dataset_id','Subject.subject_type','Subject.death']

sample_info = pandas.read_csv('../HD_mRNASeq_sample_info.csv', comment='#')

files = pandas.read_csv(fn, sep=',', comment='#', usecols=col_names)

files.reindex(columns=['Dataset.dataset_id','Subject.subject_type','Subject.death']).to_csv('sample_info_py.csv',sep=',', index=False,header=False)

#files.reindex(columns=['Dataset.dataset_id']).to_csv('sample_info_pytest.csv',sep=',', index=False,header=False)

#files.reindex(columns=['Subject.subject_type','Subject.death']).to_csv('sample_design_py.csv',sep='~', index=False,header=False)

#fn = os.path.abspath('sample_info.csv')


