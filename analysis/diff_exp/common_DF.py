import pandas
import os

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
col_names = ['Dataset.dataset_id','Subject.subject_type','Subject.death']

sample_info = pandas.read_csv('../HD_mRNASeq_sample_info.csv', comment='#')

files = pandas.read_csv(fn, sep=',', comment='#', usecols=col_names)

# To get sample name, subject type and age of death all in one file
files.reindex(columns=['Dataset.dataset_id','Subject.subject_type','Subject.death']).to_csv('sample_info_py.csv',sep=',', index=False,header=False)

# To get sample name only
#files.reindex(columns=['Dataset.dataset_id']).to_csv('sample_info_pytest.csv',sep=',', index=False,header=False)

# To get subject type and age of death only
#files.reindex(columns=['Subject.subject_type','Subject.death']).to_csv('sample_design_py.csv',sep='~', index=False,header=False)


# Makes new normalized salmon count file with no header
fl = os.path.abspath('../../samples/all_salmon_quant.tsv')

with open(fl,'r') as f:
    with open("normalized_salmon_wo_header.tsv",'w') as f1:
        next(f) # skip header line
        for line in f:
            f1.write(line)


