import pandas as pd
import os

created = 'samplepoly.tsv'
fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
#col_names = ['Dataset.dataset_id']
col_names = ['Dataset.dataset_id','Dataset.protocol']
#files = pd.read_csv(fn, sep=',', comment='#', names=col_names)
files = pd.read_csv(fn, sep=',', comment='#', usecols=col_names)
files.to_csv('protocol.tsv', sep='\t')
#print(files)

# with open(created, 'w')as o: 
#     with open(files, 'r') as f:
#         for rows in f:
#             if 'poly-A' in rows:
#                 listpoly.append(rows)
#         print(listpoly[0])
        


#Poly_A = [ _ for _ in files.protocol if _.startswith('poly-A')]

#Ribo_dep = [ _ for _ in dataset_ids if _.startswith('Ribo-depleted')]

