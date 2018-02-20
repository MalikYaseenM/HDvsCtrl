import pandas as pd
import os

# HD mRNASeq sample info
fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
# Counts matrix
fl = os.path.abspath('../../samples/all_salmon_quant.tsv')

# Read sample info
samples = pd.read_csv(fn, sep=",", comment='#')
# column15: datasetid, 2: subject type, 5: death age
sample_info = samples.iloc[:,[14,1,4]]
# Pulls only BA9 samples
BA9 = sample_info[sample_info['Dataset.dataset_id'].str.contains("BA9")].copy()
# Change column names from . to _
BA9.columns = ["Data_id","Subject_type","Subject_death"]

# To get HD or control
control_ids = [ _ for _ in BA9['Data_id'] if _.startswith('C')]
HD_ids = [ _ for _ in BA9['Data_id'] if _.startswith('H')]

# Pulls only BA9 samples from counts file
df = pd.read_csv(fl, sep='\t', comment='#')
df.drop([col for col in df.columns if 'CAP' in col],axis=1,inplace=True)
# Filtering
df['avg_control'] = df[control_ids].mean(axis=1)
df['avg_HD'] = df[HD_ids].mean(axis=1)
df = df[(df.avg_control > 5) & (df.avg_HD > 5)]
df = df.drop('avg_control', axis=1)
df = df.drop('avg_HD', axis=1)
# Creates new file with only BA9 samples
df.to_csv(os.path.abspath("../../samples/BA9_salmon_filter.csv"), index=False)

# For sample_info design
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.drop(0)
sample_i = sample_i.rename(columns = {0:'Data_id'})
df_new = pd.merge(sample_i,BA9, on='Data_id')
df_new.to_csv(os.path.abspath("../../samples/BA9_info_design.csv"), index=False)
