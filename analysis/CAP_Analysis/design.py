import pandas as pd
import os

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')

fl = os.path.abspath('../../samples/all_salmon_quant_rrna.tsv')

# Read sample info
samples = pd.read_csv(fn, sep=",", comment='#')
# column15: datasetid, 2: subject type, 5: death age
sample_info = samples.iloc[:,[14,1,4]]
# Pulls only CAP samples
sample_info = sample_info[sample_info['Dataset.dataset_id'].str.contains("CAP")]
# To get HD or control
control_ids = [ _ for _ in sample_info['Dataset.dataset_id'] if _.startswith('C')]
HD_ids = [ _ for _ in sample_info['Dataset.dataset_id'] if _.startswith('H')]

# Pulls only CAP samples from counts file
df = pd.read_csv(fl, sep='\t', comment='#')
df.drop([col for col in df.columns if 'BA9' in col],axis=1,inplace=True)

# Filtering, dropping any rows with a 0
df = df[(df != 0).all(1)]

# Filtering, drop those with control mean < 5 and HD mean <5
df['avg_control'] = df[control_ids].mean(axis=1)
df['avg_HD'] = df[HD_ids].mean(axis=1)
df = df[(df.avg_control > 5) & (df.avg_HD > 5)]
df = df.drop('avg_control', axis=1)
df = df.drop('avg_HD', axis=1)
# Creates new file with only CAP samples
df.to_csv("CAP_salmon_filter.csv", index=False)

# For sample_info design
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.drop(0)
sample_i = sample_i.rename(columns = {0:'Dataset.dataset_id'})
df_new = pd.merge(sample_i,sample_info, on='Dataset.dataset_id')
df_new.to_csv("CAP_info_design.csv", index=False)
