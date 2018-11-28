import pandas as pd
import os

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
fl = os.path.abspath('../../samples/all_salmon_quant.tsv')

# Read salmon counts and drop H_0014_BA9_mRNASeq sample
df = pd.read_csv(fl, sep="\t")
#df = df.drop('H_0014_BA9_mRNASeq', axis=1)

# Read sample info file
sample_info = pd.read_csv(fn, comment='#')
# column15: datasetid, 2: subject type, 5: death age
samples = sample_info.iloc[:,[14,1,4]]
# Rename columns
samples.columns = ["Data_id", "Subject_type", "Subject_death"]
# Drop H_0014_BA9_mRNASeq sample
#samples = samples[samples['Data_id']!='H_0014_BA9_mRNASeq']
# All sample ids
ids = samples['Data_id'].tolist()

# To get HD or control means
dataset_ids = samples['Data_id'].tolist()
control_ids = [ _ for _ in dataset_ids if _.startswith('C')]
HD_ids = [ _ for _ in dataset_ids if _.startswith('H')]

##################### Filtering #####################
df['avg_control'] = df[control_ids].mean(axis=1)
df['avg_HD'] = df[HD_ids].mean(axis=1)
df = df[(df.avg_control > 5) | (df.avg_HD > 5)]
df = df.drop('avg_control', axis=1)
df = df.drop('avg_HD', axis=1)
df.to_csv(os.path.abspath("../../samples/Analysis_Results/all_filter_w_H_0014.csv"), index=False)

##################### For sample_info design ######################
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.drop(0)
sample_i = sample_i.rename(columns = {0:'Data_id'})
df_new = pd.merge(sample_i,samples, on='Data_id')
df_new.to_csv(os.path.abspath("../../samples/Analysis_Results/all_info_design_H0014.csv"), index=False)

