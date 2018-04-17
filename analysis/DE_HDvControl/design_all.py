import pandas as pd
import os

# Last edit: April 16 2018
# Samples to drop: H_0014_BA9_mRNASeq (outlier) & all ribo-depleted

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
fl = os.path.abspath('../../samples/all_salmon_quant.tsv')

# Read salmon counts
df = pd.read_csv(fl, sep="\t")

# Read sample info file
sample_info = pd.read_csv(fn, comment='#')
sample_info = sample_info[(sample_info['Dataset.protocol'] != 'TruSeq Ribo-depleted') & (sample_info['Dataset.dataset_id'] != 'H_0014_BA9_mRNASeq')]

# column15: datasetid, 2: subject type, 5: death age
samples = sample_info.iloc[:,[14,1,4]]
# Rename columns
samples.columns = ["Data_id", "Subject_type", "Subject_death"]

# To get HD or control means
dataset_ids = samples['Data_id'].tolist()
control_ids = [ _ for _ in dataset_ids if _.startswith('C')]
HD_ids = [ _ for _ in dataset_ids if _.startswith('H')]

##################### Filtering #####################
# Drops unwanted samples (i.e. H_0014 & ribo depleted)
df = df[list(df)[:1] + dataset_ids]
# Drops group count with low means
df = df[(df[control_ids].mean(axis=1) > 5) | (df[HD_ids].mean(axis=1) > 5)]

df.to_csv(os.path.abspath("../../samples/Analysis_Results/all_filter.csv"), index=False)

##################### For sample_info design ######################
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.rename(columns={0:'Data_id'}).drop(sample_i.index[0])
df_new = pd.merge(sample_i,samples, on='Data_id')
df_new.to_csv(os.path.abspath("../../samples/Analysis_Results/all_info_design.csv"), index=False)

