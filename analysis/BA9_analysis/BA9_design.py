import pandas as pd
import os

# Last edit: 07/5/2018
# Samples to drop: H_0014_BA9_mRNASeq (outlier) & all ribo-depleted

# HD mRNASeq sample info
fl = os.path.abspath('../HD_mRNASeq_sample_info.csv')
# Counts matrix
fq = os.path.abspath('../../samples/all_salmon_quant.tsv')
fn = os.path.abspath('../../samples/Analysis_Results/all_norm.csv')

# Read sample info
sample_info = pd.read_csv(fl, sep=",", comment='#')
# Drop samples
sample_info = sample_info[(sample_info['Dataset.protocol'] != 'TruSeq Ribo-depleted') & (sample_info['Dataset.dataset_id'] != 'H_0014_BA9_mRNASeq')]

# column15: datasetid, 2: subject type, 5: death age
samples = sample_info.iloc[:,[14,1,4]]
# Change column names from . to _
samples.columns = ["Data_id","Subject_type","Subject_death"]
# Drop CAP and asymptomatic samples
samples = samples[samples['Data_id'].str.contains('BA9') & (samples['Subject_type'] != 'HDpos')]

# To get HD or control
control_ids = [ _ for _ in samples['Data_id'] if _.startswith('C')]
HD_ids = [ _ for _ in samples['Data_id'] if _.startswith('H')]

# Pulls only symptomatic BA9 samples from counts file
df = pd.read_csv(fq, sep='\t', comment='#')
cols = list(df)[:1] + control_ids + HD_ids
df = df[cols]

# Sep 20 edit: Filter base mean by 10
# Filtering, drop those with control mean < 10 or HD mean < 10
df = df[(df[control_ids].mean(axis=1) > 10) | (df[HD_ids].mean(axis=1) > 10)]

# Jul5 edit: Drop row with all 0s
df = df[df.sum(axis=1) != 0]
# Jul 24 edit: Drop rows with > 50% zeros on each group
df = df[((df[control_ids] == 0).astype(int).sum(axis=1) <= len(control_ids)/2) | ((df[HD_ids] == 0).astype(int).sum(axis=1) <= len(HD_ids)/2)]

# Creates new file with only BA9 samples
df.to_csv(os.path.abspath("../../samples/Analysis_Results/BA9_filter.csv"), index=False)

###################### Counts from norm  ###########################
dn = pd.read_csv(fn, sep=",", comment='#')
dn = dn[cols]
#dn = dn[(dn[control_ids].mean(axis=1) > 5) | (dn[HD_ids].mean(axis=1) > 5)]

dn = dn[dn.sum(axis=1) != 0]
dn = dn[((dn[control_ids] == 0).astype(int).sum(axis=1) <= len(control_ids)/2) | ((dn[HD_ids] == 0).astype(int).sum(axis=1) <= len(HD_ids)/2)]

dn.to_csv(os.path.abspath("../../samples/Analysis_Results/BA9_from_norm.csv"),index=False)

# For sample_info design
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.rename(columns={0:'Data_id'}).drop(sample_i.index[0])
df_new = pd.merge(sample_i, samples , on='Data_id')

df_new.to_csv(os.path.abspath("../../samples/Analysis_Results/BA9_info_design.csv"), index=False)
