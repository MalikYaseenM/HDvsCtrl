import pandas as pd
import os

#Last edit: 04/16/2018
# Samples to drop: H_0014_BA9_mRNASeq (outlier) & all ribo-depleted
fl = os.path.abspath('../../HD_mRNASeq_sample_info.csv')
fq = os.path.abspath('../../../samples/all_salmon_quant.tsv')
fn = os.path.abspath('../../../samples/Analysis_Results/all_norm.csv')

# Read sample info
df = pd.read_csv(fq, sep="\t")
sample_info = pd.read_csv(fl, sep=",", comment='#')
dn = pd.read_csv(fn, sep=",", comment='#')

sample_info = sample_info[(sample_info['Dataset.protocol'] != 'TruSeq Ribo-depleted') & (sample_info['Dataset.dataset_id'] != 'H_0014_BA9_mRNASeq')]

# column15: datasetid, 2: subject type, 5: death age
samples = sample_info.iloc[:,[14,1,4]].copy()
samples.columns = ["Data_id","Subject_type","Subject_death"]

# Get BA9 samples only
samples = samples[samples['Data_id'].str.contains('BA9')]
# Get Control and HDPos only
samples = samples[samples['Subject_type'].str.contains('HDpos') | samples['Subject_type'].str.contains('Control')]

# Get ages of HDpos (3 samples)
age = samples.groupby(['Subject_type']).get_group('HDpos')['Subject_death'].tolist()
# Get approximate ages +- 3 years
age = age + [x+y for x in age for y in [1,2,3]] + [x-y for x in age for y in [1,2,3]]
# Control BA9 samples with approximate ages with HDPos BA9 samples
samples = samples[samples['Subject_death'].isin(age)]

# Used to get average counts for each groups
control = samples.groupby(['Subject_type']).get_group('Control')['Data_id'].tolist()
HD_pos = samples.groupby(['Subject_type']).get_group('HDpos')['Data_id'].tolist()

###################### Filtering ###########################
# Drops other samples not in samples list
cols = list(df)[:1] + samples['Data_id'].tolist()
df.drop([col for col in df.columns if col not in cols],axis=1, inplace=True)

# Drops rows if mean of HD && HDpos  < 5
df = df[(df[control].mean(axis=1) > 5) | (df[HD_pos].mean(axis=1) > 5)]
df.to_csv(os.path.abspath("../../../samples/Analysis_Results/as_symp_C_filter.csv"),index=False)

###################### Taking normalized counts from all_norm.csv file ###########################
# Drop other samples not in samples list
dn.drop([col for col in dn.columns if col not in cols],axis=1, inplace=True)
# Drops rows if mean HD && HDPos < 5 [filtered again because there might be 0s in the normalized metadata]
dn = dn[(dn[control].mean(axis=1) > 5) | (dn[HD_pos].mean(axis=1) > 5)]

dn.to_csv(os.path.abspath("../../../samples/Analysis_Results/as_symp_C_from_norm.csv"), index=False)

###################### Creating info design file ###########################
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.rename(columns={0:'Data_id'}).drop(sample_i.index[0])
df_new = pd.merge(sample_i, samples, on='Data_id')

df_new.to_csv(os.path.abspath("../../../samples/Analysis_Results/as_symp_C_info_design.csv"),index=False)
