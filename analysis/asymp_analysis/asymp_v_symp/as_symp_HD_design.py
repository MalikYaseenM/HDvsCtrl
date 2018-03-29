import pandas as pd
import os

#Last edit: 03/28/2018

fl = os.path.abspath('../../HD_mRNASeq_sample_info.csv')
fq = os.path.abspath('../../../samples/all_salmon_quant.tsv')
fn = os.path.abspath('../../../samples/Analysis_Results/all_norm.csv')

# Read sample info
df = pd.read_csv(fq, sep="\t")
sample_info = pd.read_csv(fl, sep=",", comment='#')
dn = pd.read_csv(fn, sep=",", comment='#')

# column15: datasetid, 2: subject type, 5: death age
samples = sample_info.iloc[:,[14,1,4]].copy()
samples.columns = ["Data_id","Subject_type","Subject_death"]

# Get HD BA9 samples only
samples = samples[samples['Subject_type'].str.contains('HD') & samples['Data_id'].str.contains('BA9') & (samples['Data_id'] != 'H_0014_BA9_mRNASeq')]
# Get ages of HDpos (3 samples)
age = samples.groupby(['Subject_type']).get_group('HDpos')['Subject_death'].tolist()
# Take approximate ages
age = age + [48, 50, 79, samples.groupby(['Subject_type']).get_group('HD')['Subject_death'].min(), samples.groupby(['Subject_type']).get_group('HD')['Subject_death'].max()]
# HD BA9 samples with approximate ages with HDPos BA9 samples
samples = samples[samples['Subject_death'].isin(age)]

# Used to get average counts for each groups
HD = samples.groupby(['Subject_type']).get_group('HD')['Data_id'].tolist()
HD_pos = [x for x in samples['Data_id'] if x not in HD]

###################### Filtering ###########################
# Drops other samples not in samples list
cols = list(df)[:1] + samples['Data_id'].tolist()
df.drop([col for col in df.columns if col not in cols],axis=1, inplace=True)

# Drops rows if mean of HD && HDpos  < 5
df = df[(df[HD].mean(axis=1) > 5) | (df[HD_pos].mean(axis=1) > 5)]
df.to_csv(os.path.abspath("../../../samples/Analysis_Results/as_symp_HD_filter.csv"),index=False)

###################### Taking normalized counts from all_norm.csv file ###########################
# Drop other samples not in samples list
dn.drop([col for col in dn.columns if col not in cols],axis=1, inplace=True)
# Drops rows if mean HD && HDPos < 5 [filtered again because there might be 0s in the normalized metadata]
dn = dn[(dn[HD].mean(axis=1) > 5) | (dn[HD_pos].mean(axis=1) > 5)]

dn.to_csv(os.path.abspath("../../../samples/Analysis_Results/as_symp_HD_from_norm.csv"), index=False)

###################### Creating info design file ###########################
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.rename(columns={0:'Data_id'}).drop(sample_i.index[0])
df_new = pd.merge(sample_i, samples, on='Data_id')

df_new.to_csv(os.path.abspath("../../../samples/Analysis_Results/as_symp_HD_info_design.csv"),index=False)
