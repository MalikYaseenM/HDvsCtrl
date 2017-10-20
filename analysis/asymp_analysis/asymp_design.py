import pandas as pd
import os

# HD mRNASeq sample info
fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
# Counts matrix
fl = os.path.abspath('../../samples/all_salmon_quant_rrna.tsv')

# Read data
samples = pd.read_csv(fn, sep=",", comment='#')
df = pd.read_csv(fl, sep='\t', comment='#')

# Append the 2 batches for asymp (HDpos) datasets, 5 samples total
batch1 = samples[samples['Dataset.batch']=="mRNA_asymHD"]
batch2 = samples[samples['Dataset.batch']=="mRNA_HD3"]
batch = batch1.append(batch2)

# Takes the desired columns only, Dataset id, Subject type and subject death
sample_info = batch.iloc[:,[14,1,4]].copy()
sample_info.columns = ["Data_id","Subject_type","Subject_death"]
# Makes list of BA9 and CAP column ids
batch_id = sample_info['Data_id'].tolist()
batch_BA9 = [s for s in batch_id  if "BA9" in s]
batch_CAP = [s for s in batch_id  if "CAP" in s]

# Changes the Subject.subject_type to BA9 or CAP
sample_info.loc[:,'Subject_type'] = sample_info["Data_id"].map(lambda x: "BA9" if "BA9" in x else "CAP")

######################### Filtering #######################
col = list(df)
cols = col[:1] + batch_id

# Drop rows if column names not in batch_id
df.drop([col for col in df.columns if col not in cols],axis=1,inplace=True)

# RAW FILTERING STEP 1: Drop rows with any 0s
df = df[(df != 0).all(1)]

# Find average
df['avg_BA9'] = df[batch_BA9].mean(axis=1)
df['avg_CAP'] = df[batch_CAP].mean(axis=1)

# RAW FILTERING STEP 2: Drops rows if avg_BA < 5 or avg_CAP < 5
df = df[(df.avg_BA9 > 5) | (df.avg_CAP > 5)]
# Drops columns
df = df.drop('avg_BA9', axis=1)
df = df.drop('avg_CAP', axis=1)

# Creates new file with only CAP samples
df.to_csv("asymp_salmon_filter.csv", index=False)

######################### Merge for design info  #######################
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.drop(0)
sample_i = sample_i.rename(columns = {0:'Data_id'})
idesign = pd.merge(sample_i,sample_info, on='Data_id')
idesign.to_csv("asymp_info_design.csv", index=False)
