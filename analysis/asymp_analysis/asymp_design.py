import pandas as pd
import os

# HD mRNASeq sample info
fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
# Counts matrix
fl = os.path.abspath('../../samples/all_salmon_quant_rrna.tsv')

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
cols_BA9 = col[:1] + batch_BA9
cols_CAP = col[:1] + batch_CAP

# Makes 2 new dataframes to drop 0s
df_BA9 = salmon.drop([col for col in salmon.columns if col not in cols_BA9],axis=1, inplace=False)
df_CAP = salmon.drop([col for col in salmon.columns if col not in cols_CAP],axis=1, inplace=False)
df_BA9 = df_BA9[(df_BA9 != 0).all(1)]
df_CAP = df_CAP[(df_CAP != 0).all(1)]

# Merge the 2 dataframes by index
df_batch = pd.concat([df_BA9, df_CAP], axis=1)
# Drop duplicate gene_id column
df_batch = df_batch.loc[:,~df_batch.columns.duplicated()]

# Drop rows with na values
df_batch = df_batch.dropna(axis=0, how='any')
# Find average
df_batch['avg_BA9'] = df_batch[batch_BA9].mean(axis=1)
df_batch['avg_CAP'] = df_batch[batch_CAP].mean(axis=1)

# Drops rows if avg_BA < 5 or avg_CAP < 5
df_batch = df_batch[(df_batch.avg_BA9 > 5) | (df_batch.avg_CAP > 5)]
# Drops columns
df_batch = df_batch.drop('avg_BA9', axis=1)
df_batch = df_batch.drop('avg_CAP', axis=1)

# Creates new file with only CAP samples
df_batch.to_csv("asymp_salmon_filter.csv", index=False)

sample_i = pd.DataFrame(df_batch.columns)
sample_i = sample_i.drop(0)
sample_i = sample_i.rename(columns = {0:'Data_id'})
idesign = pd.merge(sample_i,sample_info, on='Data_id')
idesign.to_csv("asymp_info_design.csv", index=False)
