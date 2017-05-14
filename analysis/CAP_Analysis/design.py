import pandas as pd
import os

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')

fl = os.path.abspath('../diff_exp/all_salmon_norm.csv')

# Read sample info
samples = pd.read_csv(fn, sep=",", comment='#')
# column15: datasetid, 2: subject type, 5: death age
sample_info = samples.iloc[:,[14,1,4]]
# Pulls only CAP samples
sample_info = sample_info[sample_info['Dataset.dataset_id'].str.contains("CAP")]

# Pulls only CAP samples from counts file
df = pd.read_csv(fl, sep=',', comment='#')
df.drop([col for col in df.columns if 'BA9' in col],axis=1,inplace=True)
# Creates new file with only CAP samples
df.to_csv("CAP_salmon_norm.csv", index=False)


# For sample_info design
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.drop(0)
sample_i = sample_i.rename(columns = {0:'Dataset.dataset_id'})
df_new = pd.merge(sample_i,sample_info, on='Dataset.dataset_id')
df_new.to_csv("CAP_info_design.csv", index=False)
