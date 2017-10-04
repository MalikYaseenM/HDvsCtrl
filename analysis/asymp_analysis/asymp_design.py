import pandas as pd
import os

# HD mRNASeq sample info
fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
# Counts matrix
fl = os.path.abspath('../../samples/all_salmon_quant_rrna.tsv')

# Read sample info
samples = pd.read_csv(fn, sep=",", comment='#')
# column15: datasetid, 2: subject type, 5: death age
sample_info = samples.iloc[:,[14,1,4]]
# Pulls only HDPos (Asymptomatic) samples
asymp = sample_info[sample_info['Subject.subject_type'].str.contains("HDpos")]

# Changes the Subject.subject_type to BA9 or CAP
pd.options.mode.chained_assignment = None 
asymp.loc[:,'Subject.subject_type'] = asymp["Dataset.dataset_id"].map(lambda x: "BA9" if "BA9" in x else "CAP" if "CAP" in x else "")

# To get BA9 or CAP
BA9_ids = asymp['Dataset.dataset_id'][asymp['Dataset.dataset_id'].str.contains("BA9")].tolist()
CAP_ids = asymp['Dataset.dataset_id'][asymp['Dataset.dataset_id'].str.contains("CAP")].tolist()
hdpos = asymp['Dataset.dataset_id'].tolist()

# Pulls only HDPos samples from counts matrix
df = pd.read_csv(fl, sep='\t', comment='#')
df.drop([col for col in df.columns if col not in hdpos],axis=1,inplace=True)

# Filtering, dropping any rows with a 0
df = df[(df != 0).all(1)]

# Filtering, drop those with CAP mean < 5 and BA9 mean <5
df['avg_BA9'] = df[BA9_ids].mean(axis=1)
df['avg_CAP'] = df[CAP_ids].mean(axis=1)
df = df[(df.avg_BA9 > 5) & (df.avg_CAP > 5)]
df = df.drop('avg_BA9', axis=1)
df = df.drop('avg_CAP', axis=1)
# Creates new file with only CAP samples
df.to_csv("asymp_salmon_filter.csv", index=False)

# For sample info design
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.rename(columns = {0:'Dataset.dataset_id'})
df_new = pd.merge(sample_i,asymp, on='Dataset.dataset_id')
df_new.to_csv("asymp_info_design.csv", index=False)
