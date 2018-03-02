import pandas as pd
import os

#Last edit: 02/21/2018

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
fl = os.path.abspath('../../samples/Analysis_Results/all_norm.csv')

# Read sample info
df = pd.read_csv(fl)
sample_info = pd.read_csv(fn, sep=",", comment='#')

# column15: datasetid, 2: subject type, 5: death age
samples = sample_info.iloc[:,[14,1,4]].copy()
samples.columns = ["Data_id","Subject_type","Subject_death"]
# Get HDPos sample IDs only
hdpos = samples['Data_id'][samples['Subject_type'].str.contains("HDpos")].tolist()
hdpos_BA9 = [_ for _ in hdpos if "BA9" in _]
hdpos_CAP = [_ for _ in hdpos if "CAP" in _]

# Drops samples that's not HDPos
col = list(df)
cols = col[:1] + hdpos_BA9 + hdpos_CAP
df = df[cols]
df = df[(df[hdpos_BA9].mean(axis=1) > 5) | (df[hdpos_CAP].mean(axis=1) > 5)]

df.to_csv(os.path.abspath("../../samples/Analysis_Results/asymp_from_norm.csv"),index=False)

