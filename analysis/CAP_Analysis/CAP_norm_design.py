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
cap = samples['Data_id'][samples['Data_id'].str.contains("CAP")].tolist()
hcap = [ _ for _ in cap if _.startswith('H')]
ccap = [ _ for _ in cap if _.startswith('C')]

# Drops samples that's not HDPos
col = list(df)
cols = col[:1] + hcap + ccap
df = df[cols]
df = df[(df[hcap].mean(axis=1) > 5) | (df[ccap].mean(axis=1) > 5)]

df.to_csv(os.path.abspath("../../samples/Analysis_Results/CAP_from_norm.csv"),index=False)
