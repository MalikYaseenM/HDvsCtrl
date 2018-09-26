import pandas as pd
import os

# Last edit: 07/5/2018
# Jul5 edit: Filter by row if row is all 0s

fl = os.path.abspath('../HD_mRNASeq_sample_info.csv')
fq = os.path.abspath('../../samples/all_salmon_quant.tsv')
fn = os.path.abspath('../../samples/Analysis_Results/all_norm.csv')

# Read files
sample_info = pd.read_csv(fl, sep=",", comment='#')
df = pd.read_csv(fq, sep='\t', comment='#')
dn = pd.read_csv(fn, sep=",", comment='#')

# column15: datasetid, 2: subject type, 5: death age
samples = sample_info.iloc[:,[14,1,4]]
# Change column names from . to _
samples.columns = ["Data_id","Subject_type","Subject_death"]

# Pulls only CAP samples
samples = samples[samples['Data_id'].str.contains("CAP")]

###################### Filtering ###########################
# To get HD or control ids
control_ids = [ _ for _ in samples['Data_id'] if _.startswith('C')]
HD_ids = [ _ for _ in samples['Data_id'] if _.startswith('H')]

# Pulls only CAP samples from counts file
cols = list(df)[:1] + control_ids + HD_ids
df = df[cols]

# Sep 24 edit: Filtering, drop those with control mean < 10 or HD mean < 10
df = df[(df[control_ids].mean(axis=1) > 10) | (df[HD_ids].mean(axis=1) > 10)]

# Filter out row with more than 3 zeros
df = df[(df == 0).astype(int).sum(axis=1) < 3 ]

# Creates new file with only CAP samples
df.to_csv(os.path.abspath("../../samples/Analysis_Results/CAP_filter.csv"), index=False)

###################### Counts from norm  ###########################
dn = dn[cols]
#dn = dn[(dn[control_ids].mean(axis=1) > 5) | (dn[HD_ids].mean(axis=1) > 5)]
dn = dn[(dn == 0).astype(int).sum(axis=1) < 3 ]
dn.to_csv(os.path.abspath("../../samples/Analysis_Results/CAP_from_norm.csv"),index=False)

###################### Creating info design file ###########################
# For sample_info design
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.rename(columns={0:'Data_id'}).drop(sample_i.index[0])
df_new = pd.merge(sample_i,samples, on='Data_id')
df_new.to_csv(os.path.abspath("../../samples/Analysis_Results/CAP_info_design.csv"), index=False)
