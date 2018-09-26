import pandas as pd
import os

# Last edit: 07/24/2018

fl = os.path.abspath('../HD_mRNASeq_sample_info.csv')
fq = os.path.abspath('../../samples/all_salmon_quant.tsv')
fn = os.path.abspath('../../samples/Analysis_Results/all_norm.csv')

# Read sample info
df = pd.read_csv(fq, sep="\t")
dn = pd.read_csv(fn, sep=",", comment='#')
sample_info = pd.read_csv(fl, sep=",", comment='#')

# column15: datasetid, 2: subject type, 5: death age
samples = sample_info.iloc[:,[14,1,4]].copy()
samples.columns = ["Data_id","Subject_type","Subject_death"]
# Get HDPos sample IDs only
hdpos = samples['Data_id'][samples['Subject_type'].str.contains("HDpos")].tolist()

###################### Filtering ###########################
BA9 = [ _ for _ in hdpos if "BA9" in _]
CAP = [ _ for _ in hdpos if "CAP" in _]

# Drops samples that's not HDPos
cols = list(df)[:1] + BA9 + CAP
df = df[cols]

# Sep 24 edit: Drops rows if mean of CAP | BA9 < 10
df = df[(df[BA9].mean(axis=1) > 10) | (df[CAP].mean(axis=1) > 10)]

# Jul 5 edit: Drop row with all zeros
df = df[df.sum(axis=1) != 0]
# Jul 24 edit: Drop if BA9 group have more than 1 zeros
df = df[(df[BA9] == 0).astype(int).sum(axis=1) != 2]

df.to_csv(os.path.abspath("../../samples/Analysis_Results/asymp_filter.csv"),index=False)

###################### Counts from norm  ###########################
dn = dn[cols]
#dn = dn[(dn[BA9].mean(axis=1) > 5) | (dn[CAP].mean(axis=1) > 5)]

dn = dn[dn.sum(axis=1) != 0]
dn = dn[(dn[BA9] == 0).astype(int).sum(axis=1) != 2]

dn.to_csv(os.path.abspath("../../samples/Analysis_Results/asymp_from_norm.csv"),index=False)
###################### Creating info design file ###########################
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.rename(columns = {0:'Data_id'})
sample_i = pd.merge(sample_i,samples, on='Data_id')
sample_i.loc[:,'Subject_type'] = sample_i["Data_id"].map(lambda x:"BA9" if "BA9" in x else "CAP")

sample_i.to_csv(os.path.abspath("../../samples/Analysis_Results/asymp_info_design.csv"),index=False)

