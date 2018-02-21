import pandas as pd
import os

#Last edit: 02/21/2018

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
fl = os.path.abspath('../../samples/all_salmon_quant.tsv')

# Read sample info
df = pd.read_csv(fl, sep="\t")
sample_info = pd.read_csv(fn, sep=",", comment='#')

# column15: datasetid, 2: subject type, 5: death age
samples = sample_info.iloc[:,[14,1,4]].copy()
samples.columns = ["Data_id","Subject_type","Subject_death"]
# Get HDPos sample IDs only
hdpos = samples['Data_id'][samples['Subject_type'].str.contains("HDpos")].tolist()

###################### Filtering ###########################
BA9 = [ _ for _ in hdpos if "BA9" in _]
CAP = [ _ for _ in hdpos if "CAP" in _]

# Drops samples that's not HDPos
col = list(df)
cols = col[:1] + hdpos
df.drop([col for col in df.columns if col not in cols],axis=1,inplace=True)

# Drops rows if mean of CAP && BA9 < 5
df = df[(df[BA9].mean(axis=1) > 5) | (df[CAP].mean(axis=1) > 5)]

df.to_csv(os.path.abspath("../../samples/Analysis_Results/asymp_filter.csv"),index=False)

###################### Creating info design file ###########################
sample_i = pd.DataFrame(hdpos)
sample_i = sample_i.rename(columns = {0:'Data_id'})
sample_i = pd.merge(sample_i,samples, on='Data_id')
sample_i.loc[:,'Subject_type'] = sample_i["Data_id"].map(lambda x:"BA9" if "BA9" in x else "CAP")

sample_i.to_csv(os.path.abspath("../../samples/Analysis_Results/asymp_info_design.csv"),index=False)

