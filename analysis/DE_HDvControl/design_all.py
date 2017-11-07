import pandas as pd
import os

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')
fl = os.path.abspath('../../samples/all_salmon_quant_rrna.tsv')

# Read salmon counts and drop H_0014_BA9_mRNASeq sample
df = pd.read_csv(fl, sep="\t")
df = df.drop('H_0014_BA9_mRNASeq', axis=1)

# Read sample info file
sample_info = pd.read_csv(fn, comment='#')
# column15: datasetid, 2: subject type, 5: death age
samples = sample_info.iloc[:,[14,1,4]]
# Rename columns
samples.columns = ["Data_id", "Subject_type", "Subject_death"]
# Drop H_0014_BA9_mRNASeq sample
samples = samples[samples['Data_id']!='H_0014_BA9_mRNASeq']

##################### For sample_info design ######################
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.drop(0)
sample_i = sample_i.rename(columns = {0:'Data_id'})
df_new = pd.merge(sample_i,samples, on='Data_id')
df_new.to_csv("all_info_design.csv", index=False)

