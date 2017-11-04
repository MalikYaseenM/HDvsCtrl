import pandas as pd
import os

fn = os.path.abspath('../HD_mRNASeq_sample_info.csv')

fl = os.path.abspath('../../samples/all_salmon_quant_rrna.tsv')

#df = pd.read_csv("all_salmon_quant_rrna.tsv", sep="\t")
df = pd.read_csv(fl, sep="\t")

sample_info = pd.read_csv(fn, comment='#')
# column15: datasetid, 2: subject type, 5: death age
samples = sample_info.iloc[:,[14,1,4]]

# For sample_info design
sample_i = pd.DataFrame(df.columns)
sample_i = sample_i.drop(0)
sample_i = sample_i.rename(columns = {0:'Dataset.dataset_id'})
df_new = pd.merge(sample_i,samples, on='Dataset.dataset_id')
df_new.to_csv("sample_info_design.csv", index=False)
