import pandas as pd, os
# Import the file with batch data
dir_path = os.path.abspath('../../samples/GTEx/results/phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt')

df = pd.read_csv(dir_path, sep='\t', comment='#')

#only get the BA9 and CAU samples
df = df.loc[df['SMTSD'].isin(['Brain - Frontal Cortex (BA9)', 'Brain - Caudate (basal ganglia)'])]
#SAMPID full name, SMTSD info on BA9 or CAU, SMNABTCH batch info
df = df[['SAMPID','SMTSD', 'SMNABTCH']]

# create a key column which I can merge our info file on
columns = df['SAMPID'].tolist()
BA9 =  ['R10A', 'R10B', 'R10a', 'R10b']
CAU =  ['R5A', 'R5B', 'R5a', 'R5b']
key_col = []

for cols in columns:
    split = cols.split('-')
    if split[3] in BA9:
        key_col.append(split[0]+'-'+split[1]+'_BA9_mRNASeq')        
    elif split[3] in CAU:
        key_col.append(split[0]+'-'+split[1]+'_CAU_mRNASeq')
df['KEY'] = key_col


our_samples = os.path.abspath('../../samples/GTEx/results/new_GTEX_info.tsv')
our_df = pd.read_csv(our_samples, sep='\t')

merge_df = df.merge(our_df, left_on='KEY', right_on='identifier', how='inner')
final_df = merge_df[["identifier", "SMNABTCH", "brain_region", "sample_name", "SAMPID", "Run_s", "abspath"]]
final_df.to_csv(os.path.abspath('../../samples/GTEx/results/GTEx_info.tsv'), sep='\t', index=False)
os.rename(os.path.abspath('../../samples/GTEx/results/new_GTEX_info.tsv'), os.path.abspath('../../samples/GTEx/results/OLD_GTEX_info.tsv'))
